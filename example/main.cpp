#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <istream>
#include <ostream>
#include <iostream>
#include <algorithm>
#include <random>

static auto const decimation = 2;
static auto const candidates = 96;
static auto const termination = 200;

using namespace std;

struct rgb {float red, green, blue;};
struct img {int width, height; vector<rgb> pixels;};
struct site {float x, y; rgb color;};

img read(string const &name) {
    ifstream file{name, ios::in | ios::binary};
    auto result = img{0, 0, {}};
    if (file.get() != 'P' || file.get() != '6')
        return result;
    auto skip = [&](){
        while (file.peek() < '0' || '9' < file.peek())
            if (file.get() == '#')
                while (file.peek() != '\r' && file.peek() != '\n')
                    file.get();
    };
     auto maximum = 0;
     skip(); file >> result.width;
     skip(); file >> result.height;
     skip(); file >> maximum;
     file.get();
     for (auto pixel = 0; pixel < result.width * result.height; ++pixel) {
         auto red = file.get() * 1.0f / maximum;
         auto green = file.get() * 1.0f / maximum;
         auto blue = file.get() * 1.0f / maximum;
         result.pixels.emplace_back(rgb{red, green, blue});
     }
     return result;
 }

 float evaluate(img const &target, vector<site> &sites) {
     auto counts = vector<int>(sites.size());
     auto variance = vector<rgb>(sites.size());
     for (auto &site : sites)
         site.color = rgb{0.0f, 0.0f, 0.0f};
     for (auto y = 0; y < target.height; y += decimation)
         for (auto x = 0; x < target.width; x += decimation) {
             auto best = 0;
             auto closest = 1.0e30f;
             for (auto index = 0; index < sites.size(); ++index) {
                 float distance = ((x - sites[index].x) * (x - sites[index].x) +
                                   (y - sites[index].y) * (y - sites[index].y));
                 if (distance < closest) {
                     best = index;
                     closest = distance;
                 }
             }
             ++counts[best];
             auto &pixel = target.pixels[y * target.width + x];
             auto &color = sites[best].color;
             rgb delta = {pixel.red - color.red,
                          pixel.green - color.green,
                          pixel.blue - color.blue};
             color.red += delta.red / counts[best];
             color.green += delta.green / counts[best];
             color.blue += delta.blue / counts[best];
             variance[best].red += delta.red * (pixel.red - color.red);
             variance[best].green += delta.green * (pixel.green - color.green);
             variance[best].blue += delta.blue * (pixel.blue - color.blue);
         }
     auto error = 0.0f;
     auto count = 0;
     for (auto index = 0; index < sites.size(); ++index) {
         if (!counts[index]) {
             auto x = min(max(static_cast<int>(sites[index].x), 0), target.width - 1);
             auto y = min(max(static_cast<int>(sites[index].y), 0), target.height - 1);
             sites[index].color = target.pixels[y * target.width + x];
         }
         count += counts[index];
         error += variance[index].red + variance[index].green + variance[index].blue;
     }
     return 10.0f * log10f(count * 3 / error);
 }

 void write(string const &name, int const width, int const height, vector<site> const &sites) {
     ofstream file{name, ios::out};
     file << width << " " << height << endl;
     for (auto const &site : sites)
         file << site.x << " " << site.y << " "
              << site.color.red << " "<< site.color.green << " "<< site.color.blue << endl;
 }

 int main(int argc, char **argv) {
     auto rng = mt19937{random_device{}()};
     auto uniform = uniform_real_distribution<float>{0.0f, 1.0f};
     auto target = read(argv[1]);
     auto sites = vector<site>{};
     for (auto point = atoi(argv[2]); point; --point)
         sites.emplace_back(site{
             target.width * uniform(rng),
             target.height * uniform(rng)});
     auto greatest = 0.0f;
     auto remaining = termination;
     for (auto step = 0; remaining; ++step, --remaining) {
         auto best_candidate = sites;
         auto best_psnr = 0.0f;
         #pragma omp parallel for
         for (auto candidate = 0; candidate < candidates; ++candidate) {
             auto trial = sites;
             #pragma omp critical
             {
                 trial[step % sites.size()].x = target.width * (uniform(rng) * 1.2f - 0.1f);
                 trial[step % sites.size()].y = target.height * (uniform(rng) * 1.2f - 0.1f);
             }
             auto psnr = evaluate(target, trial);
             #pragma omp critical
             if (psnr > best_psnr) {
                 best_candidate = trial;
                 best_psnr = psnr;
             }
         }
         sites = best_candidate;
         if (best_psnr > greatest) {
             greatest = best_psnr;
             remaining = termination;
             write(argv[3], target.width, target.height, sites);
         }
         cout << "Step " << step << "/" << remaining
              << ", PSNR = " << best_psnr << endl;
     }
     return 0;
 }