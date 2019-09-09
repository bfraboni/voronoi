
#include <cstdio>
#include <string>
#include <algorithm>
#include <cmath>

#ifdef GK_MACOS
#include <SDL2_image/SDL_image.h>
#else
#include <SDL2/SDL_image.h>
#endif

#include "image_io.h"


Image read_image( const char *filename )
{
    // importer le fichier en utilisant SDL_image
    SDL_Surface *surface= IMG_Load(filename);
    if(surface == NULL)
    {
        printf("[error] loading image '%s'... sdl_image failed.\n", filename);
        return Image::error();
    }
    
    // verifier le format, rgb ou rgba
    const SDL_PixelFormat format= *surface->format;
    int width= surface->w;
    int height= surface->h;
    int channels= format.BitsPerPixel / 8;
    
    printf("loading image '%s' %dx%d %d channels...\n", filename, width, height, channels);
    
    Image image(surface->w, surface->h);
    // converti les donnees en pixel rgba, et retourne l'image, origine en bas a gauche.
    if(format.BitsPerPixel == 32)
    {
        int py= 0;
        for(int y= height -1; y >= 0; y--, py++)
        {
            Uint8 *pixel= (Uint8 *) surface->pixels + py * surface->pitch;

            for(int x= 0; x < width; x++)
            {
                Uint8 r= pixel[format.Rshift / 8];
                Uint8 g= pixel[format.Gshift / 8];
                Uint8 b= pixel[format.Bshift / 8];
                Uint8 a= pixel[format.Ashift / 8];

                image(x, y)= Color((float) r / 255.f, (float) g / 255.f, (float) b / 255.f, (float) a / 255.f);
                pixel= pixel + format.BytesPerPixel;
            }
        }
    }

    else
    {
        int py= 0;
        for(int y= height -1; y >= 0; y--, py++)
        {
            Uint8 *pixel= (Uint8 *) surface->pixels + py * surface->pitch;

            for(int x= 0; x < surface->w; x++)
            {
                Uint8 r= 0;
                Uint8 g= 0;
                Uint8 b= 0;
                if(format.BitsPerPixel >=  8) { r= pixel[format.Rshift / 8]; g= r; b= r; }      // rgb= rrr
                if(format.BitsPerPixel >= 16) { g= pixel[format.Gshift / 8]; b= 0; }    // rgb= rg0
                if(format.BitsPerPixel >= 24) { b= pixel[format.Bshift / 8]; }  // rgb

                image(x, y)= Color((float) r / 255.f, (float) g / 255.f, (float) b / 255.f);
                pixel= pixel + format.BytesPerPixel;
            }
        }
    }

    SDL_FreeSurface(surface);
    return image;
}


int write_image( const Image& image, const char *filename )
{
    if(std::string(filename).rfind(".png") == std::string::npos && std::string(filename).rfind(".bmp") == std::string::npos )
    {
        printf("[error] writing color image '%s'... not a .png / .bmp image.\n", filename);
        return -1;
    }

    // flip de l'image : Y inverse entre GL et BMP
    std::vector<Uint8> flip(image.width() * image.height() * 4);

    int p= 0;
    for(int y= 0; y < image.height(); y++)
    for(int x= 0; x < image.width(); x++)
    {
        Color color= image(x, image.height() - y -1);
        Uint8 r= (Uint8) std::min(std::floor(color.r * 255.f), 255.f);
        Uint8 g= (Uint8) std::min(std::floor(color.g * 255.f), 255.f);
        Uint8 b= (Uint8) std::min(std::floor(color.b * 255.f), 255.f);
        Uint8 a= (Uint8) std::min(std::floor(color.a * 255.f), 255.f);

        flip[p]= r;
        flip[p +1]= g;
        flip[p +2]= b;
        flip[p +3]= a;
        p= p + 4;
    }

    SDL_Surface *surface= SDL_CreateRGBSurfaceFrom((void *) &flip.front(), image.width(), image.height(),
        32, image.width() * 4,
#if 0
        0xFF000000,
        0x00FF0000,
        0x0000FF00,
        0x000000FF
#else
        0x000000FF,
        0x0000FF00,
        0x00FF0000,
        0xFF000000
#endif
    );

    int code= -1;
    if(std::string(filename).rfind(".png") != std::string::npos)
        code= IMG_SavePNG(surface, filename);
    else if(std::string(filename).rfind(".bmp") != std::string::npos)
        code= SDL_SaveBMP(surface, filename);

    SDL_FreeSurface(surface);
    if(code < 0)
        printf("[error] writing color image '%s'...\n%s\n", filename, SDL_GetError());
    return code;
}


ImageData read_image_data( const char *filename )
{
    // importer le fichier en utilisant SDL_image
    SDL_Surface *surface= IMG_Load(filename);
    if(surface == NULL)
    {
        printf("[error] loading image '%s'... sdl_image failed.\n%s\n", filename, SDL_GetError());
        return ImageData();
    }

    // verifier le format, rgb ou rgba
    SDL_PixelFormat format= *surface->format;

    int width= surface->w;
    int height= surface->h;
    int channels= format.BitsPerPixel / 8;
    
    if(channels < 3) channels= 3;
    ImageData image(width, height, channels);

    printf("loading image '%s' %dx%d %d channels...\n", filename, width, height, channels);

    // converti les donnees en pixel rgba, et retourne l'image, origine en bas a gauche.
    if(format.BitsPerPixel == 32)
    {
        int py= 0;
        for(int y= height -1; y >= 0; y--, py++)
        {
            Uint8 *pixel= (Uint8 *) surface->pixels + py * surface->pitch;

            for(int x= 0; x < width; x++)
            {
                Uint8 r= pixel[format.Rshift / 8];
                Uint8 g= pixel[format.Gshift / 8];
                Uint8 b= pixel[format.Bshift / 8];
                Uint8 a= pixel[format.Ashift / 8];

                std::size_t offset= image.offset(x, y);
                image.data[offset]= r;
                image.data[offset +1]= g;
                image.data[offset +2]= b;
                image.data[offset +3]= a;
                pixel= pixel + format.BytesPerPixel;
            }
        }
    }

    else 
    {
        int py= 0;
        for(int y= height -1; y >= 0; y--, py++)
        {
            Uint8 *pixel= (Uint8 *) surface->pixels + py * surface->pitch;

            for(int x= 0; x < surface->w; x++)
            {
                Uint8 r= 0;
                Uint8 g= 0;
                Uint8 b= 0;
                
                if(format.BitsPerPixel >=  8) { r= pixel[format.Rshift / 8]; g= r; b= r; }      // rgb= rrr
                if(format.BitsPerPixel >= 16) { g= pixel[format.Gshift / 8]; b= 0; }    // rgb= rg0
                if(format.BitsPerPixel >= 24) { b= pixel[format.Bshift / 8]; }  // rgb

                std::size_t offset= image.offset(x, y);
                image.data[offset]= r;
                image.data[offset +1]= g;
                image.data[offset +2]= b;
                pixel= pixel + format.BytesPerPixel;
            }
        }
    }

    SDL_FreeSurface(surface);
    return image;
}

int write_image_data( ImageData& image, const char *filename )
{
    if(std::string(filename).rfind(".png") == std::string::npos && std::string(filename).rfind(".bmp") == std::string::npos )
    {
        printf("[error] writing color image '%s'... not a .png / .bmp image.\n", filename);
        return -1;
    }

    if(image.size != 1)
    {
        printf("[error] writing color image '%s'... not an 8 bits image.\n", filename);
        return -1;
    }

    // flip de l'image : origine en bas a gauche
    std::vector<Uint8> flip(image.width * image.height * 4);

    int p= 0;
    for(int y= 0; y < image.height; y++)
    for(int x= 0; x < image.width; x++)
    {
        std::size_t offset= image.offset(x, image.height - y -1);
        Uint8 r= image.data[offset];
        Uint8 g= image.data[offset +1];
        Uint8 b= image.data[offset +2];
        Uint8 a= 255;
        if(image.channels > 3)
            a= image.data[offset +3];

        flip[p]= r;
        flip[p +1]= g;
        flip[p +2]= b;
        flip[p +3]= a;
        p= p + 4;
    }

    // construit la surface sdl
    SDL_Surface *surface= SDL_CreateRGBSurfaceFrom((void *) &flip.front(), image.width, image.height,
        32, image.width * 4,
#if 0
        0xFF000000,
        0x00FF0000,
        0x0000FF00,
        0x000000FF
#else
        0x000000FF,
        0x0000FF00,
        0x00FF0000,
        0xFF000000
#endif
    );

    // enregistre le fichier
    int code= -1;
    if(std::string(filename).rfind(".png") != std::string::npos)
        code= IMG_SavePNG(surface, filename);
    else if(std::string(filename).rfind(".bmp") != std::string::npos)
        code= SDL_SaveBMP(surface, filename);

    SDL_FreeSurface(surface);
    if(code < 0)
        printf("[error] writing color image '%s'...\n%s\n", filename, SDL_GetError());
    return code;
}
