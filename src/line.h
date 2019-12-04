#ifndef LINE_H
#define LINE_H

#include "image.h"

void BresenhamLine( Image& image, float x1, float y1, float x2, float y2, const Color& color )
{
    const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
    if(steep)
    {
        std::swap(x1, y1);
        std::swap(x2, y2);
    }

    if(x1 > x2)
    {
        std::swap(x1, x2);
        std::swap(y1, y2);
    }

    const float dx = x2 - x1;
    const float dy = fabs(y2 - y1);

    float error = dx / 2.0f;
    const int ystep = (y1 < y2) ? 1 : -1;
    int y = (int)y1;

    const int maxX = (int)x2;

    for(int x=(int)x1; x<maxX; x++)
    {
        if(steep)
        {
            image(y,x) = color;
        }
        else
        {
            image(x,y) = color;
        }

        error -= dy;
        if(error < 0)
        {
            y += ystep;
            error += dx;
        }
    }
}

void WuLine( const Image& image, Image& output, float x0, float y0, float x1, float y1, const Color& color ) 
{
    auto ipart = [](float x) -> int     {return int(std::floor(x));};
    auto round = [](float x) -> float   {return std::round(x);};
    auto fpart = [](float x) -> float   {return x - std::floor(x);};
    auto rfpart = [=](float x) -> float {return 1 - fpart(x);};
 
    const bool steep = abs(y1 - y0) > abs(x1 - x0);
    if (steep) 
    {
        std::swap(x0,y0);
        std::swap(x1,y1);
    }
    
    if (x0 > x1) 
    {
        std::swap(x0,x1);
        std::swap(y0,y1);
    }
 
    const float dx = x1 - x0;
    const float dy = y1 - y0;
    const float gradient = (dx == 0) ? 1 : dy/dx;
 
    int xpx11;
    float intery;
    {
        const float xend = round(x0);
        const float yend = y0 + gradient * (xend - x0);
        const float xgap = rfpart(x0 + 0.5);
        xpx11 = int(xend);
        const int ypx11 = ipart(yend);
        if (steep) 
        {
            output(ypx11, xpx11) = image(ypx11, xpx11) + rfpart(yend) * xgap * (color - image(ypx11, xpx11));
            output(ypx11 + 1, xpx11) = image(ypx11 + 1, xpx11) + fpart(yend) * xgap * (color - image(ypx11+1, xpx11));
        } 
        else 
        {
            output(xpx11, ypx11) = image(xpx11, ypx11) + rfpart(yend) * xgap * (color - image(xpx11, ypx11));
            output(xpx11, ypx11 + 1) = image(xpx11, ypx11 + 1) + fpart(yend) * xgap * (color - image(xpx11, ypx11 + 1));
        }
        intery = yend + gradient;
    }
 
    int xpx12;
    {
        const float xend = round(x1);
        const float yend = y1 + gradient * (xend - x1);
        const float xgap = rfpart(x1 + 0.5);
        xpx12 = int(xend);
        const int ypx12 = ipart(yend);
        if (steep) 
        {
            output(ypx12, xpx12) = image(ypx12, xpx12) + rfpart(yend) * xgap *(color - image(ypx12, xpx12));
            output(ypx12 + 1, xpx12) = image(ypx12 + 1, xpx12) +  fpart(yend) * xgap *(color - image(ypx12 + 1, xpx12));
        } 
        else 
        {
            output(xpx12, ypx12) = image(xpx12, ypx12) + rfpart(yend) * xgap * (color - image(xpx12, ypx12));
            output(xpx12, ypx12 + 1) = image(xpx12, ypx12 + 1) + fpart(yend) * xgap * (color - image(xpx12, ypx12 + 1));
        }
    }
 
    if (steep) 
    {
        for (int x = xpx11 + 1; x < xpx12; x++) 
        {
            output(ipart(intery), x) = image(ipart(intery), x) + rfpart(intery) * (color - image(ipart(intery), x));
            output(ipart(intery) + 1, x) = image(ipart(intery) + 1, x) +  fpart(intery) * (color - image(ipart(intery) + 1, x));
            intery += gradient;
        }
    } 
    else 
    {
        for (int x = xpx11 + 1; x < xpx12; x++) 
        {
            output(x, ipart(intery)) = image(x, ipart(intery)) + rfpart(intery) * (color - image(x, ipart(intery)));
            output(x, ipart(intery) + 1) = image(x, ipart(intery) + 1) + fpart(intery) * (color - image(x, ipart(intery) + 1));
            intery += gradient;
        }
    }
}

#endif