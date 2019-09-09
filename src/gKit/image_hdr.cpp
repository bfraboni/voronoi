
#include <string>

#include "rgbe.h"
#include "image_hdr.h"


bool is_hdr_image( const char *filename )
{
    return (std::string(filename).rfind(".hdr") != std::string::npos);
}


Image read_image_hdr( const char *filename )
{
    FILE *in= fopen(filename, "rb");
    if(in == NULL)
    {
        printf("[error] loading hdr image '%s'...\n", filename);
        return Image::error();
    }

    int width, height;
    rgbe_header_info info;
    if(RGBE_ReadHeader(in, &width, &height, &info) != RGBE_RETURN_SUCCESS)
    {
        fclose(in);
        printf("[error] loading hdr image '%s'...\n", filename);
        return Image::error();
    }

    std::vector<float> data(width*height*3, 0.f);
    if(RGBE_ReadPixels_RLE(in, &data.front(), width, height) != RGBE_RETURN_SUCCESS)
    {
        fclose(in);
        printf("[error] loading hdr image '%s'...\n", filename);
        return Image::error();
    }

    fclose(in);
    
    // 
    printf("loading hdr image '%s' %dx%d...\n", filename, width, height);
    Image image(width, height);
    
    int i= 0;
    for(int y= 0; y < height; y++)
    for(int x= 0; x < width; x++)
    {
        image(x, height - y -1)= Color(data[i], data[i+1], data[i+2]);
        i= i + 3;
    }
    
    return image;
}

int write_image_hdr( const Image& image, const char *filename )
{
    if(image == Image::error())
        return -1;
    
    FILE *out= fopen(filename, "wb");
    if(out == NULL)
    {
        printf("[error] writing hdr image '%s'...\n", filename);
        return -1;
    }

    int width= image.width();
    int height= image.height();
    if(RGBE_WriteHeader(out, width, height, NULL) != RGBE_RETURN_SUCCESS)
    {
        fclose(out);

        printf("[error] writing hdr image '%s'...\n", filename);
        return -1;
    }

    std::vector<float> data(width*height*3, 0.f);
    int i= 0;
    for(int y= 0; y < height; y++)
    for(int x= 0; x < width; x++)
    {
        Color color= image(x, height - y -1);
        data[i]= color.r;
        data[i+1]= color.g;
        data[i+2]= color.b;
        i= i + 3;
    }

    int code= RGBE_WritePixels_RLE(out, &data.front(), width, height);
    fclose(out);

    if(code != RGBE_RETURN_SUCCESS)
    {
        printf("[error] writing hdr image '%s'...\n", filename);
        return -1;
    }

    printf("writing hdr image '%s'...\n", filename);
    return 0;
}
