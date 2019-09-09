
#include <cassert>
#include <cstdio>
#include <algorithm>

#include "texture.h"
#include "image_io.h"


int miplevels( const int width, const int height )
{
    int w= width;
    int h= height;
    int levels= 1;
    while(w > 1 || h > 1)
    {
        w= std::max(1, w / 2);
        h= std::max(1, h / 2);
        levels= levels + 1;
    }
    
    return levels;
}

GLuint make_texture( const int unit, const Image& im, const GLenum texel_type )
{
    if(im == Image::error())
        return 0;

    // cree la texture openGL
    GLuint texture;
    glGenTextures(1, &texture);
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, texture);
    
    // fixe les parametres de filtrage par defaut
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // transfere les donnees dans la texture, 4 float par texel
    glTexImage2D(GL_TEXTURE_2D, 0,
        texel_type, im.width(), im.height(), 0,
        GL_RGBA, GL_FLOAT, im.buffer());
    
    // prefiltre la texture
    glGenerateMipmap(GL_TEXTURE_2D);
    return texture;
}

GLuint make_texture( const int unit, const ImageData& im, const GLenum texel_type )
{
    if(im.data.empty())
        return 0;
    
    // cree la texture openGL
    GLuint texture;
    glGenTextures(1, &texture);
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, texture);
    
    // fixe les parametres de filtrage par defaut
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    GLenum format;
    switch(im.channels)
    {
        case 1: format= GL_RED; break;
        case 2: format= GL_RG; break;
        case 3: format= GL_RGB; break;
        case 4: format= GL_RGBA; break;
        default: format= GL_RGBA; 
    }
    
    GLenum type;
    switch(im.size)
    {
        case 1: type= GL_UNSIGNED_BYTE; break;
        case 4: type= GL_FLOAT; break;
        default: type= GL_UNSIGNED_BYTE;
    }
    
    // transfere les donnees dans la texture
    glTexImage2D(GL_TEXTURE_2D, 0,
        texel_type, im.width, im.height, 0,
        format, type, im.buffer());
    
    // prefiltre la texture
    glGenerateMipmap(GL_TEXTURE_2D);
    return texture;
}


GLuint read_texture( const int unit, const char *filename, const GLenum texel_type )
{
    ImageData image= read_image_data(filename);
    return make_texture(unit, image, texel_type);
}


int screenshot( const char *filename )
{
    // recupere le contenu de la fenetre / framebuffer par defaut
    glFinish();

    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    glReadBuffer(GL_BACK);

    // recupere les dimensions de la fenetre
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    // transfere les pixels
    ImageData image(viewport[2], viewport[3], 4);
    glReadPixels(0, 0, image.width, image.height,
        GL_RGBA, GL_UNSIGNED_BYTE, image.buffer());
    // ecrit l'image
    return write_image_data(image, filename);
}

int screenshot( const char *prefix, const int id )
{
    char tmp[4096];
    sprintf(tmp,"%s%02d.png", prefix, id);
    return screenshot(tmp);
}

int capture( const char *prefix )
{
    static int id= 1;

    char tmp[4096];
    sprintf(tmp,"%s%04d.bmp", prefix, id);

    if(id % 30 == 0)
        printf("capture frame '%s'...\n", tmp);

    id++;
    return screenshot(tmp);
}

