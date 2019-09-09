
#ifndef _TEXTURE_H
#define _TEXTURE_H

#include "glcore.h"
#include "image.h"
#include "image_io.h"


//! \addtogroup openGL
///@{

//! \file 
//! texture2D openGL.


//! cree une texture a partir d'une image im. a detruire avec glDeleteTextures( ). 
//! \param texel_type permet de choisir la representation interne des valeurs de la texture.
GLuint make_texture( const int unit, const Image& im, const GLenum texel_type= GL_RGBA32F );

//! cree une texture a partir des donnees d'une image, cf image_io.h. a detruire avec glDeleteTextures( ).
//! \param texel_type permet de choisir la representation interne des valeurs de la texture.
GLuint make_texture( const int unit, const ImageData& im, const GLenum texel_type= GL_RGBA );

//! cree une texture a partir d'un fichier filename. a detruire avec glDeleteTextures( ).
//! \param texel_type permet de choisir la representation interne des valeurs de la texture.
GLuint read_texture( const int unit, const char *filename, const GLenum texel_type= GL_RGBA );

//! renvoie le nombre de mipmap d'une image width x height.
int miplevels( const int width, const int height );

//! enregistre le contenu de la fenetre dans un fichier. doit etre de type .png / .bmp
int screenshot( const char *filename );

//! enregistre le contenu de la fenetre dans un fichier numerote prefixXXX.png. id est le numero de la capture.
int screenshot( const char *prefix, const int id );

/*! capture video. enregistre le contenu de la fenetre dans un fichier prefix%04d.bmp.

pour obtenir une video 30 images par secondes, compresser avec :
avconv -r 30 -f image2 -i prefix%04d.bmp -c:v libx264 -crf 19 video.m4v

verifier que le codec x264 est installe :
avconv -codecs | grep x264

s'il n'est pas installe :
sudo apt-get install libavcodec-extra-53 
ou la version actuelle :
sudo apt-get install libavcodec-extra-*

exemple d'utilisation : cf shader_kit.cpp
\code
    static bool video= false;
    if(key_state(SDLK_RETURN))
    {
        clear_key_state(SDLK_RETURN);
        video= !video;
        
        if(video) printf("start video capture...\n");
        else      printf("stop video capture.\n");
    }
    
    if(video) capture("shader_kit");
\endcode

 */
int capture( const char *prefix );

///@}
#endif
