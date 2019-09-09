
#ifndef _IMAGE_HDR_H
#define _IMAGE_HDR_H

#include "image.h"


//! \addtogroup image utilitaires pour manipuler des images
//@{

//! \file
//! manipulation directe d'images, format .hdr

//! charge une image a partir d'un fichier .hdr. renvoie Image::error() en cas d'echec. a detruire avec image::release( ).
//! \param filemane nom de l'image .hdr a charger
Image read_image_hdr( const char *filename );

//! enregistre une image dans un fichier .hdr.
int write_image_hdr( const Image& image, const char *filename );

//! renvoie vrai si le nom de fichier se termine par .hdr.
bool is_hdr_image( const char *filename );

//@}

#endif
