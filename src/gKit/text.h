
#ifndef _TEXT_H
#define _TEXT_H 

#include "glcore.h"

#include "color.h"


//! \addtogroup application 
///@{

/*! \file
console texte minimaliste de taille fixe : 24 lignes et 128 colonnes.

permet de placer une chaine de caracteres n'importe ou dans la console. utilise la meme convention que printf().

exemple :
\code
#include "text.h"

Text console;

int init( )
{
    console= create_text()
    ...
}

int quit( )
{
    release_text(console);
    ...
}

int draw( ) 
{
    // effacer la console
    clear(console);
    // afficher
    printf(console, 0, 0, "en haut a gauche");
    printf(console, 0, 23, "en bas a gauche");
    printf(console, 0, 10, "au mileu, int= %d string= %s", 10, "coucou");
    printf(console, 0, 12, "et\n avec\n   plusieurs\n    lignes ?")
    
    // dessiner la console
    draw(console, window_width(), window_height());
    ...
}
\endcode    
 */

/*! representation d'une console texte, dimension fixe, 24 lignes de 128 colonnes.

remarque : *oui*, il y a un constructeur et un destructeur par defaut, mais *non*, ils ne creent pas les objets openGL... ils sont crees / detruits par create_text() / release_text()

pourquoi ? le cycle de vie des objets n'est *pas* le meme. pour creer ou detruire un objet openGL, il *faut* qu'un contexte openGL existe... \n
d'ou la creation et desctruction des objets openGL dans les fonctions init() et quit() de l'application, a un moment ou on controle le contexte openGL.
*/
//! \todo interface c++
struct Text
{
    Text( ) : color( White() ), font(0), program(0), vao(0), ubo(0) {}
    
    int buffer[24][128];
    Color color;        //!< couleur du texte.
    GLuint font;        //!< texture contenant les caracteres.
    GLuint program;     //!< shader pour afficher le texte.
    GLuint vao;         //!< vertex array object.
    GLuint ubo;         //!< uniform buffer object, pour transferrer le texte a afficher
};

//! cree une console. a detruire avec release_text( ).
Text create_text( );
//! detruit une console.
void release_text( Text& text );



//! efface le contenu de la console.
void clear( Text& text );
//! affiche un caractere c sur un fond background.
void print_background( Text& text, const int x, const int y, const int background, const char c );
//! affiche un caractere c sur un fond par defaut.
void print_background( Text& text, const int x, const int y, const char *message );
//! affiche un texte a la position x, y. 
void print( Text& text, const int x, const int y, const char *message );
//! affiche un texte a la position x, y sur un fond par defaut.
void printf_background( Text& text, const int x, const int y, const char *format, ... );
//! affiche un texte a la position x, y. meme utilisation que printf().
void printf( Text& text, const int x, const int y, const char *format, ... );

//! choisit une couleur par defaut pour le texte.
void default_color( Text& text, const Color& color );

//! dessine la console.
void draw( const Text& text, const int width, const int height );

///@}
#endif
