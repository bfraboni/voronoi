
#ifndef _WINDOW_H
#define _WINDOW_H


#include <SDL2/SDL.h>


//! \addtogroup application utilitaires pour creer une application
///@{

//! \file
//! squelette d'application: creation d'une fenetre, d'un contexte openGL et gestion des evenements. cf tuto1.cpp pour un exemple complet.

typedef SDL_Window *Window;

//! creation d'une fenetre pour l'application.
Window create_window( const int width, const int height );
//! destruction de la fenetre.
void release_window( Window w );

typedef SDL_GLContext Context;

//! cree et configure un contexte opengl.
Context create_context( Window window, const int major= 3, const int minor= 2 );
//! detruit le contexte openGL.
void release_context( Context context );

//! renvoie la largeur de la fenetre de l'application.
int window_width( );
//! renvoie la hauteur de la fenetre de l'application.
int window_height( );

//! renvoie l'etat d'une touche du clavier. cf la doc SDL2 pour les codes.
int key_state( const SDL_Keycode key );
//! desactive une touche du clavier.
void clear_key_state( const SDL_Keycode key );

//! renvoie le dernier evenement. touche speciales.
SDL_KeyboardEvent key_event( );
//! desactive l'evenement.
void clear_key_event( );

//! renvoie le dernier evenement. etat des boutons de la souris.
SDL_MouseButtonEvent button_event( );
//! desactive l'evenement.
void clear_button_event( );

//! renvoie le dernier evenement. etat de la molette de la souris / glisser sur le pad.
SDL_MouseWheelEvent wheel_event( );
//! desactive l'evenement.
void clear_wheel_event( );

//! renvoie le dernier evenement. saisie de texte.
SDL_TextInputEvent text_event( );
//! desactive l'evenement.
void clear_text_event( );


//! renvoie le temps ecoule depuis le lancement de l'application, en millisecondes.
float global_time( );
//! renvoie le temps ecoule depuis la derniere frame, en millisecondes.
float delta_time( );

//! fonction principale. gestion des evenements et appel de la fonction draw de l'application.
int run( Window window, int (*draw)( void ) );

int last_event_count( );
bool laptop_mode( );


//! fonction interne de gestion d'evenements.
int events( Window window );

//! renvoie le chemin(path) vers le fichier filename après l'avoir chercher par rapport à l'executable ou au répertoire père de l'executable
const char* smart_path(const char* filename);

///@}
#endif
