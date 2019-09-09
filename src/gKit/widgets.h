
#ifndef _WIDGETS_H
#define _WIDGETS_H

#include "text.h"


//! \addtogroup application 
///@{

/*! \file
interface graphique minimaliste.
les elements sont disposes sur une grille de la droite vers la gauche.
pour changer de ligne, utiliser begin_line( ).

les elements interactifs renvoient vrai lorsqu'ils changent d'etat. 
par exemple, un bouton renvoie vrai lorsque l'on clique dessus. son etat (selectionne ou pas) est stocke dans une variable passee en parametre.

dans l'exemple suivant, des widgets differents sont affiches selon qu'un bouton est selectionne ou pas.
\code
#include "widgets.h"

Widgets widgets;
int click;

int init( )
{
    widgets= create_widgets();
    click= 0;
    ...
}

int quit( )
{
    release_widgets(widgets);
    ...
}

int draw( )
{
    // description des elements de l'interface, a chaque fois qu'il faut l'afficher.
    // pas de stockage, pas de callbacks...
    
    begin(widgets);
        begin_line(widgets);
        label(widgets, "texte");
        label(widgets, "texte formate %d", 10);
        
        button(widgets, "cliquez ici", click);
        if(click)
        {
            begin_line(widgets);
            label(widgets, "ahah tu as clique...");
        }
    end(widgets);
    
    // afficher l'interface
    draw(widgets, window_width(), window_height());
    ...
}
\endcode
*/

//! representation d'une interface graphique minimaliste.
//! \todo interface c++
struct Widgets
{
    Widgets( ) : console(), px(0), py(0), focus(0), fx(0), fy(0), mb(0), mx(0), my(0), wx(0), wy(0), key(0), mod(0) {}
    
    Text console;       //!< affichage des elements de l'interface.
    int px, py;         //!< placement du prochain widget
    
    int focus;          //!< focus
    int fx, fy;         //!< position du focus
    
    int mb;             //!< click
    int mx, my;         //!< position du click
    int wx, wy;         //!< scroll
    
    int key;            //!< touche
    unsigned int mod;   //!< touches supplementaires, alt, ctrl, etc.
};


//! cree une interface graphique. a detruire avec release_widgets( ).
Widgets create_widgets( );
//! detruit l'interface graphique.
void release_widgets( Widgets& widgets );

//! debut de la description des elements de l'interface graphique.
void begin( Widgets& widgets );

//! place les prochains elements sur une nouvelle ligne.
void begin_line( Widgets& widgets );

//! cree un texte. meme fonctionnement que printf().
void label( Widgets& widgets, const char *format, ... );

//! cree un bouton. renvoie true si le bouton a change d'etat.
//! \param text legende du bouton,
//! \param status etat du bouton, 1 selectionne, 0 non selectionne. doit etre initialise a 0 avant la premiere utilisation.
bool button( Widgets& widgets, const char *text, int& status );
//! cree un radio bouton. selectionne une seule option parmis une liste. renvoie true si le bouton a change d'etat.
//! \param text legende du bouton,
//! \param option position du bouton dans le groupe, 
//! \param status option/bouton selectionne. doit etre initialise avec l'option par defaut avant la premiere utilisation.
bool select( Widgets& widgets, const char *text, const int option, int& status );

//! cree une zone de texte scrollable.
//! \param height nombre de lignes.
//! \param text contenu a afficher.
//! \param begin_line premiere ligne du texte a afficher. doit etre initialise avant la premiere utilisation.
void text_area( Widgets& w, const int height, const char *text, int& begin_line );

/*! cree un champ texte editable. renvoie true si le contenu a change.
    \param text_size taille du champ texte et de la chaine text, zero inclus. si text_size fait 10 caracteres, on ne pourra saisir que 9 caracteres.
    \param text contenu du champ, doit etre initialise avant la premiere utilisation. 

\code
char valeur[32]= "";

int draw( )
{
    begin(widgets);
    edit(widgets, sizeof(valeur), valeur);
    end(widgets);
    
    draw(widgets, window_width(), window_height());
}
\endcode
*/
bool edit( Widgets& widgets, const int text_size, char *text );

//! valeur editable par increment.
bool value( Widgets& widgets, const char *label, int& value, const int value_min, const int value_max, const int value_step );
//! valeur editable par increment.
bool value( Widgets& widgets, const char *label, float& value, const float value_min, const float value_max, const float value_step );

//! termine la description des elements de la ligne.
void end_line( Widgets& widgets );

//! termine la description des elements de l'interface graphique.
void end( Widgets& widgets );

//! choisit une couleur par defaut pour le texte.
void default_color( Widgets& widgets, const Color& color );

//! affiche les elements decrits entre begin() et end().
void draw( Widgets& widgets, const int width, const int height );

///@}
#endif
