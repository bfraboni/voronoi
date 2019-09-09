
#ifndef _APP_CAMERA_H
#define _APP_CAMERA_H

#include "app.h"
#include "orbiter.h"


//! \addtogroup application utilitaires pour creer une application, avec un orbiter manipulable a la souris

//! \file
/*! squelette d'application: creation d'une fenetre, d'un contexte openGL et gestion des evenements.
    tuto7.cpp et tuto8.cpp presentent un exemple simple d'utilisation.

    la class App expose les fonctionnalites de window.h, elles sont juste presentees differemment.
    les fonctions globales de window.h sont toujours utilisables (a part run() qui est remplace par App::run()).
*/

//! classe application.
class AppCamera : public App
{
public:
    //! constructeur, dimensions de la fenetre et version d'openGL.
    AppCamera( const int width, const int height, const int major= 3, const int minor= 3 );
    virtual ~AppCamera( );

    //! a deriver pour creer les objets openGL. renvoie -1 pour indiquer une erreur, 0 sinon.
    virtual int init( ) = 0;
    //! a deriver pour detruire les objets openGL. renvoie -1 pour indiquer une erreur, 0 sinon.
    virtual int quit( ) = 0;

    //! a deriver et redefinir pour animer les objets en fonction du temps. 
    virtual int update_camera( const float time, const float delta ) { return 0; };
    //! a deriver pour afficher les objets. renvoie 1 pour continuer, 0 pour fermer l'application.
    virtual int render( ) = 0;
    
    const Orbiter& camera( ) const { return m_camera; }
    Orbiter& camera( ) { return m_camera; }

protected:
    virtual int update( const float time, const float delta );

    Orbiter m_camera;
};


#endif // _APP_CAMERA_H
