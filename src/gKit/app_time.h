
#ifndef _APP_TIME_H
#define _APP_TIME_H

#include "glcore.h"
#include "app.h"
#include "text.h"


//! \addtogroup application utilitaires pour creer une application

//! \file
//! classe application, avec mesure integree du temps d'execution cpu et gpu.
class AppTime : public App
{
public:
    //! constructeur, dimensions de la fenetre et version d'openGL.
    AppTime( const int width, const int height, const int major= 3, const int minor= 3 );
    virtual ~AppTime( );

    //! a deriver pour creer les objets openGL.
    virtual int init( ) = 0;
    //! a deriver pour detruire les objets openGL.
    virtual int quit( ) = 0;

    //! a deriver et redefinir pour animer les objets en fonction du temps.
    using App::update;

    //! a deriver pour afficher les objets.
    virtual int render( ) = 0;

    //! execution de l'application.
    int run( );

protected:
    Text m_console;
    GLuint m_time_query;
};


#endif // _APP_H
