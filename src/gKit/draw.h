
#ifndef _DRAW_H
#define _DRAW_H

#include "mesh.h"
#include "orbiter.h"
#include "window.h"


//! \addtogroup objet3D
///@{

//! \file
//! dessine un objet du point de vue d'une camera

//! dessine l'objet avec les transformations model, vue et projection.
void draw( Mesh& m, const Transform& model, const Transform& view, const Transform& projection );
//! applique une texture a la surface de l'objet. ne fonctionne que si les coordonnees de textures sont fournies avec tous les sommets de l'objet.
void draw( Mesh& m, const Transform& model, const Transform& view, const Transform& projection, const GLuint texture );

//! dessine l'objet avec les transformations vue et projection, definies par la camera. model est la transformation identite.
void draw( Mesh& m, const Orbiter& camera );
//! dessine l'objet avec une transformation model. les transformations vue et projection sont celles de la camera
void draw( Mesh& m, const Transform& model, const Orbiter& camera );

//! dessine l'objet avec les transformations vue et projection. model est l'identite. applique une texture a la surface de l'objet. ne fonctionne que si les coordonnees de textures sont fournies avec tous les sommets de l'objet.
void draw( Mesh& m, const Orbiter& camera, const GLuint texture );
//! dessine l'objet avec une transformation model. les transformations vue et projection sont celles de la camera. applique une texture a la surface de l'objet. ne fonctionne que si les coordonnees de textures sont fournies avec tous les sommets de l'objet.
void draw( Mesh& m, const Transform& model, const Orbiter& camera, const GLuint texture );

//! dessine l'objet avec un shader program "specifique"
void draw( Mesh& m, const GLuint program, const bool use_position= true, const bool use_texcoord= true, const bool use_normal= true, const bool use_color= true );

/*! representation des options / parametres d'un draw.
    permet de donner tous les parametres d'un draw de maniere flexible.

    exemple :
    \code
    Mesh objet= { ... };

    DrawParam param;
    param.light(Point(0, 20, 0), Red());
    param.camera(orbiter);
    param.draw(objet);
    \endcode

    ou de maniere encore plus compacte :
    \code
    DrawParam().light(Point(0, 20, 0), Red()).model(m).camera(orbiter).draw(objet);
    \endcode
    les parametres peuvent etre decrits dans un ordre quelconque, mais DrawParam::draw() doit etre appele en dernier.
 */
class DrawParam
{
public:
    //! constructeur par defaut.
    DrawParam( ) : m_model(), m_view(), m_projection(),
        m_use_light(false), m_light(), m_light_color(),
        m_use_texture(false), m_texture(0),
        m_use_alpha_test(false), m_alpha_min(0.3f)
    {}

    //! modifie la transformation model utilisee pour afficher l'objet.
    DrawParam& model( const Transform& m ) { m_model= m; return *this; }
    //! modifie la transformation view utilisee pour afficher l'objet.
    DrawParam& view( const Transform& m ) { m_view= m; return *this; }
    //! modifie la transformation projection utilisee pour afficher l'objet.
    DrawParam& projection( const Transform& m ) { m_projection= m; return *this; }

    //! utilise les transformations view et projection definies par une camera.
    DrawParam& camera( const Orbiter& o ) { m_view= o.view(); m_projection= o.projection((float) window_width(), (float) window_height(), 45); return *this; }
    //! utilise les transformations view et projection definies par une camera. parametres explicites de la projection.
    DrawParam& camera( const Orbiter& o, const float width, const float height, const float fov ) { m_view= o.view(); m_projection= o.projection(width, height, fov); return *this; }
    //! eclaire l'objet avec une source ponctuelle, de position p et de couleur c.
    DrawParam& light( const Point& p, const Color& c= White() ) { m_use_light= true; m_light= p; m_light_color=c; return *this; }
    //! plaque une texture a la surface de l'objet.
    DrawParam& texture( const GLuint t ) { m_use_texture= true; m_texture= t; return *this; }

    //! texture semi transparente, si l'alpha du texel est plus petit que alpha_min, le pixel est transparent. desactive aussi les calculs d'eclairage.
    DrawParam& alpha( const float a=0.f ) { m_use_alpha_test= (a>0.f); m_alpha_min= a; return *this; }

    //! Use light: on/off
    DrawParam& lighting(bool use_light=true) { m_use_light=use_light;  return *this; }

    //! dessine l'objet avec l'ensemble des parametres definis.
    void draw( Mesh& mesh ) const;

    //! renvoie la position de la lumière
    const Point& light() const { return m_light; }

protected:
    Transform m_model;
    Transform m_view;
    Transform m_projection;

    bool m_use_light;
    Point m_light;
    Color m_light_color;

    bool m_use_texture;
    GLuint m_texture;

    bool m_use_alpha_test;
    float m_alpha_min;
};

//! dessine l'objet avec l'ensemble des parametres definis.
void draw( Mesh& mesh, const DrawParam& param );

///@}
#endif
