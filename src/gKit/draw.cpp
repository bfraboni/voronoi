
#include "draw.h"
#include "window.h"


void draw( Mesh& m, const Transform& model, const Transform& view, const Transform& projection, const GLuint texture )
{
    m.draw(model, view, projection, /* use light */ false, Point(), White(), /* use texture */ true, texture, /* use alpha test */ false, 0);
}

void draw( Mesh& m, const Transform& model, const Transform& view, const Transform& projection )
{
    m.draw(model, view, projection, /* use light */ false, Point(), White(), /* use texture */ false, 0, /* use alpha test */ false, 0);
}

void draw( Mesh& m, const Transform& model, const Orbiter& camera, const GLuint texture )
{
    // recupere les transformations
    Transform view= camera.view();
    Transform projection= camera.projection((float) window_width(), (float) window_height(), 45);
    
    // affiche l'objet
    draw(m, model, view, projection, texture);
}

void draw( Mesh& m, const Orbiter& camera, const GLuint texture )
{
    draw(m, Identity(), camera, texture);
}

void draw( Mesh& m, const Orbiter& camera )
{
    draw(m, Identity(), camera, 0);
}

void draw( Mesh& m, const Transform& model, const Orbiter& camera )
{
    draw(m, model, camera, 0);
}

void draw( Mesh& m, const GLuint program, const bool use_position, const bool use_texcoord, const bool use_normal, const bool use_color )
{
    m.draw(program, use_position, use_texcoord, use_normal, use_color);
}

void DrawParam::draw( Mesh& mesh ) const
{
    mesh.draw(m_model, m_view, m_projection, m_use_light, m_light, m_light_color, m_use_texture, m_texture, m_use_alpha_test, m_alpha_min);
}

void draw( Mesh& m, const DrawParam& param )
{
    param.draw(m);
}
