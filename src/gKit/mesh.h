
#ifndef _MESH_H
#define _MESH_H

#include <vector>
#include <unordered_map>

#include "glcore.h"

#include "vec.h"
#include "mat.h"
#include "color.h"


//! \addtogroup objet3D utilitaires pour manipuler des objets 3d
///@{

/*! \file 
representation d'un objet.

Mesh propose plusieurs manieres de decrire un maillage.
par exemple, un triangle :
\code
Mesh m(GL_TRIANGLES);

// coordonnnees des 3 sommets d'un triangle
Point a= { ... }
Point b= { ... }
Point c= { ... }
m.vertex( a );
m.vertex( b );
m.vertex( c );
\endcode

il est aussi possible de definir d'autres attributs : la couleur du sommet, sa normale, ses coordonnees de texture, ce sont 
les fonctions Mesh::color(), Mesh::normal() et Mesh::texcoord(). on peut decrire ces informations de maniere assez compacte :
\code
Mesh m(GL_TRIANGLES);

m.color(Red()).vertex(a);
m.color(Green()).vertex(b);
m.color(Blue()).vertex(c);
\endcode
et :
\code
Mesh m(GL_TRIANGLES);

m.color(Red()).normal(na).vertex(a);
m.color(Green()).normal(nb).vertex(b);
m.color(Blue()).normal(nc).vertex(c);
\endcode

vertex() doit etre utilise en dernier. 
pour permettre de decrire un triangle de couleur uniforme, par exemple, au lieu d'ecrire :
\code
Mesh m(GL_TRIANGLES);

m.color(Red()).vertex(a);
m.color(Red()).vertex(b);
m.color(Red()).vertex(c);
\endcode
il est plus simple d'ecrire :
\code
Mesh m(GL_TRIANGLES);

m.color(Red());
m.vertex(a);
m.vertex(b);
m.vertex(c);
\endcode
et cette solution permet aussi de decrire plusieurs triangles partageant leurs sommets, par exemple :
\code
Mesh m(GL_TRIANGLES);

// insere 4 sommets dans le maillage et conserve leurs indices
unsigned int a= m.vertex(Point(...));
unsigned int b= m.vertex(Point(...));
unsigned int c= m.vertex(Point(...));
unsigned int d= m.vertex(Point(...));

// decrit 2 triangles avec les indices des 4 sommets
m.triangle(a, b, c);
m.triangle(a, c, d);
\endcode
*/

//! representation d'une matiere.
struct Material
{
    Color diffuse;      //!< couleur diffuse
    Color specular;     //!< couleur du reflet
    Color emission;     //!< pour une source de lumiere
    float ns;           //!< exposant pour les reflets blinn-phong
    
    Material( ) : diffuse(0.8f, 0.8f, 0.8f), specular(Black()), emission(), ns(0) {}
};

//! representation d'un triangle.
struct TriangleData
{
    vec3 a, b, c;       //!< positions       
    vec3 na, nb, nc;    //!< normales
    vec2 ta, tb, tc;    //!< texcoords
};


//! representation d'un objet / maillage.
class Mesh
{
public:
    //! \name construction.
    //@{
    //! constructeur par defaut.
    Mesh( ) : m_positions(), m_texcoords(), m_normals(), m_colors(), m_indices(), m_state_map(), m_state(0),
        m_color(White()), m_primitives(GL_POINTS), m_vao(0), m_buffer(0), m_index_buffer(0), m_program(0), m_update_buffers(false) {}
    
    //! constructeur.
    Mesh( const GLenum primitives ) : m_positions(), m_texcoords(), m_normals(), m_colors(), m_indices(), m_state_map(), m_state(0),
        m_color(White()), m_primitives(primitives), m_vao(0), m_buffer(0), m_index_buffer(0), m_program(0), m_update_buffers(false) {}
    
    //! construit les objets openGL.
    int create( const GLenum primitives );
    //! detruit les objets openGL.
    void release( );
    //@}
    
    //! \name description des attributs des sommets. 
    //@{
    //! definit la couleur du prochain sommet.
    Mesh& color( const vec4& c );
    //! definit la couleur du prochain sommet.
    Mesh& color( const Color& c ) { return color(vec4(c.r, c.g, c.b, c.a)); }
    //! definit la couleur du prochain sommet.
    Mesh& color( const float r, const float g, const float b, const float a= 1) { return color(Color(r, g, b, a)); }
    
    //! definit la normale du prochain sommet.
    Mesh& normal( const vec3& n );
    //! definit la normale du prochain sommet.
    Mesh& normal( const Vector& n ) { return normal(vec3(n)); }
    //! definit la normale du prochain sommet.
    Mesh& normal( const float x, const float y, const float z ) { return normal(vec3(x, y, z)); }
    
    //! definit les coordonnees de texture du prochain sommet.
    Mesh& texcoord( const vec2& uv );
    //! definit les coordonnees de texture du prochain sommet.
    Mesh& texcoord( const float x, const float y ) { return texcoord(vec2(x, y)); }
    
    //! insere un sommet de position p, et ses attributs (s'ils sont definis par color(), texcoord(), normal()), dans l'objet. renvoie l'indice du sommet.
    unsigned int vertex( const vec3& p );
    //! insere un sommet de position p, et ses attributs (s'ils sont definis par color(), texcoord(), normal()), dans l'objet. renvoie l'indice du sommet.
    unsigned int vertex( const Point& p ) { return vertex(vec3(p)); }
    //! insere un sommet de position p, et ses attributs (s'ils sont definis par color(), texcoord(), normal()), dans l'objet. renvoie l'indice du sommet.
    unsigned int vertex( const float x, const float y, const float z ) { return vertex(vec3(x, y, z)); }
    //@}

    //! \name description de triangles indexes.
    //@{
    /*! insere un triangle.  a, b, c sont les indices des sommets deja inseres dans l'objet. ne fonctionne pas avec les strips et les fans.
    \code
    Mesh m(GL_TRIANGLES);
    unsigned int a= m.vertex( Point(ax, ay, az) );
    unsigned int b= m.vertex( Point(bx, by, bz) );
    unsigned int c= m.vertex( Point(cx, cy, cz) );
    m.triangle(a, b, c);
    \endcode
    */
    Mesh& triangle( const unsigned int a, const unsigned int b, const unsigned int c );
    
    /*! insere un triangle, a, b, c sont les indices des sommets deja inseres dans l'objet, en comptant en partant du dernier. ne fonctionne pas avec les strips et les fans.
    \code
    Mesh m(GL_TRIANGLES);
    m.vertex( Point(ax, ay, az) );
    m.vertex( Point(bx, by, bz) );
    m.vertex( Point(cx, cy, cz) );
    m.triangle_last(-3, -2, -1);
    \endcode
    */    
    Mesh& triangle_last( const int a, const int b, const int c );
    
    //! demarre un nouveau strip. a utiliser avec un objet composes de GL_TRIANGLE_STRIP, doit aussi fonctionner avec GL_TRIANGLE_FAN, GL_LINE_STRIP, GL_LINE_LOOP, etc.
    Mesh& restart_strip( );
    //@}
    
    //! \name modification des attributs des sommets.
    //@{
    //! modifie la couleur du sommet d'indice id.
    Mesh& color( const unsigned int id, const vec4& c );
    //! modifie la couleur du sommet d'indice id.
    Mesh& color( const unsigned int id, const Color& c ) { return color(id, vec4(c.r, c.g, c.b, c.a)); }
    //! modifie la couleur du sommet d'indice id.
    Mesh& color( const unsigned int id, const float r, const float g, const float b, const float a= 1) { return color(id, vec4(r, g, b, a)); }

    //! modifie la normale du sommet d'indice id.
    Mesh& normal( const unsigned int id, const vec3& n );
    //! modifie la normale du sommet d'indice id.
    Mesh& normal( const unsigned int id, const Vector& n ) { return normal(id, vec3(n)); }
    //! modifie la normale du sommet d'indice id.
    Mesh& normal( const unsigned int id, const float x, const float y, const float z ) { return normal(id, vec3(x, y, z)); }
    
    //! modifie les coordonnees du sommet d'indice id.
    Mesh& texcoord( const unsigned int id, const vec2& uv );
    //! modifie les coordonnees du sommet d'indice id.
    Mesh& texcoord( const unsigned int id, const float x, const float y ) { return texcoord(id, vec2(x, y)); }
    
    //! modifie la position du sommet d'indice id.
    void vertex( const unsigned int id, const vec3& p );
    //! modifie la position du sommet d'indice id.
    void vertex( const unsigned int id, const Point& p ) { vertex(id, vec3(p)); }
    //! modifie la position du sommet d'indice id.
    void vertex( const unsigned int id, const float x, const float y, const float z ) { vertex(id, vec3(x, y, z)); }
    //@}
    
    //! \name description des matieres.
    //@{
    //! ajoute une description de matiere. et renvoie son indice.
    unsigned int mesh_material( const Material& m );
    //! ajoute un ensemble de description de matieres.
    void mesh_materials( const std::vector<Material>& m );
    
    //! renvoie le nombre de matieres.
    int mesh_material_count( ) const;
    //! renvoie la description de matiere d'indice id.
    const Material& mesh_material( const unsigned int id ) const;
    //! renvoie l'ensemble des matieres.
    const std::vector<Material>& mesh_materials( ) const;
    //! renvoie les indices des matieres des triangles.
    const std::vector<unsigned int>& materials( ) const;
    
    //! definit la matiere du prochain triangle. id est l'indice d'une matiere ajoutee par mesh_material() ou mesh_materials( ). ne fonctionne que pour les primitives GL_TRIANGLES, indexees ou pas.
    Mesh& material( const unsigned int id );
    //@}
    
    //! \name description des triangles d'un maillage.
    //@{
    //! renvoie la matiere d'un triangle.
    const Material &triangle_material( const unsigned int id ) const;
    //! renvoie le nombre de triangles.
    int triangle_count( ) const;
    //! renvoie un triangle.
    TriangleData triangle( const unsigned int id ) const;
    //@}
    
    //! renvoie min et max les coordonnees des extremites des positions des sommets de l'objet (boite englobante alignee sur les axes, aabb).
    void bounds( Point& pmin, Point& pmax );
    
    //! renvoie la couleur par defaut du mesh, utilisee si les sommets n'ont pas de couleur associee.
    Color default_color( ) const { return m_color; }
    //! modifie la couleur par defaut, utilisee si les sommets n'ont pas de couleur associee.
    Mesh& default_color( const Color& color );
    
    //! \name manipulation des buffers d'attributs.
    //@{
    //! renvoie le nombre de sommets.
    int vertex_count( ) const { return (int) m_positions.size(); }
    //! renvoie le nombre d'indices de sommets.
    int index_count( ) const { return (int) m_indices.size(); }
    
    //! renvoie l'adresse de la position du premier sommet. permet de construire les vertex buffers openGL. par convention, la position est un vec3, 3 GL_FLOAT.
    const float *vertex_buffer( ) const { return &m_positions.front().x; }
    //! renvoie la longueur (en octets) du vertex buffer.
    std::size_t vertex_buffer_size( ) const { return m_positions.size() * sizeof(vec3); }

    //! renvoie l'adresse de la normale du premier sommet. par convention, la normale est un vec3, 3 GL_FLOAT.
    const float *normal_buffer( ) const { return &m_normals.front().x; }
    //! renvoie la longueur (en octets) du normal buffer.
    std::size_t normal_buffer_size( ) const { return m_normals.size() * sizeof(vec3); }
    
    //! renvoie l'adresse des coordonnees de textures du premier sommet. par convention, c'est un vec2, 2 GL_FLOAT.
    const float *texcoord_buffer( ) const { return &m_texcoords.front().x; }
    //! renvoie la taille (en octets) du texcoord buffer.
    std::size_t texcoord_buffer_size( ) const { return m_texcoords.size() * sizeof(vec2); }
    
    //! renvoie l'adresse de la couleur du premier sommet. par convention, la couleur est un vec4, 4 GL_FLOAT.
    const float *color_buffer( ) const { return &m_colors.front().x; }
    //! renvoie la taille (en octets) du color buffer.
    std::size_t color_buffer_size( ) const { return m_colors.size() * sizeof(vec4); }
    
    //! renvoie l'adresse du premier indice du premier triangle. par convention c'est un uint, 1, GL_UNSIGNED_INT.
    const void *index_buffer( ) const { return &m_indices.front(); }
    //! renvoie la taille (en octets) de l'index buffer.
    std::size_t index_buffer_size( ) const { return m_indices.size() * sizeof(unsigned int); }
    
    //
    const std::vector<vec3>& positions( ) const { return m_positions; }
    const std::vector<vec2>& texcoords( ) const { return m_texcoords; }
    const std::vector<vec3>& normals( ) const { return m_normals; }
    const std::vector<vec4>& colors( ) const { return m_colors; }
    const std::vector<unsigned int>& indices( ) const { return m_indices; }
    //@}
    
    //! renvoie le type de primitives.
    GLenum primitives( ) const { return m_primitives; }
    
    //! \name gestion d'erreur.
    //@{
    /*! sentinelle pour la gestion d'erreur lors du chargement d'un fichier.
    exemple :
    \code
    Mesh mesh= read_mesh("data/bigguy.obj");
    if(mesh == Mesh::error())
        return "erreur de chargement";
    \endcode
    */
    static Mesh& error( )
    {
        static Mesh mesh;
        return mesh;
    }
    
    bool operator== ( const Mesh& m ) const
    {
        return (this == &m);
    }
    //@}
    
    /*! dessine l'objet avec les transformations model, vue et projection.\n
        eventuellement eclaire l'objet avec une source de lumiere.\n
        eventuellement applique une texture sur la surface de l'objet. 
    
        \param use_light si vrai light et light_color doivent etre definis,
        \param use_texture si vrai texture doit est l'identifiant d'une texture openGL.
    */
    void draw( const Transform& model, const Transform& view, const Transform& projection, 
        const bool use_light, const Point& light, const Color& light_color, 
        const bool use_texture, const GLuint texture,
        const bool use_alpha_test, const float alpha_min );
    
    //! dessine l'objet avec un shader fourni par l'application. les uniforms du shader doivent deja etre configure, cf transformations...
    void draw( const GLuint program, const bool use_position= true, const bool use_texcoord= true, const bool use_normal= true, const bool use_color= true );
    
    //! construit les buffers et le vertex array object necessaires pour dessiner l'objet avec openGL. utilitaire. detruit par release( ).\n
    GLuint create_buffers( const bool use_texcoord= true, const bool use_normal= true, const bool use_color= true );
    
protected:    
    /*! construit un shader program configure.
    \param use_texcoord force l'utilisation des coordonnees de texture
    \param use_normal force l'utilisation des normales
    \param use_color force l'utilisation des couleurs 
    \param use_light force l'utilisation d'un source de lumiere 
    \param use_alpha_test force l'utilisation d'un test de transparence, cf utilisation d'une texture avec un canal alpha
     */
    GLuint create_program( const bool use_texcoord= true, const bool use_normal= true, const bool use_color= true, const bool use_light= false, const bool use_alpha_test= false );
    
    //! modifie les buffers openGL, si necessaire.
    int update_buffers( const bool use_texcoord, const bool use_normal, const bool use_color );
    
    //
    std::vector<vec3> m_positions;
    std::vector<vec2> m_texcoords;
    std::vector<vec3> m_normals;
    std::vector<vec4> m_colors;
    
    std::vector<unsigned int> m_indices;

    std::vector<Material> m_materials;
    std::vector<unsigned int> m_triangle_materials;

    std::unordered_map<unsigned int, GLuint> m_state_map;
    unsigned int m_state;

    Color m_color;
    
    GLenum m_primitives;
    GLuint m_vao;
    GLuint m_buffer;
    GLuint m_index_buffer;
    GLuint m_program;
    
    bool m_update_buffers;
};

///@}
#endif
