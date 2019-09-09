
#include <cstdio>
#include <cassert>
#include <string>
#include <algorithm>

#include "vec.h"
#include "mesh.h"

#include "program.h"
#include "uniforms.h"

#include "window.h"


int Mesh::create( const GLenum primitives )
{
    m_primitives= primitives;
    return 0;
}

void Mesh::release( )
{
    glDeleteVertexArrays(1, &m_vao);
    glDeleteBuffers(1, &m_buffer);
    glDeleteBuffers(1, &m_index_buffer);

    // detruit tous les shaders crees...
    for(auto it= m_state_map.begin(); it != m_state_map.end(); ++it)
        if(it->second > 0)
            release_program(it->second);
}

// definit les attributs du prochain sommet
Mesh& Mesh::default_color( const Color& color )
{
    m_color= color;
    return *this;
}

Mesh& Mesh::color( const vec4& color )
{
    m_colors.push_back(color);
    return *this;
}

Mesh& Mesh::normal( const vec3& normal )
{
    m_normals.push_back(normal);
    return *this;
}

Mesh& Mesh::texcoord( const vec2& uv )
{
    m_texcoords.push_back(uv);
    return *this;
}

// insere un nouveau sommet
unsigned int Mesh::vertex( const vec3& position )
{
    m_positions.push_back(position);

    // copie les autres attributs du sommet, uniquement s'ils sont definis
    if(m_texcoords.size() > 0 && m_texcoords.size() != m_positions.size())
        m_texcoords.push_back(m_texcoords.back());
    if(m_normals.size() > 0 && m_normals.size() != m_positions.size())
        m_normals.push_back(m_normals.back());
    if(m_colors.size() > 0 && m_colors.size() != m_positions.size())
        m_colors.push_back(m_colors.back());

    // copie la matiere courante, uniquement si elle est definie
    if(m_triangle_materials.size() > 0 && (size_t) triangle_count() > m_triangle_materials.size())
        m_triangle_materials.push_back(m_triangle_materials.back());
    
    unsigned int index= (unsigned int) m_positions.size() -1;
    // construction de l'index buffer pour les strip
    switch(m_primitives)
    {
        case GL_LINE_STRIP:
        case GL_LINE_LOOP:
        case GL_TRIANGLE_STRIP:
        case GL_TRIANGLE_FAN:
            m_indices.push_back(index);
            break;
        default:
            break;
    }

    // renvoie l'indice du sommet
    return index;
}

// update attributes
Mesh& Mesh::color( const unsigned int id, const vec4& c )
{
    assert(id < m_colors.size());
    m_update_buffers= true;
    m_colors[id]= c;
    return *this;
}

Mesh& Mesh::normal( const unsigned int id, const vec3& n )
{
    assert(id < m_normals.size());
    m_update_buffers= true;
    m_normals[id]= n;
    return *this;
}

Mesh& Mesh::texcoord( const unsigned int id, const vec2& uv )
{
    assert(id < m_texcoords.size());
    m_update_buffers= true;
    m_texcoords[id]= uv;
    return *this;
}

void Mesh::vertex( const unsigned int id, const vec3& p )
{
    assert(id < m_positions.size());
    m_update_buffers= true;
    m_positions[id]= p;
}

//
Mesh& Mesh::triangle( const unsigned int a, const unsigned int b, const unsigned int c )
{
    assert(a < m_positions.size());
    assert(b < m_positions.size());
    assert(c < m_positions.size());
    m_indices.push_back(a);
    m_indices.push_back(b);
    m_indices.push_back(c);
    return *this;
}

Mesh& Mesh::triangle_last( const int a, const int b, const int c )
{
    assert(a < 0);
    assert(b < 0);
    assert(c < 0);
    m_indices.push_back((int) m_positions.size() + a);
    m_indices.push_back((int) m_positions.size() + b);
    m_indices.push_back((int) m_positions.size() + c);
    return *this;
}

Mesh& Mesh::restart_strip( )
{
    m_indices.push_back(~0u);   // ~0u plus grand entier non signe representable
#if 1
    glPrimitiveRestartIndex(~0u);
    glEnable(GL_PRIMITIVE_RESTART);
#else
    glEnable(GL_PRIMITIVE_RESTART_FIXED_INDEX); // n'existe pas sur mac ?!
#endif
    return *this;
}


unsigned int Mesh::mesh_material( const Material& m )
{
    m_materials.push_back(m);
    return (unsigned int) m_materials.size() -1;
}

void Mesh::mesh_materials( const std::vector<Material>& m )
{
    m_materials= m;
}

int Mesh::mesh_material_count( ) const
{
    return (int) m_materials.size();
}

const Material& Mesh::mesh_material( const unsigned int id ) const
{
    assert((size_t) id < m_materials.size());
    return m_materials[id];
}

Mesh& Mesh::material( const unsigned int id )
{
    m_triangle_materials.push_back(id);
    return *this;
}

const Material &Mesh::triangle_material( const unsigned int id ) const
{
    assert((size_t) id < m_triangle_materials.size());
    assert((size_t) m_triangle_materials[id] < m_materials.size());
    return m_materials[m_triangle_materials[id]];
}

const std::vector<Material>& Mesh::mesh_materials( ) const
{
    return m_materials;
}

const std::vector<unsigned int>& Mesh::materials( ) const
{
    return m_triangle_materials;
}

int Mesh::triangle_count( ) const
{
    if(m_primitives != GL_TRIANGLES)
        return 0;
    
    if(m_indices.size() > 0)
        return (int) m_indices.size() / 3;
    else
        return (int) m_positions.size() / 3;
}

TriangleData Mesh::triangle( const unsigned int id ) const
{
    unsigned int a, b, c;
    if(m_indices.size() > 0)
    {
        assert((size_t) id*3+2 < m_indices.size());
        a= m_indices[id*3];
        b= m_indices[id*3 +1];
        c= m_indices[id*3 +2];
    }
    else
    {
        assert((size_t) id*3+2 < m_positions.size());
        a= id*3;
        b= id*3 +1;
        c= id*3 +2;
    }
    
    TriangleData triangle;
    triangle.a= m_positions[a];
    triangle.b= m_positions[b];
    triangle.c= m_positions[c];
    
    if(m_normals.size() == m_positions.size())
    {
        triangle.na= m_normals[a];
        triangle.nb= m_normals[b];
        triangle.nc= m_normals[c];
    }
    else
    {
        // calculer la normale geometrique
        Vector ab= Point(m_positions[b]) - Point(m_positions[a]);
        Vector ac= Point(m_positions[c]) - Point(m_positions[a]);
        Vector n= normalize(cross(ab, ac));
        triangle.na= vec3(n);
        triangle.nb= vec3(n);
        triangle.nc= vec3(n);
    }
    
    if(m_texcoords.size() == m_positions.size())
    {
        triangle.ta= m_texcoords[a];
        triangle.tb= m_texcoords[b];
        triangle.tc= m_texcoords[c];
    }
    else
    {
        // coordonnees barycentriques des sommets, convention p(u, v)= w*a + u*b + v*c, avec w= 1 - u -v
        triangle.ta= vec2(0, 0);        // w= 1
        triangle.tb= vec2(1, 0);        // w= 0
        triangle.tc= vec2(0, 1);        // w= 0
    }
    
    return triangle;
}

void Mesh::bounds( Point& pmin, Point& pmax )
{
    if(m_positions.size() < 1)
        return;

    pmin= Point(m_positions[0]);
    pmax= pmin;

    for(unsigned int i= 1; i < (unsigned int) m_positions.size(); i++)
    {
        vec3 p= m_positions[i];
        pmin= Point( std::min(pmin.x, p.x), std::min(pmin.y, p.y), std::min(pmin.z, p.z) );
        pmax= Point( std::max(pmax.x, p.x), std::max(pmax.y, p.y), std::max(pmax.z, p.z) );
    }
}

GLuint Mesh::create_buffers( const bool use_texcoord, const bool use_normal, const bool use_color )
{
    if(m_positions.size() == 0)
        return 0;

#if 1
    if(m_texcoords.size() > 0 && m_texcoords.size() < m_positions.size() && use_texcoord)
        printf("[oops] mesh: invalid texcoords array...\n");
    if(m_normals.size() > 0 && m_normals.size() < m_positions.size() && use_normal)
        printf("[oops] mesh: invalid normals array...\n");
    if(m_colors.size() > 0 && m_colors.size() < m_positions.size() && use_color)
        printf("[oops] mesh: invalid colors array...\n");
#endif
    
    glGenVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);
    
    // determine la taille du buffer pour stocker tous les attributs et les indices
    size_t size= vertex_buffer_size() + texcoord_buffer_size() + normal_buffer_size() + color_buffer_size();
    // allouer le buffer
    glGenBuffers(1, &m_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, m_buffer);
    glBufferData(GL_ARRAY_BUFFER, size, nullptr, GL_STATIC_DRAW);
    
    // transferer les attributs et configurer le format de sommet (vao)
    size_t offset= 0;
    size= vertex_buffer_size();
    glBufferSubData(GL_ARRAY_BUFFER, offset, size, vertex_buffer());
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const void *) offset);
    glEnableVertexAttribArray(0);
    
    if(m_texcoords.size() == m_positions.size() && use_texcoord)
    {
        offset= offset + size;
        size= texcoord_buffer_size();
        glBufferSubData(GL_ARRAY_BUFFER, offset, size, texcoord_buffer());
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (const void *) offset);
        glEnableVertexAttribArray(1);
    }
    
    if(m_normals.size() == m_positions.size() && use_normal)
    {
        offset= offset + size;
        size= normal_buffer_size();
        glBufferSubData(GL_ARRAY_BUFFER, offset, size, normal_buffer());
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, (const void *) offset);
        glEnableVertexAttribArray(2);
    }
    
    if(m_colors.size() == m_positions.size() && use_color)
    {
        offset= offset + size;
        size= color_buffer_size();
        glBufferSubData(GL_ARRAY_BUFFER, offset, size, color_buffer());
        glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 0, (const void *) offset);
        glEnableVertexAttribArray(3);
    }
    
    // allouer l'index buffer
    if(index_buffer_size())
    {
        glGenBuffers(1, &m_index_buffer);    
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_index_buffer);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, index_buffer_size(), index_buffer(), GL_STATIC_DRAW);
    }

    m_update_buffers= false;
    
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    return m_vao;
}

int Mesh::update_buffers( const bool use_texcoord, const bool use_normal, const bool use_color )
{
    assert(m_vao > 0);
    assert(m_buffer > 0);
    if(!m_update_buffers)
        return 0;
    
    glBindBuffer(GL_ARRAY_BUFFER, m_buffer);
    
    size_t offset= 0;
    size_t size= vertex_buffer_size();
    glBufferSubData(GL_ARRAY_BUFFER, offset, size, vertex_buffer());
    
    if(m_texcoords.size() == m_positions.size() && use_texcoord)
    {
        offset= offset + size;
        size= texcoord_buffer_size();
        glBufferSubData(GL_ARRAY_BUFFER, offset, size, texcoord_buffer());
    }
    
    if(m_normals.size() == m_positions.size() && use_normal)
    {
        offset= offset + size;
        size= normal_buffer_size();
        glBufferSubData(GL_ARRAY_BUFFER, offset, size, normal_buffer());
    }
    
    if(m_colors.size() == m_positions.size() && use_color)
    {
        offset= offset + size;
        size= color_buffer_size();
        glBufferSubData(GL_ARRAY_BUFFER, offset, size, color_buffer());
    }
    
    m_update_buffers= false;
    return 1;
}


GLuint Mesh::create_program( const bool use_texcoord, const bool use_normal, const bool use_color, const bool use_light, const bool use_alpha_test )
{
    std::string definitions;

    if(m_texcoords.size() == m_positions.size() && use_texcoord)
        definitions.append("#define USE_TEXCOORD\n");
    if(m_normals.size() == m_positions.size() && use_normal)
        definitions.append("#define USE_NORMAL\n");
    if(m_colors.size() == m_positions.size() && use_color)
        definitions.append("#define USE_COLOR\n");
    if(use_light)
        definitions.append("#define USE_LIGHT\n");
    if(use_texcoord && use_alpha_test)
        definitions.append("#define USE_ALPHATEST\n");

    //~ printf("--\n%s", definitions.c_str());
    bool use_mesh_color= (m_primitives == GL_POINTS || m_primitives == GL_LINES || m_primitives == GL_LINE_STRIP || m_primitives == GL_LINE_LOOP);
    if(!use_mesh_color)
        m_program= read_program( smart_path("data/shaders/mesh.glsl"), definitions.c_str());
    else
        m_program= read_program( smart_path("data/shaders/mesh_color.glsl"), definitions.c_str());
    return m_program;
}


void Mesh::draw( const Transform& model, const Transform& view, const Transform& projection,
    const bool use_light, const Point& light, const Color& light_color,
    const bool use_texture, const GLuint texture,
    const bool use_alpha_test, const float alpha_min )
{
    bool use_texcoord= (m_texcoords.size() == m_positions.size() && texture > 0);
    bool use_normal= (m_normals.size() == m_positions.size());
    bool use_color= (m_colors.size() == m_positions.size());
    
    // etape 1 : construit le program en fonction des attributs du mesh et des options choisies
    unsigned int key= 0;
    if(use_texcoord) key= key | 1;
    if(use_normal) key= key | 2;
    if(use_color) key= key | 4;
    if(use_texture) key= key | 8;
    if(use_light) key= key | 16;
    if(use_alpha_test) key= key | 32;

    if(m_state != key)
        // recherche un shader deja compile pour ce type de draw
        m_program= m_state_map[key];

    if(m_program == 0)
    {
        // pas de shader pour ce type de draw
        create_program(use_texcoord, use_normal, use_color, use_light, use_alpha_test);
        program_print_errors(m_program);

        // conserver le shader
        m_state_map[key]= m_program;
    }

    // conserve la config du shader selectionne.
    assert(m_program != 0);
    m_state= key;

    glUseProgram(m_program);
    program_uniform(m_program, "mesh_color", default_color());

    Transform mv= view * model;
    Transform mvp= projection * mv;

    program_uniform(m_program, "mvpMatrix", mvp);
    program_uniform(m_program, "mvMatrix", mv);
    if(use_normal)
        program_uniform(m_program, "normalMatrix", mv.normal()); // transforme les normales dans le repere camera.

    // utiliser une texture, elle ne sera visible que si le mesh a des texcoords...
    if(texture && use_texcoord && use_texture)
        program_use_texture(m_program, "diffuse_color", 0, texture);

    if(use_light)
    {
        program_uniform(m_program, "light", view(light));       // transforme la position de la source dans le repere camera, comme les normales
        program_uniform(m_program, "light_color", light_color);
    }
    
    if(use_alpha_test)
        program_uniform(m_program, "alpha_min", alpha_min);
    
    // etape  2 : cree les buffers et le vao
    if(m_vao == 0)
        create_buffers(true, true, true);
    
    assert(m_vao != 0);
    if(m_update_buffers)
        update_buffers(true, true, true);
    
    glBindVertexArray(m_vao);
    
    // etape 3 : dessiner
    if(m_indices.size() > 0)
        glDrawElements(m_primitives, (GLsizei) m_indices.size(), GL_UNSIGNED_INT, 0);
    else
        glDrawArrays(m_primitives, 0, (GLsizei) m_positions.size());
}

void Mesh::draw( const GLuint program, const bool use_position, const bool use_texcoord, const bool use_normal, const bool use_color )
{
    if(program == 0)
    {
        printf("[oops]  no program... can't draw !!");
        return;
    }
    
    if(m_vao == 0)
        // cree les buffers demandes, inclus toujours position
        create_buffers(use_texcoord, use_normal, use_color);
    assert(m_vao != 0);
    
    if(m_update_buffers)
        update_buffers(use_texcoord, use_normal, use_color);
    
    glBindVertexArray(m_vao);
    
#ifndef GK_RELEASE
    // verifie que le program est selectionne
    GLuint current;
    glGetIntegerv(GL_CURRENT_PROGRAM, (GLint *) &current);
    if(current != program)
    {
    #ifdef GL_VERSION_4_3
        {
            char label[1024];
            glGetObjectLabel(GL_PROGRAM, program, sizeof(label), nullptr, label);
            
            printf("[oops] program( %u '%s' ): not active... can't draw !!\n", program, label); 
        }
    #else
        printf("[oops] program( %u ): not active... can't draw !!", program); 
    #endif
    }
    
    // verifie que les attributs necessaires a l'execution du shader sont presents dans le mesh...

    // etape 1 : recuperer le nombre d'attributs
    GLint n= 0;
    glGetProgramiv(program, GL_ACTIVE_ATTRIBUTES, &n);
    
    // etape 2 : recuperer les infos de chaque attribut
    char name[1024];
    for(int index= 0; index < n; index++)
    {
        GLint glsl_size;
        GLenum glsl_type;
        glGetActiveAttrib(program, index, sizeof(name), nullptr, &glsl_size, &glsl_type, name);

        GLint location= glGetAttribLocation(program, name);
        if(location == 0)       // attribut position necessaire a l'execution du shader
        {
            if(!use_position || !vertex_buffer_size())
                printf("[oops]  no position '%s' attribute in mesh... can't draw !!\n", name);
            if(glsl_size != 1 || glsl_type != GL_FLOAT_VEC3)
                printf("[oops]  attribute '%s' is not declared as a vec3... can't draw !!\n", name);
        }
        else if(location == 1)  // attribut texcoord necessaire 
        {
            if(!use_texcoord || !texcoord_buffer_size())
                printf("[oops]  no texcoord '%s' attribute in mesh... can't draw !!\n", name);
            if(glsl_size != 1 || glsl_type != GL_FLOAT_VEC2)
                printf("[oops]  attribute '%s' is not declared as a vec2... can't draw !!\n", name);
        }
        else if(location == 2)  // attribut normal necessaire
        {
            if(!use_normal || !normal_buffer_size())
                printf("[oops]  no normal '%s' attribute in mesh... can't draw !!\n", name);
            if(glsl_size != 1 || glsl_type != GL_FLOAT_VEC3)
                printf("[oops]  attribute '%s' is not declared as a vec3... can't draw !!\n", name);
        }
        else if(location == 3)  // attribut color necessaire
        {
            if(!use_color || !color_buffer_size())
                printf("[oops]  no color '%s' attribute in mesh... can't draw !!\n", name);
            if(glsl_size != 1 || glsl_type != GL_FLOAT_VEC4)
                printf("[oops]  attribute '%s' is not declared as a vec4... can't draw !!\n", name);
        }
    }
#endif
    
    if(m_indices.size() > 0)
        glDrawElements(m_primitives, (GLsizei) m_indices.size(), GL_UNSIGNED_INT, 0);
    else
        glDrawArrays(m_primitives, 0, (GLsizei) m_positions.size());
}
