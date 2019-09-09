
#ifndef _PROGRAM_H
#define _PROGRAM_H

#include <string>

#include "glcore.h"


//! \addtogroup openGL utilitaires openGL
///@{

//! \file
//! shader program openGL

//! cree un shader program. a detruire avec release_program( ).\n
//! charge un seul fichier, les shaders sont separes par \#ifdef VERTEX_SHADER / \#endif et \#ifdef FRAGMENT_SHADER / \#endif.\n
//! renvoie l'identifiant openGL du program et le program est selectionne (cf glUseProgram( )).
//! \param filename nom du fichier source.
//! \param definitions chaine de caracteres pouvant comporter plusieurs lignes "#define what value\n".
GLuint read_program( const char *filename, const char *definitions= "" );

//! detruit les shaders et le program.
int release_program( const GLuint program );

//! recharge les sources et recompile un shader program.
//! \param program shader program a modifier
//! \param filename nom du fichier source a charger
//! \param definitions cf read_program
int reload_program( const GLuint program, const char *filename, const char *definitions= "" );

//! renvoie les erreurs de compilation.
int program_format_errors( const GLuint program, std::string& errors );

//! affiche les erreurs de compilation.
int program_print_errors( const GLuint program );

///@}
#endif
