
#ifndef _MAT_H
#define _MAT_H

#include "vec.h"


//! \addtogroup math manipulations de points, vecteur, matrices, transformations
///@{

//! \file
//! transformation de points et vecteurs

//! conversion en radians.
float radians( const float deg );
//! conversion en degres.
float degrees( const float rad );

//! representation d'une transformation, une matrice 4x4, organisee par ligne / row major.
struct Transform
{
    //! constructeur.
    Transform (
        const float t00=1.f, const float t01=0.f, const float t02=0.f, const float t03=0.f,
        const float t10=0.f, const float t11=1.f, const float t12=0.f, const float t13=0.f,
        const float t20=0.f, const float t21=0.f, const float t22=1.f, const float t23=0.f,
        const float t30=0.f, const float t31=0.f, const float t32=0.f, const float t33=1.f );
    
	//! constructeur a partir de 4 Vector colonnes, met (0, 0, 0, 1) dans la 4e ligne 
	Transform(const Vector& x, const Vector& y, const Vector& z, const Vector& w);

	//! renvoie le Vector colonne c de la matrice
	Vector operator[](int c) const;

    //! renvoie le point transforme.
    Point operator() ( const Point& p ) const;
    //! renvoie le vecteur transforme.
    Vector operator() ( const Vector& v ) const;
    //! renvoie le point/vecteur homogene transforme.
    vec4 operator() ( const vec4& v ) const;
    
    //! renvoie la composition de la transformation this et b, t = this * b. permet de transformer un point sans "ambiguite" Point q= a(b(c(p)));
    Transform operator() ( const Transform& b ) const;
    
    //! renvoie la transposee de la matrice.
    Transform transpose( ) const;
    //! renvoie l'inverse de la matrice.
    Transform inverse( ) const;
    //! renvoie la transformation a appliquer aux normales d'un objet transforme par la matrice m.
    Transform normal( ) const;  
    
    //! renvoie l'adresse de la premiere valeur de la matrice.
    const float *buffer( ) const { return &m[0][0]; }
    
    float m[4][4];
};

//! construit la transformation identite.
Transform Identity( );

//! renvoie la transposee de la matrice.
Transform Transpose( const Transform& m );
//! renvoie l'inverse de la matrice.
Transform Inverse( const Transform& m );
//! renvoie la transformation a appliquer aux normales d'un objet transforme par la matrice m.
Transform Normal( const Transform& m );

//! renvoie la matrice representant une mise a l'echelle / etirement.
Transform Scale( const float x, const float y, const float z );

//! renvoie la matrice representant une translation par un vecteur.
Transform Translation( const Vector& v );
//! renvoie la matrice representant une translation par un vecteur x y z.
Transform Translation( const float x, const float y, const float z );

//! renvoie la matrice representation une rotation de angle degree autour de l'axe X.
Transform RotationX( const float angle );
//! renvoie la matrice representation une rotation de a degree autour de l'axe Y.
Transform RotationY( const float angle );
//! renvoie la matrice representation une rotation de angle degree autour de l'axe Z.
Transform RotationZ( const float angle );
//! renvoie la matrice representation une rotation de angle degree autour de l'axe axis.
Transform Rotation( const Vector& axis, const float angle );

//! renvoie la matrice representant une transformation viewport.
Transform Viewport( const float width, const float height );
//! renvoie la matrice representant une transformation projection perspective.
Transform Perspective( const float fov, const float aspect, const float znear, const float zfar );
//! renvoie la matrice representant le placement et l'orientation d'une camera pour observer le point to.
Transform Lookat( const Point& from, const Point& to, const Vector& up );

//! renvoie la composition des transformations a et b, t= a * b.
Transform compose_transform( const Transform& a, const Transform& b );
//! renvoie la composition des transformations a et b, t = a * b.
Transform operator* ( const Transform& a, const Transform& b );

#include <iostream>

inline std::ostream& operator<<(std::ostream& o, const Transform& t)
{
	o << t.m[0][0] << " " << t.m[0][1] << " " << t.m[0][2] << " " << t.m[0][3] << " " << std::endl;
	o << t.m[1][0] << " " << t.m[1][1] << " " << t.m[1][2] << " " << t.m[1][3] << " " << std::endl;
	o << t.m[2][0] << " " << t.m[2][1] << " " << t.m[2][2] << " " << t.m[2][3] << " " << std::endl;
	o << t.m[3][0] << " " << t.m[3][1] << " " << t.m[3][2] << " " << t.m[3][3] << " " << std::endl;
	return o;
}


///@}
#endif
