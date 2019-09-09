
#ifndef _COLOR_H
#define _COLOR_H


//! \addtogroup image
///@{

//! \file
//! manipulation de couleurs

//! representation d'une couleur (rgba) transparente ou opaque.
struct Color
{
    //! constructeur par defaut.
    Color( ) : r(0.f), g(0.f), b(0.f), a(1.f) {}
    explicit Color( const float _r, const float _g, const float _b, const float _a= 1.f ) : r(_r), g(_g), b(_b), a(_a) {}
    explicit Color( const float _value ) : r(_value), g(_value), b(_value), a(1.f) {}
    
    //! cree une couleur avec les memes composantes que color, mais remplace sa composante alpha (color.r, color.g, color.b, alpha).
    Color( const Color& color, const float alpha ) : r(color.r), g(color.g), b(color.b), a(alpha) {}  // remplace alpha.
    
    float power( ) const;
    
    float r, g, b, a;
};

//! utilitaire. renvoie une couleur noire.
Color Black( );
//! utilitaire. renvoie une couleur blanche.
Color White( );
//! utilitaire. renvoie une couleur rouge.
Color Red( );
//! utilitaire. renvoie une couleur verte.
Color Green( );
//! utilitaire. renvoie une couleur bleue.
Color Blue( );
//! utilitaire. renvoie une couleur jaune.
Color Yellow( );

Color operator+ ( const Color& a, const Color& b );
Color operator- ( const Color& a, const Color& b );
Color operator- ( const Color& c );
Color operator* ( const Color& a, const Color& b );
Color operator* ( const Color& c, const float k );
Color operator* ( const float k, const Color& c );
Color operator/ ( const Color& a, const Color& b );
Color operator/ ( const float k, const Color& c );
Color operator/ ( const Color& c, const float k );

///@}
#endif
