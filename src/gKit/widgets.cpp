
#include <cstdio>
#include <cassert>
#include <cstring>
#include <vector>
#include <algorithm>

#include "window.h"
#include "widgets.h"


Widgets create_widgets( )
{
    Widgets w;
    w.console= create_text();
    w.px= 0; w.py= 0;
    w.focus= 0; w.fx= 0; w.fy= 0;
    w.mb= 0; w.mx= 0; w.my= 0;
    w.wx= 0; w.wy= 0;
    return w;
}

void release_widgets( Widgets& w )
{
    release_text(w.console);
}


void begin( Widgets& w )
{
    clear(w.console);
    w.px= 0;
    w.py= 0;

    SDL_MouseButtonEvent mouse= button_event();
    w.mb= 0;
    w.mx= mouse.x / 8;
    w.my= mouse.y / 16;
    if(mouse.state == SDL_PRESSED)
    {
        clear_button_event();
        w.mb= 1;
    }
    
    SDL_MouseWheelEvent wheel= wheel_event();
    w.wx= 0;
    w.wy= 0;
    if(wheel.x != 0 || wheel.y != 0)
    {
        clear_wheel_event();
        w.wx= wheel.x;
        w.wy= wheel.y;
    }
    
    SDL_KeyboardEvent key= key_event( );
    w.key= 0;
    w.mod= 0;
    if(key.type == SDL_KEYDOWN)
    {
        clear_key_event();
        w.key= key.keysym.sym;
        w.mod= key.keysym.mod;
        
        // filtre les touches speciales
        switch(w.key)
        {
            case SDLK_SPACE:
            case SDLK_BACKSPACE:
            case SDLK_DELETE:
            case SDLK_UP:
            case SDLK_DOWN:
            case SDLK_PAGEUP:
            case SDLK_PAGEDOWN:
            case SDLK_LEFT:
            case SDLK_RIGHT:
            case SDLK_RETURN:
                break;
            default:
                w.key= 0;
        }
    }
    
    SDL_TextInputEvent input= text_event();
    if(input.text[0] != 0)
    {
        w.key= input.text[0];
        clear_text_event();
    }
}


struct Rect
{
    int x, y;
    int w, h;
};

static
bool overlap( const Rect r, const int x, const int y )
{
    return (x >= r.x && x < r.x + r.w && y >= r.y && y < r.y + r.h);
}

static
Rect place( Widgets& w, const int width, const int height= 1 )
{
    Rect r;
    r.x= w.px;
    r.y= w.py;
    r.w= width;
    r.h= height;
    
    if(height == 1)
        // place le prochain widget a droite 
        w.px= r.x + r.w +2; // +2 marge
    else
        // place le prochain widget sur la ligne suivante
        w.py= w.py + height;

    return r;
}

static
Rect place( Widgets& w, const char *text )
{
    return place(w, (int) strlen(text));
}

void begin_line( Widgets& w )
{
    // place le prochain widget sur une nouvelle ligne
    w.px= 0;
    w.py= w.py +1;
}

void end_line( Widgets& widgets )
{
    return;
}


void label( Widgets& w, const char *format, ... )
{
    char tmp[4096] = { 0 };
    
    va_list args;
    va_start(args, format);
    vsnprintf(tmp, sizeof(tmp), format, args);
    va_end(args);
    
    Rect r= place(w, tmp);
    print(w.console, r.x, r.y, tmp);
}

bool button( Widgets& w, const char *text, int& status )
{
    Rect r= place(w, (int) strlen(text) +2);
    
    bool change= false;
    if(w.mb > 0 && overlap(r, w.mx, w.my))
    {
        change= true;
        status= (status + 1) % 2;
    }

    char tmp[128];
    sprintf(tmp, "%c %s", (status != 0) ? 22 : 20, text);  // strlen(text) + 2

    print(w.console, r.x, r.y, tmp);
    return change;
}

bool select( Widgets& w, const char *text, const int option, int& status )
{
    Rect r= place(w, (int) strlen(text) +2);
    
    bool change= false;
    if(w.mb > 0 && overlap(r, w.mx, w.my))
    {
        change= true;
        status= option;
    }

    char tmp[128];
    sprintf(tmp, "%c %s", (status == option) ? 4 : 3, text);  // strlen(text) + 2

    print(w.console, r.x, r.y, tmp);
    return change;
}

bool value( Widgets& w, const char *label, int& value, const int value_min, const int value_max, const int value_step )
{
    char tmp[128];
    sprintf(tmp, "%s: %4d", label, value);
    
    Rect r= place(w, (int) strlen(tmp));

    bool change= false;
    if(w.mb > 0 && overlap(r, w.mx, w.my))
        change= true;
    if(overlap(r, w.mx, w.my))
    {
        if(w.wy != 0)
            value= value + w.wy * value_step;
        if(w.key == SDLK_UP)
            value= value + value_step;
        if(w.key == SDLK_PAGEUP)
            value= value + value_step * 10;
        if(w.key == SDLK_DOWN)
            value= value - value_step;
        if(w.key == SDLK_PAGEDOWN)
            value= value - value_step * 10;
        
        if(value < value_min)
            value= value_min;
        if(value > value_max)
            value= value_max;
        
        sprintf(tmp, "%s: ", label);
        int l= (int) strlen(tmp);
        print(w.console, r.x, r.y, tmp);
        
        sprintf(tmp, "%4d", value);
        print_background(w.console, r.x + l, r.y, tmp);
    }
    else
        print(w.console, r.x, r.y, tmp);
    
    return change;
}

bool value( Widgets& w, const char *label, float& value, const float value_min, const float value_max, const float value_step )
{
    char tmp[128];
    sprintf(tmp, "%s: %.6f", label, value);
    
    Rect r= place(w, (int) strlen(tmp));

    bool change= false;
    if(w.mb > 0 && overlap(r, w.mx, w.my))
        change= true;
    if(overlap(r, w.mx, w.my))
    {
        if(w.wy != 0)
            value= value + w.wy * value_step;
        if(w.key == SDLK_UP)
            value= value + value_step;
        if(w.key == SDLK_PAGEUP)
            value= value + value_step * 10;
        if(w.key == SDLK_DOWN)
            value= value - value_step;
        if(w.key == SDLK_PAGEDOWN)
            value= value - value_step * 10;
        
        if(value < value_min)
            value= value_min;
        if(value > value_max)
            value= value_max;
        
        sprintf(tmp, "%s: ", label);
        int l= (int) strlen(tmp);
        print(w.console, r.x, r.y, tmp);
        
        sprintf(tmp, "%.6f", value);
        print_background(w.console, r.x + l, r.y, tmp);
    }
    else
        print(w.console, r.x, r.y, tmp);
    
    return change;
}

void text_area( Widgets& w, const int height, const char *text, int& begin_line )
{
    Rect r= place(w, 128, height);
    
    if(overlap(r, w.mx, w.my))
    {
        if(w.wy != 0)
            begin_line= begin_line - w.wy;
        if(w.key == SDLK_PAGEUP)
            begin_line= begin_line - height;
        if(w.key == SDLK_PAGEDOWN)
            begin_line= begin_line + height;
        if(w.key == SDLK_SPACE)
            begin_line= begin_line + height / 2;
        if(w.key == SDLK_UP)
            begin_line= begin_line - 1;
        if(w.key == SDLK_DOWN)
            begin_line= begin_line + 1;
    }
    
    // compter les lignes
    int n= 1;
    for(int i= 0; text[i] != 0; i++)
        if(text[i] == '\n')
            n++;
    // ligne de depart, pout que tout le texte reste affiche
    if(begin_line + height > n)
        begin_line= n - height;
    if(begin_line < 1)
        begin_line= 1;
    // retrouve le debut de la ligne
    int line= 1;
    int offset= 0;
    for(int i= 0; text[i] != 0; i++)
    {
        if(text[i] == '\n')
        {
            line++;
            if(line == begin_line)
                offset= i;
        }
    }
    // affiche le texte
    print(w.console, r.x, r.y, text + offset);
}

bool edit( Widgets& w, int text_size, char *text )
{
    assert(text_size > 1);
    int size= std::min((int) strlen(text), text_size -2);
    Rect r= place(w, text_size -1);
    
    // focus
    bool change= false;
    if(w.mb > 0)
    {
        if(overlap(r, w.mx, w.my))
        {
            w.focus= 1;
            w.fx= w.mx;
            w.fy= w.my;
            
            if(w.fx >= r.x + size)
                w.fx= r.x + size;
        }
        
        else
            change= (w.focus > 0);      // click en dehors de la zone editable
    }
    
    // edition
    bool focus= overlap(r, w.fx, w.fy);
    if(focus && w.key > 0)
    {
        int c= w.fx - r.x;
        assert(c < text_size -1);
        
        if(w.key == SDLK_BACKSPACE)
        {
            w.fx--;     // curseur a gauche
            for(int i= c -1; i >= 0 && i+1 < text_size; i++) text[i]= text[i+1];
        }
        else if(w.key == SDLK_DELETE)
        {
            // curseur ne bouge pas
            for(int i= c; i+1 < text_size; i++) text[i]= text[i+1];
        }
        else if(w.key == SDLK_LEFT)
        {
            w.fx--;     // curseur a gauche
        }
        else if(w.key == SDLK_RIGHT)
        {
            w.fx++;     // curseur a droite
        }
        else if(w.key == SDLK_RETURN)
        {
            w.focus= 0;
            change= true;
        }
        else
        {
            w.fx++;     // curseur a droite
            for(int i= text_size -1; i > c; i--) text[i]= text[i -1];
            text[c]= w.key;

            if(size < text_size - 2)
                size++;
            text[size+1]= 0;
        }
        
        // verifier que le curseur reste dans la zone editable
        if(w.fx < r.x)
            w.fx= r.x;
        if(w.fx >= r.x + size)
            w.fx= r.x + size;
    }
    
    int i= 0;
    char tmp[128];
    for(; text[i] != 0; i++) tmp[i]= text[i];   // copie les caracteres
    for(; i < text_size; i++) tmp[i]= ' ';      // complete avec des espaces
    tmp[text_size -1]= 0;                       // termine avec un 0

    print_background(w.console, r.x, r.y, tmp);
    if(focus)
        print_background(w.console, w.fx, w.fy, text[w.fx - r.x], 1);
    
    return change;
}

void end( Widgets& w )
{
    return;
}

void default_color( Widgets& w, const Color& color )
{
    default_color(w.console, color);
}

void draw( Widgets& w, const int width, const int height )
{
    draw(w.console, width, height);
}


