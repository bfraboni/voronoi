
#include "app_camera.h"

AppCamera::AppCamera( const int width, const int height, const int major, const int minor ) : App(width, height, major, minor), m_camera() {}

AppCamera::~AppCamera( ) {}

int AppCamera::update( const float time, const float delta )
{
    // recupere les mouvements de la souris
    int mx, my;
    unsigned int mb= SDL_GetRelativeMouseState(&mx, &my);
    int mousex, mousey;
    SDL_GetMouseState(&mousex, &mousey);
    
    // deplace la camera
    if(mb & SDL_BUTTON(1))
        m_camera.rotation(mx, my);      // tourne autour de l'objet
    else if(mb & SDL_BUTTON(3))
        m_camera.translation((float) mx / (float) window_width(), (float) my / (float) window_height()); // deplace le point de rotation
    else if(mb & SDL_BUTTON(2))
        m_camera.move(mx);           // approche / eloigne l'objet
    
    SDL_MouseWheelEvent wheel= wheel_event();
    if(wheel.y != 0)
    {
        clear_wheel_event();
        m_camera.move(8.f * wheel.y);  // approche / eloigne l'objet
    }
    
    return 0;
}
