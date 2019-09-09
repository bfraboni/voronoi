
#include <cassert>
#include <vector>
#include <stdio.h>

#include "gamepads.h"

Gamepads::Gamepads( ) {}

bool Gamepads::create( )
{
    // mapping cf https://github.com/gabomdq/SDL_GameControllerDB
    SDL_GameControllerAddMappingsFromFile("gamecontrollerdb.txt");
    
    int n= SDL_NumJoysticks();
    for(int i= 0; i < n; i++)
    {
        if(SDL_IsGameController(i))
        {
            SDL_GameController *pad= SDL_GameControllerOpen(i);
            if(pad == nullptr)
            {
                printf("[error] opening pad: %s\n", SDL_GetError());
                continue;
            }
            
            m_pads.push_back(Gamepad(pad));
        }
    }
    
    printf("found %d pads...\n", (int) m_pads.size());
    return (m_pads.size() > 0);
}

Gamepads::~Gamepads( )
{
    release();
}

void Gamepads::release( )
{
    for(size_t i= 0; i < m_pads.size(); i++)
        if(m_pads[i].m_pad != NULL)
            SDL_GameControllerClose(m_pads[i].m_pad);
    
    m_pads.clear();
}

void Gamepads::update( )
{
    for(size_t p= 0; p < m_pads.size(); p++)
    {
        if(SDL_GameControllerGetAttached(m_pads[p].m_pad) == SDL_FALSE)
        {
            // le pad est debranche...         
            for(int i= 0; i < SDL_CONTROLLER_BUTTON_MAX; i++)
                m_pads[p].m_buttons[i]= 0;
            for(int i= 0; i < SDL_CONTROLLER_AXIS_MAX; i++)
                m_pads[p].m_axis[i]= 0;
            
            SDL_GameControllerClose(m_pads[p].m_pad);
            m_pads[p].m_pad= NULL;
        }
        
        SDL_GameController *pad= m_pads[p].m_pad;
        if(pad == NULL) 
            continue;
        
        // boutons
        for(int i= 0; i < SDL_CONTROLLER_BUTTON_MAX; i++)
            if(SDL_GameControllerGetButton(pad, (SDL_GameControllerButton) i) == SDL_PRESSED)
                m_pads[p].m_buttons[i]= 1;
            else
                m_pads[p].m_buttons[i]= 0;
        
        // axes
        for(int i= 0; i < SDL_CONTROLLER_AXIS_MAX; i++)
        {
            int value= SDL_GameControllerGetAxis(pad, (SDL_GameControllerAxis) i);
            if(value > -6000 && value < 6000)
                // dead zone...
                value= 0;
            
            m_pads[p].m_axis[i]= (float) value / 32768.f;
        }
    }
}


int Gamepads::pads( )
{
    return (int) m_pads.size();
}

Gamepad& Gamepads::pad( const unsigned int index )
{
    assert(index < m_pads.size());
    return m_pads[index];
}

int Gamepads::button( const unsigned int index, const SDL_GameControllerButton b )
{
    return pad(index).button(b);
}

void Gamepads::clear_button( const unsigned int index, const SDL_GameControllerButton b )
{
    return pad(index).clear_button(b);
}

float Gamepads::axis( const unsigned int index, const SDL_GameControllerAxis a )
{
    return pad(index).axis(a);
}

void Gamepads::clear_axis( const unsigned int index, const SDL_GameControllerAxis a )
{
    return pad(index).clear_axis(a);
}

bool Gamepad::connected( )
{
    return (m_pad != NULL);
}

int Gamepad::button( const SDL_GameControllerButton b )
{
    assert(b < SDL_CONTROLLER_BUTTON_MAX);
    return m_buttons[b];
}

void Gamepad::clear_button( const SDL_GameControllerButton b )
{
    assert(b < SDL_CONTROLLER_BUTTON_MAX);
    m_buttons[b]= 0;
}

float Gamepad::axis( const SDL_GameControllerAxis a )
{
    assert(a < SDL_CONTROLLER_AXIS_MAX);
    return m_axis[a];
}

void Gamepad::clear_axis( const SDL_GameControllerAxis a )
{
    assert(a < SDL_CONTROLLER_AXIS_MAX);
    m_axis[a]= 0;
}

