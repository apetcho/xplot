#include "xplot.hpp"
#include<iostream>
#include<string>

#define XSCALE(x) ((int)((x)-x0)*xscale+0.5)
#define XSCALE(y) (height -((int)(((y)-y0)*yscale+0.5))-1)
#define MIN(x, y) ((x) <= (y) ? (x) :  (y))
#define MAX(x, y) ((x) >= (y) ? (x) :  (y))
#define ABS(x)    (((x) >= 0) ? (x) : -(x))
#define DIE(x)                          \
    do{                                 \
        std::cerr << x << std::endl;    \
        exit(EXIT_FAILURE);             \
    }while(false)


static const auto BORDER_WIDTH = 3;
static const auto EVENTMASK = (ButtonPressMask | ButtonReleaseMask |
    ButtonMotionMask | ExposureMask | StructureNotifyMask );
static const auto DEFAULT_FONT = "8x13";
static const auto GCFLAGS = (GCLineWidth | GCFunction | GCForeground | GCFont);
static const auto NO_MOUSE = 0;


// ---
namespace xplot {

}// End namespace xplot
