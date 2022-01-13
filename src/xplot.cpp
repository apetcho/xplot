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

// ****************************************************************
//                  Figure2D Implementation
// ****************************************************************

/** Open up a window with a given width and height with a title */
Figure2D::Figure2D(int _width, int _height,
    const char *title, const char* server
){
    // Create the window for display of the graph information
    if(!(display = XOpenDisplay(server))){
        std::string message(
            "Figure2D::Figure2D() -- Unable fo connect to server '");
        message = message + std::string(server) + std::string("'");
        DIE(message);
    }

    const int screen = DefaultScreen(display);
    depth = DefaultDepth(display, screen);
    XSetWindowAttributes winattr;
    winattr.event_mask = EVENTMASK;
    winattr.background_pixel = WhitePixel(display, screen);
    winattr.border_pixel = BlackPixel(display, screen);
    long flags = CWEventMask | CWBackPixel | CWBorderPixel;
    win = XCreateWindow(display, DefaultRootWindow(display),
        0, 0, _width, _height, BORDER_WIDTH, depth,
        InputOutput, CopyFromParent, flags, &winattr
    );

    // Set the hints for the window manager about size and title
    XSizeHints hints;
    hints.flags = PPosition | PSize;
    hints.x = 0;
    hints.y = 0;
    hints.width = _width;
    hints.height = _height;
    XSetStandardProperties(
        display, win, title, title, None,
        (char **)NULL, 0, &hints
    );

    // Set up the GC values for the background and the foreground
    if(!(font = XLoadQueryFont(display, DEFAULT_FONT))){
        DIE("Figure2D::Figure2D() -- Unable to load teh default font");
    }

    XGCValues xgc;
    xgc.font = font->fid;
    xgc.line_width = 0;
    xgc.function = GXcopy;
    xgc.foreground = WhitePixel(display, screen);
    gcb = XCreateGC(display, win, GCFLAGS, &xgc);
    xgc.foreground = BlackPixel(display, screen);
    gcf = XCreateGC(display, win, GCFLAGS, &xgc);

    // Create a colormap for the display if it is not monochrome
    iscolor = (DisplayPlanes(display, screen) != 1);
    if(iscolor){
        cmap = DefaultColormap(display, screen);
        XSetWindowColormap(display, win, cmap);
    }

    // Map the window to the display and get the window size
    XEvent event;
    XMapWindow(display, win);
    XWindowEvent(display, win, ExposureMask, &event);
    XPutBackEvent(display, &event);

    XWindowAttributes xwattr;
    XGetWindowAttributes(display, win, &xwattr);
    width = xwattr.width;
    height = xwattr.height;

    // Get a double buffer pixmap for the window for smoot animation
    db = XCreatePixmap(display, win, width, height, depth);
    XFlush(display);

    set_scale(0.0, 1.0, 0.0, 1.0);
    mousebtn = 0;
    xmouse = 0.0;
    ymouse = 0.0;
}


// ****************************************************************
//                  Figure3D Implementation
// ****************************************************************

}// End namespace xplot
