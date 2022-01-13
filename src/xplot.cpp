#include "xplot.hpp"
#include<iostream>
#include<string>

#define XSCALE(x) ((int)((x)-x0)*xscale+0.5)
#define YSCALE(y) (height -((int)(((y)-y0)*yscale+0.5))-1)
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

/** Destruct of Figure2D class. Deallocates all resources. */
Figure2D::~Figure2D(){
    XUnmapWindow(display, win);
    XDestroyWindow(display, win);
    XFreePixmap(display, db);
    XFreeFont(display, font);
    XFreeGC(display, gcf);
    XFreeGC(display, gcb);
    XCloseDisplay(display);
}

/** Set up the scaling for the coordinate system. */
void Figure2D::set_scale(double u0, double u1, double v0, double v1){
    x0 = u0;
    x1 = u1;
    y0 = v0;
    y1 = v1;
    xscale = ((double)width)/(x1 - x0);
    yscale = ((double)height)/(y1 - y0);
}

void Figure2D::set_scale(const Point2D& minp, const Point2D& maxp){
    x0 = minp.x;
    x1 = maxp.x;
    y0 = minp.y;
    y1 = maxp.y;
    xscale = ((double)width)/(x1 - x0);
    yscale = ((double)height)/(y1 - y0);
}

/** Allocate a color and return its pixel value. */
unsigned long Figure2D::allocate_color(const char *color){
    unsigned long pixel = 0;
    if(iscolor){
        XColor match;
        XColor exact;
        if(XAllocNamedColor(display, cmap, color, &match, &exact)){
            pixel = match.pixel;
        }
    }

    return pixel;
}


/** Set the background cursor to a given color. */
unsigned long Figure2D::set_background_color(const char *color){
    unsigned long pixel = 0;
    if(iscolor){
        XColor match;
        XColor exact;
        if(XAllocNamedColor(display, cmap, color, &match, &exact)){
            XSetForeground(display, gcb, match.pixel);
            pixel = match.pixel;
        }
    }

    return pixel;
}

/** Set foreground cursor to a given color */
unsigned long Figure2D::set_foreground_color(const char *color){
    unsigned long pixel = 0;
    if(iscolor){
        XColor match;
        XColor exact;
        if(XAllocNamedColor(display, cmap, color, &match, &exact)){
            XSetForeground(display, gcf, match.pixel);
            pixel = match.pixel;
        }
    }

    return pixel;
}

/** Set new font. If the new font does not exist, keep the old one. */
void Figure2D::set_font(const char* reqfont){
    XFontStruct *newfont;
    if(newfont = XLoadQueryFont(display, reqfont)){
        XFreeFont(display, font);
        font = newfont;
        XGCValues xgc;
        xgc.font = font->fid;
        XChangeGC(display, gcf, GCFont, &xgc);
    }
}


/** Draw a single point on the display */
void Figure2D::draw_point(double x, double y){
    XDrawPoint(display, db, gcf, XSCALE(x), YSCALE(y));
}

void Figure2D::draw_point(const Point2D& point){
    auto xp = point.x;
    auto yp = point.y;
    XDrawPoint(display, db, gcf, XSCALE(xp), YSCALE(yp));
}



// ****************************************************************
//                  Figure3D Implementation
// ****************************************************************

}// End namespace xplot
