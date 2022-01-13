#include "xplot.hpp"
#include<iostream>
#include<cstring>
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

/** Constants use in 3D context specifically */
static const auto MOUSE_BUTTON = 1;
static const auto XYROT_SCALE = -360.0;
static const auto PHI_SCALE = 360.0;


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

/** Draw an array of points on the display */
void Figure2D::draw_points(double points[][2], int n){
    XPoint *xpts = new XPoint[n];
    for(int i=0; i < n; i++){
        xpts[i].x = XSCALE(points[i][0]);
        xpts[i].y = YSCALE(points[i][1]);
    }
    XDrawPoints(display, db, gcf, xpts, n, CoordModeOrigin);
    delete [] xpts;
}

void Figure2D::draw_points(const std::vector<Point2D>& points){
    const auto n = points.size();
    XPoint *xpts = new XPoint[n];
    int i = 0;
    for(auto point: points){
        xpts[i].x = XSCALE(point.x);
        xpts[i].y = YSCALE(point.y);
        i += 1;
    }
    XDrawPoints(display, db, gcf, xpts, n, CoordModeOrigin);
    delete [] xpts;
}

/** Draw a blob of a given radius */
void Figure2D::draw_blob(double x, double y, double radius){
    int xr = (int)(xscale*radius+0.5);
    int yr = (int)(yscale*radius+0.5);
    XDrawPoint(display, db, gcf, XSCALE(x), YSCALE(y));
    XFillArc(display, db, gcf, XSCALE(x)-xr, YSCALE(y)-yr,
        2*xr, 2*yr, 0, 64*360
    );
}

void Figure2D::draw_blob(const Point2D& point, double radius){
    auto xx = point.x;
    auto yy = point.y;
    int xr = (int)(xscale*radius+0.5);
    int yr = (int)(yscale*radius+0.5);
    XDrawPoint(display, db, gcf, XSCALE(xx), YSCALE(yy));
    XFillArc(display, db, gcf, XSCALE(xx)-xr, YSCALE(yy)-yr,
        2*xr, 2*yr, 0, 64*360
    );
}

/** Draw a rectangle with corners at (u0, v0) and (u1, v1) */
void Figure2D::draw_box(double u0, double v0, double u1, double v1){
    draw_line(u0, v0, u0, v1);
    draw_line(u0, v0, u1, v0);
    draw_line(u1, v1, u0, v1);
    draw_line(u1, v1, u1, v0);
}

void Figure2D::draw_box(const Point2D& p0, const Point2D& p1){
    auto u0 = p0.x;
    auto v0 = p0.y;
    auto u1 = p1.x;
    auto v1 = p1.y;
    draw_line(u0, v0, u0, v1);
    draw_line(u0, v0, u1, v0);
    draw_line(u1, v1, u0, v1);
    draw_line(u1, v1, u1, v0);
}

/** Draw a line between two points (u0, v0) and (u1, v1) */
void Figure2D::draw_line(double u0, double v0, double u1, double v1){
    auto xmin = XSCALE(u0);
    auto xmax = XSCALE(u1);
    auto ymin = YSCALE(v0);
    auto ymax = YSCALE(v1);
    XDrawLine(display, db, gcf, xmin, ymin, xmax, ymax);
}

void Figure2D::draw_line(const Point2D& start, const Point2D& end){
    auto xmin = XSCALE(start.x);
    auto xmax = XSCALE(end.x);
    auto ymin = YSCALE(start.y);
    auto ymax = YSCALE(end.y);
    XDrawLine(display, db, gcf, xmin, ymin, xmax, ymax);
}

/** Write a string on the display using the default font. */
void Figure2D::write_string(double x, double y, const char *text){
    auto xpos = XSCALE(x);
    auto ypos = YSCALE(y);
    XDrawString(display, db, gcf, xpos, ypos, text, strlen(text));
}

void Figure2D::write_string(const Point2D& pos, const char *text){
    auto xpos = XSCALE(pos.x);
    auto ypos = YSCALE(pos.y);
    XDrawString(display, db, gcf, xpos, ypos, text, strlen(text));
}


/** Diplay an image at a given coordinates. The image is flipped such that
 * (0, 0) is in the lower left corner.*/
void Figure2D::put_image(
    double u0, double v0, double u1, double v1,
    int imw, int imh, char *data
){
    auto xmin = XSCALE(u0);
    auto xmax = XSCALE(u1);
    auto ymin = YSCALE(v0);
    auto ymax = YSCALE(v1);
    auto xwidth = ABS(xmax - xmin) + 1;
    auto xheight = ABS(ymax - ymin) + 1;
    char *image = new char[xwidth*xheight];

    for(int j=0; j < xheight; j++){
        int dy = (imh - 1 - (j*imh)/xheight)*imw;
        int xj = j * xwidth;
        for(int i=0; i < xwidth; i++){
            int dx = (i*imw)/xwidth;
            image[xj+i] = data[dy+dx];
        }
    }

    Screen *screen = DefaultScreenOfDisplay(display);
    Visual *visual = DefaultVisualOfScreen(screen);
    unsigned depth = DefaultDepthOfScreen(screen);
    XImage *ximg = XCreateImage(
        display, visual, depth, ZPixmap,
        0, image, xwidth, xheight, 8, 0
    );
    int ll = MIN(xmin, xmax);
    int ur = MIN(ymin, ymax);
    XPutImage(display, db, gcf, ximg, 0, 0, ll, ur, xwidth, xheight);
    XDestroyImage(ximg);
    delete [] image;
}

void Figure2D::put_image(
    const Point2D& llcnr, const Point2D& urcnr,
    int imw, int imh, char *data
){
    auto u0 = llcnr.x;
    auto u1 = urcnr.x;
    auto v0 = llcnr.y;
    auto v1 = urcnr.y;
    put_image(u0, v0, u1, v1, imw, imh, data);
}

/** Handle all events which affect the window. */
int Figure2D::process_event(){
    XEvent event;
    XFlush(display);
    XSync(display, false);
    int event_processed = 0;

    while(XCheckWindowEvent(display, win, EVENTMASK, &event)){
        event_processed = 1;
        switch(event.type){
        case ButtonPress:
            mousebtn = event.xbutton.button;
            xmouse = ((double)event.xbutton.x)/xscale+x0;
            ymouse = ((double)event.xbutton.y)/xscale+y0; // XXX yscale?
            break;
        case ButtonRelease:
            mousebtn = NO_MOUSE;
            break;
        case EnterNotify:
            break;
        case LeaveNotify:
            break;
        case MotionNotify:
            xmouse = ((double)event.xbutton.x)/xscale+x0;
            ymouse = ((double)event.xbutton.y)/xscale+y0;
            break;
        case Expose:
            XCopyArea(display, db, win, gcf, 0, 0, width, height, 0, 0);
            XFlush(display);
            break;
        case ConfigureNotify:
            if((event.xconfigure.width != width ||
                event.xconfigure.height != height))
            {
                width = event.xconfigure.width;
                height = event.xconfigure.height;
                XFreePixmap(display, db);
                db = XCreatePixmap(display, win, width, height, depth);
                clear();
                xscale = ((double)width)/(x1 - x0);
                yscale = ((double)height)/(y1 - y0);
            }
            break;
        case MapNotify:
            break;
        case UnmapNotify:
            break;
        case DestroyNotify:
            break;
        default:
            break;
        }// switch
    }

    return event_processed;
}


// ****************************************************************
//                  Figure3D Implementation
// ****************************************************************

/** Figure3D constructor: create the display and set up default coordinate
 * transform. */
Figure3D::Figure3D(int width, int height,
    const char *title, const char*server
) : Figure2D(width, height, title, server), rotating(false)
{
    Figure2D::set_scale(0.0, 1.0, 0.0, 1.0);
    set_scale(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    set_eye(0.0, 0.0);
}

}// End namespace xplot
