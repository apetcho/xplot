#ifndef _SIMPLE_XPLOT_H
#define _SIMPLE_XPLOT_H
#include<vector>

extern "C"{
#include<X11/X.h>
#include<X11/Xlib.h>
#include<X11/Xutil.h>
}

namespace xplot{
// ---
struct Point2D{
    double x;
    double y;
};


struct Point3D: public Point2D {
    double z;
};


// 2D Plot
class Figure2D{
private:
    /* Display information needed to manage the display */
    Display *display;
    Window win;
    GC gcb;
    GC gcf;
    XFontStruct *font;
    Pixmap db;
    Colormap cmap;
    int depth;
    int iscolor; // bool?
    int width;
    int height;
    double x0, x1;
    double y0, y1;
    double xscale;
    int mousebtn;
    double xmouse, ymouse;

    int process_event();

public:
    Figure2D(int width, int height, const char *title, const char *server=NULL);

    /****/
    void set_scale(double u0, double u1, double v0, double v1);
    void set_scale(Point2D& src, Point2D& dst);
    void set_font(const char *font);

    // ---
    inline void set_linewidth(int lw){
        XSetLineAttributes(display, gcf, lw, LineSolid, CapRound, JoinRound);
    }

    // ---
    unsigned long allocate_color(const char *color);
    unsigned long set_background_color(const char *color);
    inline void set_backgound_color(unsigned long pixel){
        XSetForeground(display, gcf, pixel);
    }

    void draw_point(double x, double y);
    void draw_point(const Point2D& point);
    void draw_points(double points[][2], int n);
    void draw_points(const std::vector<Point2D>& points);
    void draw_blob(double x, double y, double radius);
    void draw_blob(const Point2D& point, double radius);
    void draw_box(double u0, double v0, double u1, double v1);
    void draw_box(const Point2D& llcnr, const Point2D& urcnr);
    void draw_line(double u0, double v0, double u1, double v);
    void draw_line(const Point2D& start, const Point2D& end);
    void put_image(double u0, double v0, double u1, double v1,
        int width, int height, char *data);
    void put_image(const Point2D& llcnr, const Point2D& urcnr, int width,
        int height, char *data);
    void write_string(double x, double y, const char *text);
    void write_string(const Point2D& position, const char *text);
    // ---
    inline void clear(){
        XFillRectangle(display, db, gcb, 0, 0, width, height);
    }

    // ---
    inline void show(){
        XCopyArea(display, db, win, gcf, 0, 0, width, height, 0, 0);
        XFlush(display);
    }

    // ---
    inline void close(){
        XUnmapWindow(display, win);
        XFlush(display);
    }

    // ---
    inline int checkwin(){
        return process_event();
    }

    // ---
    inline int checkbutton(){
        process_event();
        return mousebtn;
    }

    inline int get_position(double& x, double& y){
        process_event();
        x = xmouse;
        y = ymouse;
        return mousebtn;
    }

    inline int get_position(Point2D& point){
        process_event();
        point.x = xmouse;
        point.y = ymouse;
        return mousebtn;
    }

    // ---
    ~Figure2D();

}; // Figure2D


// 3D Plot
class Figure3D : public Figure2D {
private:
    double uorg, vorg, worg;
    double uscl, vscl, wscl;
    double len;
    double m11, m12, m13;
    double m21, m22, m23;
    int rotating;
    double mx, my;
    double oldxyrot, orldphi;

private:
    void xlate(double& x, double& y, double u, double v, double w);
    void xlate(Point2D& point2d, const Point3D& point3d);

public:
    Figure3D(int width, int height, const char *title, const char* server=NULL);

    // ---
    void set_scale(double u0, double u1, double v0, double v1,
        double w0, double w1);
    void set_scale(Point3D& src, Point3D& dst);

    // ---
    void set_eye(double xyrot, double phi);
    int rotate_mouse();

    // ---
    void draw_point(double u, double v, double w);
    void draw_point(const Point3D& point);
    // ---
    void draw_points(double points[][3], int n);
    void draw_points(const std::vector<Point3D>& point3);

    // ---
    void draw_blob(double u, double v, double w, double radius);
    void draw_blob(const Point3D& point, double radius);

    // ---
    void draw_box(
        double u0, double v0, double w0,
        double u1, double v1, double w1
    );
    void draw_box(const Point3D& ollcnr, const Point3D& iurcnr);

    // ---
    void draw_line(
        double u0, double v0, double w0,
        double u1, double v1, double w1
    );
    void draw_line(const Point3D& start, const Point3D& end);

    // ---
    void draw_axes(double len);

    // --
    inline ~Figure3D(){}

}; // End Figure3D

} // End name space xplot
#endif
