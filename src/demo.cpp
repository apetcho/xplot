#include<iostream>
#include "xplot.hpp"

#define EXIT_BUTTON_NUMBER  2

void plot2d();
void plot3d();

int main(){
    plot2d();
    plot3d();

    return EXIT_SUCCESS;
}

/** 2D XPlot graphics demo. */
void plot2d(){
    using namespace xplot;
    Figure2D figure(500, 500, "A 2D scene");
    figure.set_background_color("green");
    figure.clear();
    figure.set_foreground_color("blue");
    figure.draw_line(0.1, 0.1, 0.9, 0.9);
    figure.set_foreground_color("yellow");
    figure.draw_line(0.9, 0.1, 0.1, 0.9);
    figure.show();

    while(!figure.checkbutton()){ figure.checkwin(); }

    while(!figure.checkbutton()){figure.checkwin();}

    figure.close();
    
}

/** 3D XPlot graphics demo. */
void plot3d(){
    using namespace xplot;
    Figure3D figure(500, 500, "A 3D scene");
    figure.set_background_color("black");
    while(figure.checkbutton() != EXIT_BUTTON_NUMBER){
        figure.clear();
        figure.set_foreground_color("blue");
        figure.draw_box(0.1, 0.1, 0.1, 0.9, 0.9, 0.9);
        figure.set_foreground_color("green");
        figure.draw_line(0.5, 0.5, 0.5, 0.9, 0.5, 0.5);
        figure.set_foreground_color("yellow");
        figure.draw_line(0.5, 0.5, 0.5, 0.5, 0.9, 0.5);
        figure.set_foreground_color("red");
        figure.draw_line(0.5, 0.5, 0.5, 0.5, 0.5, 0.9);
        figure.show();
        figure.rotate_mouse();
    }

    while(figure.checkbutton() == EXIT_BUTTON_NUMBER){
        figure.close();
    }
}
