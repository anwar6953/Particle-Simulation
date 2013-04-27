//
//	Geometry.h
#ifndef _Geometry_h
#define _Geometry_h

#include <vector>
#include "ColorAndVector.h"

using namespace std;

// *****************************
// forward Declaration
// *****************************
class Viewport;
class Shape;
class sphere;
class plane;
class KDtree;

// *****************************
// Viewport prototype
// *****************************
class Viewport {
 public:
    int w, h; // width and height
};

// *****************************
// Shape prototype
// *****************************
class Shape {
 public:
    virtual int myType() = 0;
};

// *****************************
// plane prototype
// *****************************
class plane : public Shape {
 public:
    plane();
    plane(float, float, float, float);
    plane(Vect3, Vect3, Vect3, Vect3);
	plane(float, float, float, float, Vect3);
    plane(Vect3, Vect3, Vect3, Vect3, Vect3);
    Vect3 center;
    Vect3 pt1,pt2,pt3,pt4;
    Vect3 n;
	Vect3 color;
    float a,b,c,d;
    float isRect;
    void render();
    int myType();
};

// *****************************
// sphere prototype
// *****************************
class sphere : public Shape {
 public:
    sphere();
    sphere(Vect3, Vect3, float);
    sphere(Vect3, Vect3, float, float);
	sphere(Vect3, Vect3, float, Vect3);
    sphere(Vect3, Vect3, float, float, Vect3);
    Vect3 pos;
    Vect3 vel;
    float r;
    float m;
    int collideWithIndex;
	Vect3 color;
    vector<Shape *> collideWith;
    void render();
    bool intersect(sphere);
    bool intersect(plane);
    void move();
    void drag();
    int myType();
 private:
    void init(Vect3, Vect3, float, float);
	void init(Vect3, Vect3, float, float, Vect3);
};

// *****************************
// KDtree prototype
// *****************************
class KDtree {
 public:
    KDtree();
    KDtree(Vect3 &, Vect3 &);
    Vect3 UL;
    Vect3 LR;
 private:
    void init(Vect3, Vect3);
};

#endif
