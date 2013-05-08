//
//	Geometry.h
#ifndef _Geometry_h
#define _Geometry_h

#include <vector>
#include <list>
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

struct queueNode {
    queueNode * prevKdLink; //these are for KDnodes to link
    queueNode * nextKdLink;
    queueNode * next; //these are for ActiveSpheres
    queueNode * prev;
    int serial;
    queueNode();
};


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
    plane(Vect3, Vect3, Vect3, Vect3, bool);
    plane(float, float, float, float, Vect3);
    plane(Vect3, Vect3, Vect3, Vect3, Vect3, bool);
    Vect3 center;
    Vect3 pt1,pt2,pt3,pt4;
    Vect3 apt1,apt2,apt3,apt4;
    Vect3 n;
    Vect3 color;
    float a,b,c,d;
    float isRect;
    bool approx;
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
    int serial;
    Vect3 pos;
    Vect3 vel;
    Vect3 cumForce; // Force that is acting on sphere in present timestep
    Vect3 momentum; // Momentum of Sphere
    vector<Vect3> localForces;
    float r;
    float m;
    Vect3 color;
    void render();
    bool intersect(sphere);
    bool intersect(sphere, float);
    bool intersect(plane);
    void move();
    void drag();
    int myType();
    void copy(sphere);
    KDtree * KDnode;
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
    KDtree(Vect3, Vect3);
    //    ~KDtree();
    Vect3 UL;
    Vect3 LR;
    bool isLeaf;
    int divX;
    int divY;
    int divZ;
    KDtree * nextX;
    KDtree * prevX;
    KDtree * nextY;
    KDtree * prevY;
    KDtree * nextZ;
    KDtree * prevZ;
    char axisSplit;
    KDtree * leftChild;
    KDtree * rightChild;
    int leafCount;
    queueNode * sphereHead; //the linked list implementation
    queueNode * sphereTail;
    vector<sphere *> localSpheres;
    float getHypotenuse();
    void constructTree(float baseHypotenuse, char axis, KDtree * root);
    void getDivisions();
    void constructWeb();
    void printMe(int depth);
    KDtree * getNode(Vect3 & point);
    void clearNodes();
    void render();
    void fullRender();
 private:
    void init(Vect3, Vect3);
    //    void destroy(KDtree * child, int d);
};

// *****************************
// function prototypes
// *****************************
void appendToFile(string fnameParam, string toAppend);
void bindLeaf(KDtree * primary, KDtree * secondary, char type);

/* The following are used in the linked list implementation 
for spheres that belong to a node.
They make use of the queueNode apparatus for quick de-referencing. */
void deLinkSphere(sphere * front, queueNode * tail);
void linkSphere(sphere * front, queueNode * tail);
void appendLinkSphere(KDtree * node, queueNode * tail);
void insertLinkSphereNext(queueNode * currentFront, queueNode * insertNext);
void removeLinkSphere(KDtree *node, queueNode * removeMe);


KDtree * turnHandle(KDtree * current, char sign, char axis);
char flipSign (char sign);
void renderNode(KDtree * node, sphere *);
KDtree * visitEdge(KDtree * current, char sign1, char axis1, char sign2, char axis2, void (*intersectNode)(KDtree *, sphere * ), bool, sphere * );
void visitCorner(KDtree * current, char sign1, char axis1, char sign2, char axis2, char sign3, char axis3, void (*interfaceNode)(KDtree *, sphere * ), sphere * );
void nodeNeighborTest(sphere * sph, void (*intersectNode)(KDtree *, sphere * ));
void recoverNav(bool * , char *, char * );
void doNothing(KDtree * node, sphere *);
sphere * removeSpherePtr(sphere * sph);
void intersectNode(KDtree * node, sphere *);

// *****************************
// as4.cpp function port prototypes
// *****************************
/* need to make collide visible by Geometry.cpp*/
void collide(sphere& s1, sphere& s2);
#endif
