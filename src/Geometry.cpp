#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>
#include <string>
#include <sstream>

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include "Geometry.h"
#include "ColorAndVector.h"

#define sphAcc 20 //higher number => prettier spheres.
#define thresholdForBounce 0 //higher number => bounces happen sooner.

using namespace std;

extern bool saveToFile;
extern int counter;
extern int prevCounter;

// *****************************
// forward Declaration
// *****************************
void appendToFile(string fnameParam, string toAppend);

// *****************************
// plane Implementation
// *****************************
plane::plane(){
    a = 0;
    b = 0;
    c = 0;
    d = 0;
    n = Vect3();
    pt = Vect3();
}    
plane::plane(float ap, float bp, float cp, float dp){
    a = ap;
    b = bp;
    c = cp;
    d = dp;
    n = Vect3();
    pt = Vect3();
}

// *****************************
// sphere Implementation
// *****************************
sphere::sphere(){
    pos = Vect3();
    vel = Vect3();
    r = 0;
    m = 1;
}    
sphere::sphere(Vect3 posp, Vect3 velp, float rp){
    pos = posp;
    vel = velp;
    r = rp;
    m = 1;
}
sphere::sphere(Vect3 posp, Vect3 velp, float rp, float massp){
    pos = posp;
    vel = velp;
    r = rp;
    m = massp;
}
void sphere::render(){
    glTranslatef(pos.x,pos.y,pos.z);
    
    // glutWireSphere(r,sphAcc,sphAcc);
    glutSolidSphere(r,sphAcc,sphAcc);
    
    glTranslatef(-pos.x,-pos.y,-pos.z);
    
    if(saveToFile){
        if (counter != prevCounter){
            appendToFile("test1","EOF\n");
            prevCounter = counter;
        }
        string str = "";
        ostringstream ss;
        ss << pos.x << " " << pos.y << " " << pos.z << "\n";
        appendToFile("test1",ss.str());
    }
}
bool sphere::intersect(plane p){
    float sumabcsquared = pow(p.a,2.0f) + pow(p.b,2.0f) + pow(p.c,2.0f);
    float xo = pos.x;
    float yo = pos.y;
    float zo = pos.z;
    float firstTerm = p.a*xo + p.b*yo + p.c*zo + p.d;
    float xc = xo - p.a * (firstTerm) / (sumabcsquared);
    float yc = yo - p.b * (firstTerm) / (sumabcsquared);
    float zc = zo - p.c * (firstTerm) / (sumabcsquared);

    float d = abs(p.a*xo + p.b*yo + p.c*zo + p.d) / pow(sumabcsquared,0.5f);
    
    if (r < d){ return false; }
    else { return true; }
}
bool sphere::intersect(sphere s2){
    float sumOfRadii = r + s2.r;
    float distBetweenRadii = (pos - s2.pos).getNorm();
    float surfaceDistance = distBetweenRadii - sumOfRadii;
    if (surfaceDistance < thresholdForBounce) { return true; }
    else { return false; }
}
void sphere::move(){
    pos = pos + 1*vel;
}
void sphere::drag(){
void sphere::drag(){
    /*float dragCoef = 2;
    float force = 0.5 * vel * vel * dragCoef;
    float acc = (force / m);
    vel = vel - Vect3(acc, acc, acc);*/
    // vel = vel * (1 / vel.getNorm());
    float selfVelX = vel.x;
    float selfVelY = vel.y;
    float selfVelZ = vel.z;
    
    float pi=acos(-1.0E0);
    float dragCoef = 0.47;
    float force = 0.5 * vel * vel * dragCoef * (1 / (pi * r * r));
    float acc = (force / m);
    
    float accX = sqrt((vel.x * vel.x) / (vel * vel)) * force;
    float accY = sqrt((vel.y * vel.y) / (vel * vel)) * force;
    float accZ = sqrt((vel.z * vel.z) / (vel * vel)) * force;
    
    if (vel.x >= 0) {
        vel.x = max(0.0f, vel.x - accX);
    } else {
        vel.x = min(0.0f, vel.x + accX);
    }
    
    if (vel.y >= 0) {
        vel.y = max(0.0f, vel.y - accY);
    } else {
        vel.y = min(0.0f, vel.y + accY);
    }
    
    if (vel.z >= 0) {
        vel.z = max(0.0f, vel.z - accZ);
    } else {
        vel.z = min(0.0f, vel.z + accZ);
    }
    
}

// *****************************
// KDtree Implementation
// *****************************
KDtree::KDtree(void) {
    KDtree::init(Vect3(1, 1, 1), Vect3(-1, -1, -1)); 
}
KDtree::KDtree(Vect3 & upperLeft, Vect3 & lowerRight) {
    KDtree::init(upperLeft, lowerRight);
}
void KDtree::init(Vect3 & upperLeft, Vect3 & lowerRight) {
    this->UL = Vect3(upperLeft.x, upperLeft.y, upperLeft.z);
    this->LR = Vect3(lowerRight.x, lowerRight.y, lowerRight.z);
}
