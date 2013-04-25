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
}    
plane::plane(float ap, float bp, float cp, float dp){
    a = ap;
    b = bp;
    c = cp;
    d = dp;
    n = Vect3(a,b,c);
	isRect = 0;
}
plane::plane(Vect3 p1, Vect3 p2, Vect3 p3, Vect3 p4){
	Vect3 tmpNormal = normalize((p2-p1) ^ (p4-p1));
    a = tmpNormal.x;
    b = tmpNormal.y;
    c = tmpNormal.z;
    n = tmpNormal;
	isRect = 1;
    pt1 = p1;
    pt2 = p2;
    pt3 = p3;
    pt4 = p4;
    d = -a * pt1.x - b*pt1.y - c*pt1.z;
	center = (pt1 + pt3) * 0.5;

}
void plane::render(){
    glBegin(GL_QUADS);

    glVertex3f(pt1.x,pt1.y,pt1.z);

    glVertex3f(pt2.x,pt2.y,pt2.z);

    glVertex3f(pt3.x,pt3.y,pt3.z);

    glVertex3f(pt4.x,pt4.y,pt4.z);
	
    glNormal3f(n.x,n.y,n.z);
    
    glEnd();
}
int plane::myType() {
    return 1;
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
    Vect3 intPt;
    if (r < d){ return false; }
    else { 
		if (p.isRect){
			intPt = Vect3(xc,yc,zc);
			//
			float magnitude = sqrt(r*r-(intPt - pos)*(intPt - pos));
			intPt = intPt + magnitude * normalize(p.center-intPt);
			//
			bool b1 = (p.pt2-p.pt1)*(intPt-p.pt1)>0;
			bool b2 = (p.pt4-p.pt1)*(intPt-p.pt1)>0;
			bool b3 = (p.pt2-p.pt3)*(intPt-p.pt3)>0;
			bool b4 = (p.pt4-p.pt3)*(intPt-p.pt3)>0;
			if (!b1 || !b2 || !b3 || !b4)
				return false;
		}
	}
	float mag = vel.getNorm();
    float d2 = normalize(-1*vel) * (p.n);
    Vect3 normal = p.n;
    if (d2 < 0){
		normal = normal * -1;
		d2 = normalize(-1*vel)*(p.n*-1);
    }
    vel = mag*normalize(normalize(vel) + 2*d2*(normal));
	
	if (p.isRect){
		Vect3 diff = pos - Vect3(xc,yc,zc);
		diff = (r-diff.getNorm()) * normalize(diff);
		pos = pos + diff;
	}
	return true;
    //move sphere OUT of plane, if necessary.
}
bool sphere::intersect(sphere s2){
    float sumOfRadii = r + s2.r;
    float distBetweenRadii = (pos - s2.pos).getNorm();
    float surfaceDistance = distBetweenRadii - sumOfRadii;
    if (surfaceDistance < thresholdForBounce) { return true; }
    else { return false; }
}
void sphere::move() {
    pos = pos + 1*vel;
}
void sphere::drag() {
    /*float dragCoef = 2;
    float force = 0.5 * vel * vel * dragCoef;
    float acc = (force / m);
    vel = vel - Vect3(acc, acc, acc);*/
    // vel = vel * (1 / vel.getNorm());
    float selfVelX = vel.x;
    float selfVelY = vel.y;
    float selfVelZ = vel.z;
    
    float pi=acos(-1.0E0);
    float dragCoef = 0.07;
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
int sphere::myType() {
    return 0;
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
void KDtree::init(Vect3 upperLeft, Vect3 lowerRight) {
    this->UL = Vect3(upperLeft.x, upperLeft.y, upperLeft.z);
    this->LR = Vect3(lowerRight.x, lowerRight.y, lowerRight.z);
}
