#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>
#include <string>
#include <sstream>
#include <cstdlib>

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
extern string fname;
extern string globalToAppend;
extern int counter;
extern int prevCounter;
extern float timeStp;
extern float R;
extern float dragCoef;

// *****************************
// forward Declaration
// *****************************
/*
void appendToFile(string fnameParam, string toAppend);
void bindLeaf(KDtree * primary, KDtree * secondary, char type);
KDtree * turnHandle(KDtree * current, char sign, char axis);
char flipSign (char sign);
void renderNode(KDtree * node);
void visitEdge(KDtree * current, char sign1, char axis1, char sign2, char axis2, void (*interfaceNode)(KDtree *));
*/
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
	float c1 = ((float)rand())/RAND_MAX;
	float c2 = ((float)rand())/RAND_MAX;
	float c3 = ((float)rand())/RAND_MAX;
	color = Vect3(c1,c2,c3);
}
plane::plane(float ap, float bp, float cp, float dp, Vect3 cl){
    a = ap;
    b = bp;
    c = cp;
    d = dp;
    n = Vect3(a,b,c);
	isRect = 0;
	color = cl;
}
Vect3 floorPt(Vect3 pt){
	float newx;
	float newy;
	float newz;
	if (pt.x < 0)
		newx = -floor(-pt.x);
	else
		newx = floor(pt.x);
	if (pt.y < 0)
		newy = -floor(-pt.y);
	else
		newy = floor(pt.y);
	if (pt.z < 0)
		newz = -floor(-pt.z);
	else
		newz = floor(pt.z);
	return Vect3(newx,newy,newz);
}
plane::plane(Vect3 p1, Vect3 p2, Vect3 p3, Vect3 p4, bool approx){
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
	if (approx){
	apt1 = floorPt(pt1);
	apt2 = floorPt(pt2);
	apt3 = floorPt(pt3);
	apt4 = floorPt(pt4);
	}
	else {
	apt1 = pt1;
	apt2 = pt2;
	apt3 = pt3;
	apt4 = pt4;
	}
    d = -a * pt1.x - b*pt1.y - c*pt1.z;
	center = (pt1 + pt3) * 0.5;
	
	float c1 = ((float)rand())/RAND_MAX;
	float c2 = ((float)rand())/RAND_MAX;
	float c3 = ((float)rand())/RAND_MAX;
	color = Vect3(c1,c2,c3);
}
plane::plane(Vect3 p1, Vect3 p2, Vect3 p3, Vect3 p4, Vect3 cl, bool approx){
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
	if (approx){
	apt1 = floorPt(pt1);
	apt2 = floorPt(pt2);
	apt3 = floorPt(pt3);
	apt4 = floorPt(pt4);
	}
	else {
	apt1 = pt1;
	apt2 = pt2;
	apt3 = pt3;
	apt4 = pt4;
	}
    d = -a * pt1.x - b*pt1.y - c*pt1.z;
	center = (pt1 + pt3) * 0.5;

	color = cl;
}
void plane::render(){
	if (isRect){
		glBegin(GL_QUADS);
	
		glVertex3f(apt1.x,apt1.y,apt1.z);
	
		glVertex3f(apt2.x,apt2.y,apt2.z);
	
		glVertex3f(apt3.x,apt3.y,apt3.z);
	
		glVertex3f(apt4.x,apt4.y,apt4.z);
		
		glNormal3f(n.x,n.y,n.z);
		
		glEnd();
	}
}
int plane::myType() {
    return 1;
}
// *****************************
// sphere Implementation
// *****************************
sphere::sphere() {
    init(Vect3(), Vect3(), 0.0f, 1.0f);
}    
sphere::sphere(Vect3 posp, Vect3 velp, float rp) {
    init(posp, velp, rp, 1.0f);
}
sphere::sphere(Vect3 posp, Vect3 velp, float rp, float massp) {
    init(posp, velp, rp, massp);
}
sphere::sphere(Vect3 posp, Vect3 velp, float rp, Vect3 cl) {
    init(posp, velp, rp, 1.0f, cl);
}
sphere::sphere(Vect3 posp, Vect3 velp, float rp, float massp, Vect3 cl) {
    init(posp, velp, rp, massp, cl);
}
void sphere::init(Vect3 center, Vect3 velocity, float radius, float mass) {
    pos = center;
    vel = velocity;
    r = radius;
    m = mass;
    collideWithIndex = 0;
	float c1 = ((float)rand())/RAND_MAX;
	float c2 = ((float)rand())/RAND_MAX;
	float c3 = ((float)rand())/RAND_MAX;
	color = Vect3(c1,c2,c3);
	// color = Vect3((rand() / RAND_MAX),(rand() / RAND_MAX),(rand() / RAND_MAX));
}
void sphere::init(Vect3 center, Vect3 velocity, float radius, float mass, Vect3 cl) {
    pos = center;
    vel = velocity;
    r = radius;
    m = mass;
    collideWithIndex = 0;
	color = cl;
}
void sphere::render(){
    glTranslatef(pos.x,pos.y,pos.z);
    
    // glutWireSphere(r,sphAcc,sphAcc);
    glutSolidSphere(r,sphAcc,sphAcc);
    
    glTranslatef(-pos.x,-pos.y,-pos.z);
    
    if(saveToFile){
        if (counter != prevCounter){
			globalToAppend += "EOF\n";
            prevCounter = counter;
        }
        ostringstream ss;
        ss << pos.x << " " << pos.y << " " << pos.z << " " << r << " " << color.x << " " << color.y << " " << color.z << "\n";
		globalToAppend += ss.str();
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
    Vect3 intPt = Vect3(xc,yc,zc);
	// glTranslatef(xc,yc,zc);
    // glutSolidSphere(1,sphAcc,sphAcc);
    // glTranslatef(-xc,-yc,-zc);
	
	float dir = normalize(vel)*(normalize(intPt-pos));
	Vect3 cVel = (dir)*(vel);
	if (dir < 0) dir = -1;
	else dir = 1;
	float t = dir * ((d-r) / cVel.getNorm());	
	bool willIntersect = 0;	
	
    if (r < d){ 
		if (t < timeStp && t > 0)
			willIntersect = 1;
		else
			return false; 
	}
    
	if (p.isRect){
		if (r > 0.08){
			if (willIntersect)
				intPt = pos + (t*vel);
				float magnitude = sqrt(r*r-(intPt - pos)*(intPt - pos));
				intPt = intPt + magnitude * normalize(p.center-intPt);
		}
		//
	if (p.isRect){
		if (r > 0.08){
			if (willIntersect)
				intPt = pos + (t*vel);
				float magnitude = sqrt(r*r-(intPt - pos)*(intPt - pos));
				intPt = intPt + magnitude * normalize(p.center-intPt);
		}
		
		// (1 -> 2)
		Vect3 crossB1 = ((p.pt1-p.pt2)^(p.pt4-p.pt2))^(p.pt1-p.pt2);
		// (2 -> 3)
		Vect3 crossB2 = ((p.pt2-p.pt3)^(p.pt1-p.pt3))^(p.pt2-p.pt3);
		// (3 -> 4)
		Vect3 crossB3 = ((p.pt3-p.pt4)^(p.pt2-p.pt4))^(p.pt3-p.pt4);
		// (4 -> 1)
		Vect3 crossB4 = ((p.pt4-p.pt1)^(p.pt3-p.pt1))^(p.pt4-p.pt1);
		
		bool b1 = crossB1 * (intPt-p.pt1)>0;
		bool b2 = crossB2 * (intPt-p.pt3)>0;
		bool b3 = crossB3 * (intPt-p.pt3)>0;
		bool b4 = crossB4 * (intPt-p.pt1)>0;
		
		if (!b1 || !b2 || !b3 || !b4)
			return false;
	}
	}
	// if (t < 1 && t > 0)
		// cout << "changing velocity." << endl;
	float mag = vel.getNorm();
    float d2 = normalize(-1*vel) * (p.n);
    Vect3 normal = p.n;
    if (d2 < 0){
		normal = normal * -1;
		d2 = normalize(-1*vel)*(p.n*-1);
    }
    vel = (max(1.0*R,1.0)) * mag * normalize(normalize(vel) + 2*d2*(normal));
	
	if (p.isRect && !willIntersect){
		Vect3 diff = pos - Vect3(xc,yc,zc);
		diff = (r-diff.getNorm()) * normalize(diff);
		pos = pos + diff;
	}
	return true;
    //move sphere OUT of plane, if necessary.
}
bool sphere::intersect(sphere s2){
    // float sumOfRadii = (r + s2.r);
    // float distBetweenCenters = (pos - s2.pos).getNorm();
	Vect3 diff = (pos - s2.pos);
    float distBetweenCenters = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
    float sumOfRadii = (r + s2.r)*(r + s2.r);
    float surfaceDistance = distBetweenCenters - sumOfRadii;
    if (surfaceDistance <= thresholdForBounce) { return true; }
    else { return false; }
}
bool sphere::intersect(sphere s2, float rSquared){
    // float sumOfRadii = (r + s2.r);
    // float distBetweenCenters = (pos - s2.pos).getNorm();
	Vect3 diff = (pos - s2.pos);
    float distBetweenCenters = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;
    float surfaceDistance = distBetweenCenters - rSquared;
    if (surfaceDistance <= thresholdForBounce) { return true; }
    else { return false; }
}
void sphere::move() {
    pos = pos + timeStp*vel;
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
    
    float force = 0.5 * vel * vel * dragCoef * (1 / (pi * r * r));
    float acc = (force / m);
    
    float accX = sqrt((vel.x * vel.x) / (vel * vel)) * force;
    float accY = sqrt((vel.y * vel.y) / (vel * vel)) * force;
    float accZ = sqrt((vel.z * vel.z) / (vel * vel)) * force;
    
    if (vel.x >= 0) {
        vel.x = max(0.0f, vel.x - accX*timeStp);
    } else {
        vel.x = min(0.0f, vel.x + accX*timeStp);
    }
    
    if (vel.y >= 0) {
        vel.y = max(0.0f, vel.y - accY*timeStp);
    } else {
        vel.y = min(0.0f, vel.y + accY*timeStp);
    }
    
    if (vel.z >= 0) {
        vel.z = max(0.0f, vel.z - accZ*timeStp);
    } else {
        vel.z = min(0.0f, vel.z + accZ*timeStp);
    }
    
}
int sphere::myType() {
    return 0;
}

// *****************************
// KDtree Implementation
// *****************************
KDtree::KDtree(void) {
    KDtree::init(Vect3(-300, 300, -300), Vect3(300, -300, 300)); 
}
KDtree::KDtree(Vect3 upperLeft, Vect3 lowerRight) {
    KDtree::init(upperLeft, lowerRight);
}
/*
KDtree::~KDtree() {
    if (leftChild != NULL) {
	destroy(leftChild, 1);
	destroy(rightChild, 1);
	//	delete leftChild;
	//delete rightChild;
    }
}
*/
void KDtree::init(Vect3 upperLeft, Vect3 lowerRight) {
    this->UL = Vect3(upperLeft.x, upperLeft.y, upperLeft.z);
    this->LR = Vect3(lowerRight.x, lowerRight.y, lowerRight.z);
    this->isLeaf = false;
    this->axisSplit = 0;
    this->leftChild = NULL;
    this->rightChild = NULL;
    this->nextX = NULL;
    this->prevX = NULL;
    this->nextY = NULL;
    this->prevY = NULL;
    this->nextZ = NULL;
    this->prevZ = NULL;
    this->divX = 1;
    this->divY = 1;
    this->divZ = 1;
    this->leafCount = 0;
}
/*
void KDtree::destroy(KDtree * child, int d) {
    if (child != NULL) {
	destroy(child->rightChild, d + 1);
	destroy(child->leftChild, d + 1);
	cout << string(d, '*') << endl;
	delete child;
	
    }
}
*/
float KDtree::getHypotenuse(void) {
    return (LR - UL).getNorm();
}
void KDtree::constructTree( float baseHypotenuse, char axis, KDtree * rootTree) {
    
    axisSplit = axis;
    if (getHypotenuse() < baseHypotenuse) {
	rootTree->leafCount += 1;
        isLeaf = true;
        return;
    }

    char nextAxis = 0;
    Vect3 rightUL, rightLR;
    Vect3 leftUL , leftLR;

    if ( axis == 'x' ) { // Divide space into 2 along X-axis.
        rightUL = Vect3((UL.x + LR.x)/2.0f, UL.y, UL.z);
        rightLR = LR;
        leftUL = UL;
        leftLR = Vect3((UL.x + LR.x)/2.0f, LR.y, LR.z);
        nextAxis = 'y';
    } else if ( axis == 'y') { // Divide space into 2 along Y-axis. 
        rightUL = Vect3(UL.x, ((UL.y + LR.y)/2.0f), UL.z);
        rightLR = LR;
        leftUL = UL;
        leftLR = Vect3(LR.x, ((UL.y + LR.y)/2.0f), LR.z);
        nextAxis = 'z';
    } else if (axis == 'z') { // Divide space into 2 along Z-axis.
        rightUL = Vect3(UL.x, UL.y, ((UL.z + LR.z)/2.0f));
        rightLR = LR;
        leftUL = UL;
        leftLR = Vect3(LR.x, LR.y, ((UL.z + LR.z)/2.0f));
        nextAxis = 'x';
    }
    rightChild =  new KDtree(rightUL, rightLR);
    leftChild  =  new KDtree(leftUL, leftLR);
    rightChild -> constructTree (baseHypotenuse, nextAxis, rootTree);
    leftChild  -> constructTree (baseHypotenuse, nextAxis, rootTree);
}

void KDtree::getDivisions() {
    KDtree * localTree = this;
    char localAxis = 0;
    while( ! localTree->isLeaf ) {
	localAxis = localTree->axisSplit;
	if      ( localAxis == 'x' ) this->divX *= 2;
	else if ( localAxis == 'y' ) this->divY *= 2;
	else if ( localAxis == 'z' ) this->divZ *= 2;
	localTree = localTree->leftChild;
    }
}

void bindLeaf(KDtree * primary, KDtree * secondary, char type) {
    switch(type) {
    case 'd':
	primary->prevY = secondary;
	secondary->nextY = primary;
	break;
    case 'r':
	primary->nextX = secondary;
	secondary->prevX = primary;
	break;
    case 'b':
	primary->nextZ = secondary;
	secondary->prevZ = primary;
	break;
    }
}

void KDtree::constructWeb() {

    getDivisions();

    float startX(UL.x), startY(UL.y), startZ(UL.z);
    float lenX(abs(LR.x - startX)), lenY(abs(LR.y - startY)), lenZ(abs(LR.z - startZ));
    float x(this->divX * 1.0f), y(this->divY * 1.0f), z(this->divZ * 1.0f);
    float incrX(lenX / x), incrY(lenY / y), incrZ(lenZ / z);

    for (float i = 0.5 * incrX; i < lenX; i += incrX)
	for (float j = 0.5 * incrY; j < lenY; j += incrY)
	    for (float k = 0.5 * incrZ; k < lenZ; k += incrZ) {

		Vect3 posOrig  (startX + i        , startY - j        , startZ + k         );
		Vect3 posBack  (startX + i        , startY - j        , startZ + k + incrZ );
		Vect3 posRight (startX + i + incrX, startY - j        , startZ + k         );
		Vect3 posDown  (startX + i        , startY - j - incrY, startZ + k         );

		KDtree * localOrig = getNode(posOrig);
		KDtree * localBack = getNode(posBack);
		KDtree * localRight = getNode(posRight);
		KDtree * localDown = getNode(posDown);

		if (k < lenZ - incrZ) bindLeaf(localOrig, localBack, 'b');
		if (i < lenX - incrX) bindLeaf(localOrig, localRight, 'r');
		if (j < lenY - incrY) bindLeaf(localOrig, localDown, 'd');
	    } //end loop
}
void KDtree::printMe(int depth) {
    cout << string(depth, '-') << "leaf? " << isLeaf << 
	"; left,right =  " << (leftChild == NULL) << ", " << (rightChild == NULL) <<
	"; upperLeft: " << UL.printMe() << ", lowerRight: " << LR.printMe() << endl;
    if (! isLeaf && ( leftChild != NULL )) {
        leftChild->printMe(depth + 1);
        rightChild->printMe(depth + 1);
    }
}
KDtree * KDtree::getNode(Vect3 & point) {
    KDtree * local = this;
    char localAxis = local->axisSplit;
    while (local->isLeaf == false) {
        // Start by assuming non-corner case
        if (localAxis == 'x') {
            if (point.x < (local->UL.x + local->LR.x) / 2.0f ) local = local->leftChild;
            else local = local->rightChild;
            localAxis = 'y';
        } else if (localAxis == 'y') {
            if (point.y < (local->UL.y + local->LR.y) / 2.0f ) local = local->rightChild;
            else local = local->leftChild;
            localAxis = 'z';
        } else if (localAxis == 'z') {
            if (point.z < (local->UL.z + local->LR.z) / 2.0f ) local = local->leftChild;
            else local = local->rightChild;
            localAxis = 'x';
        }
    } //end while
    return local;
}
void KDtree::render() {

    Vect3 pos = 0.5 * (UL + LR);
    float length = abs(UL.x - LR.x);
    glTranslatef(pos.x,pos.y,pos.z);
    glutSolidSphere(.05,sphAcc,sphAcc);
    glutWireCube( (GLdouble) length);
    glTranslatef(-pos.x,-pos.y,-pos.z);

}

// *****************************
// Sphere / KDtree interfacing utilities
// *****************************
KDtree * turnHandle(KDtree * current, char sign, char axis) {
    KDtree * nextVisit = NULL;
    switch(sign) {

    case '+':
	if      (axis == 'x') nextVisit = current->nextX;
	else if (axis == 'y') nextVisit = current->nextY;
	else if (axis == 'z') nextVisit = current->nextZ;
	break;

    case '-':
	if      (axis == 'x') nextVisit = current->prevX;
	else if (axis == 'y') nextVisit = current->prevY;
	else if (axis == 'z') nextVisit = current->prevZ;
	break;
    }
    return nextVisit;
}

char flipSign (char sign) {
    if (sign == '+') return '-';
    else             return '+';
}

void renderNode(KDtree * node) {
    node->render();
}
/*
struct Navigate {
    char sign;
    char axis;
};
*/
void visitEdge(KDtree * current, char sign1, char axis1, char sign2, char axis2, void (*interfaceNode)(KDtree *)) {
    // TODO: How to capture all 4 of the Nodes.
    KDtree * nextVisit = current;
    (*interfaceNode)(nextVisit);
    // First Step
    nextVisit = turnHandle (nextVisit, sign1, axis1);
    (*interfaceNode)(nextVisit);
    // Second Step
    nextVisit = turnHandle (nextVisit, sign2, axis2);
    (*interfaceNode)(nextVisit);
    // Third Step
    nextVisit = turnHandle (nextVisit, flipSign(sign1), axis1);
    (*interfaceNode)(nextVisit);    
}

void visitCorner(KDtree * current, char sign1, char axis1, char sign2, char axis2, char sign3, char axis3, void (*interfaceNode)(KDtree *)) {
    KDtree * edge1 = current;
    KDtree * edge2 = turnHandle(current, sign3, axis3);

    visitEdge(edge1, sign1, axis1, sign2, axis2, interfaceNode);
    visitEdge(edge2, sign1, axis1, sign2, axis2, interfaceNode);
}

void intersectNode(KDtree * node) {
    
}
