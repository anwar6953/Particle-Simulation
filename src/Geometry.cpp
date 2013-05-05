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
#include <windows.h>
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
extern KDtree * mainTree;
extern vector<sphere> listOfSpheres;
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
    float c1 = ((float)rand())/RAND_MAX;
    float c2 = ((float)rand())/RAND_MAX;
    float c3 = ((float)rand())/RAND_MAX;
    color = Vect3(c1,c2,c3);
    init(center, velocity, radius, mass, color);
}
void sphere::init(Vect3 center, Vect3 velocity, float radius, float mass, Vect3 cl) {
    pos = center;
    vel = velocity;
    momentum = mass * velocity;
    cumForce = Vect3(0.0f, 0.0f, 0.0f); // spheres initially have 0 accumulated force
    r = radius;
    m = mass;
    color = cl;
    
    KDnode = mainTree->getNode(center);
    KDnode->localSpheres.push_back(this);
}
void sphere::render(){
    //cout << r << endl;
    //cout << pos.printMe() << endl << endl;
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
void sphere::addForce(Vect3 forceToAdd) {
    localForces.push_back(forceToAdd);
}
void sphere::applyForce(vector<Vect3> & globalForces) {
    // This should happen immediately after move
    for (int i = 0; i < globalForces.size(); i++) {
        // Accumulate all global forces that act on the sphere. This should be a global variable.
        cumForce = cumForce + globalForces.at(i);
    }
    for (int j = 0; j < localForces.size(); j++) {
        // Accumulate all local forces that act on the sphere
        cumForce = cumForce + localForces.at(i);
    }
    // Velocity for next timestep
    vel = vel + (cumForce * (1/m));
    // Momentum for next timestep
    momentum = vel * m;
    // Blank the forces until reapplication for next timestep
    cumForce = Vect3();
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
    pos = pos + timeStp*vel; //update pos
    //cout << "no" << endl;
    sphere * moveMe = removeSpherePtr(this);
    //if (listOfSpheres.size() == 0) return;
    if (moveMe == NULL || listOfSpheres.size() == 0) return;
    KDtree * nextNode = mainTree->getNode(pos);
    nextNode->localSpheres.push_back(moveMe);
    //    KDnode = nextNode;
    nextNode->fullRender();
    //nextNode->fullRender();
    //nodeNeighborTest(this, renderNode);
    //cout << "yes"<< endl;
    char cnt = 0;
    bool wallIntersect [] = {false, false, false, false, false, false}; // [ negX, posX, negY, posY, negZ, posZ ];

    float margin = pos.x;
    bool nextWall = margin < KDnode->UL.x;

    if (KDnode->prevX != NULL) {
	wallIntersect[0] = nextWall;
	cnt += nextWall;
    }

    if (KDnode->nextX != NULL) {    
	nextWall = margin > KDnode->LR.x;
	wallIntersect[1] = nextWall;
	cnt += nextWall;
    }

    if (KDnode->prevY != NULL) {
	margin = pos.y;
	nextWall = margin < KDnode->LR.y;
	wallIntersect[2] = nextWall;
	cnt += nextWall;
    }

    if (KDnode->nextY != NULL) {
	nextWall = margin > KDnode->UL.y;
	wallIntersect[3] = nextWall;
	cnt += nextWall;
    }

    if (KDnode->prevZ != NULL) {
	margin = pos.z;
	nextWall = margin < KDnode->UL.z;
	wallIntersect[4] = nextWall;
	cnt += nextWall;
    }

    if (KDnode->nextZ != NULL) {
	nextWall = margin > KDnode->LR.z;
	wallIntersect[5] = nextWall;
	cnt += nextWall;
    }

    if (cnt == 0 || 1) ; //no need to update node
    else {
	sphere * moveMe = removeSpherePtr(this);
	KDtree * nextNode = NULL;
	char sign [3];
	char axis [3];
	recoverNav(&wallIntersect[0], &sign[0], &axis[0]);
	if (cnt == 1) {
	    nextNode = turnHandle(KDnode, sign[0], axis[0]);
	} else if (cnt == 2) {
	    nextNode = visitEdge(KDnode, sign[0], axis[0], sign[1], axis[1], doNothing, true, moveMe);
	} else if (cnt == 3) {
	    KDtree * edge2 = turnHandle(KDnode, sign[2], axis[2]);
	    nextNode = visitEdge(edge2, sign[0], axis[0], sign[1], axis[1], doNothing, true, moveMe);
	} else cout << "Ve need reinforcements!" << endl;
	//	renderNode(nextNode); //for testing
	nextNode->localSpheres.push_back(moveMe);
	KDnode = nextNode;
    }
    //    KDnode->render();
    //KDnode->fullRender();
    
    
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
void sphere::copy(sphere copyMe) {
    pos.x = copyMe.pos.x;
    pos.y = copyMe.pos.y;
    pos.z = copyMe.pos.z;

    vel.x = copyMe.vel.x;
    vel.y = copyMe.vel.y;
    vel.z = copyMe.vel.z;

    r = copyMe.r;
    m = copyMe.m;

    color.x = copyMe.color.x;
    color.y = copyMe.color.y;
    color.z = copyMe.color.z;
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
    this->sphereHead = NULL;
    this->sphereTail = NULL;
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
void linkSphere(sphere * front, sphere * tail) {
    front->nextKdSphere = tail;
    tail->preKdSphere = front;
}
void appendLinkSphere(KDtree * node, sphere * tail) {
    ; //Todo: implement me
}
void deLinkSphere(sphere * front, sphere * tail) {
    front->nextKdSphere = NULL;
    tail->prevKdSphere = NULL;
}
void insertLinkSphereAfter(sphere * front, sphere * insertNext) {
    linkSphere(insertNext, front->nextKdSphere);
    linkSphere(front, insertNext); //this order matters
}
void removeLinkSphere(sphere * removeMe) {
    if (removeMe->nextKdSphere == NULL) removeMe->prevKdSphere = NULL;
    else deLinkSphere(
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
    //    Vect3 diff = LR - UL;
    //float lenX(abs(diff.x)), lenY(abs(diff.y)), lenZ(abs(diff.z));

    //posZ
    glBegin(GL_LINES);
    glVertex3f(UL.x, UL.y, UL.z);
    glVertex3f(UL.x, LR.y, UL.z);

    glVertex3f(UL.x, LR.y, UL.z);
    glVertex3f(LR.x, LR.y, UL.z);

    glVertex3f(LR.x, LR.y, UL.z);
    glVertex3f(LR.x, UL.y, UL.z);

    glVertex3f(LR.x, UL.y, UL.z);
    glVertex3f(UL.x, UL.y, UL.z);
    glEnd();

    //negZ
    glBegin(GL_LINES);
    glVertex3f(UL.x, UL.y, LR.z);
    glVertex3f(UL.x, LR.y, LR.z);

    glVertex3f(UL.x, LR.y, LR.z);
    glVertex3f(LR.x, LR.y, LR.z);

    glVertex3f(LR.x, LR.y, LR.z);
    glVertex3f(LR.x, UL.y, LR.z);

    glVertex3f(LR.x, UL.y, LR.z);
    glVertex3f(UL.x, UL.y, LR.z);
    glEnd();

    //negX
    glBegin(GL_LINES);
    glVertex3f(UL.x, UL.y, UL.z);
    glVertex3f(UL.x, UL.y, LR.z);

    glVertex3f(UL.x, UL.y, LR.z);
    glVertex3f(UL.x, LR.y, LR.z);

    glVertex3f(UL.x, LR.y, LR.z);
    glVertex3f(UL.x, LR.y, UL.z);

    glVertex3f(UL.x, LR.y, UL.z);
    glVertex3f(UL.x, UL.y, UL.z);
    glEnd();

    //posX
    glBegin(GL_LINES);
    glVertex3f(LR.x, UL.y, UL.z);
    glVertex3f(LR.x, UL.y, LR.z);

    glVertex3f(LR.x, UL.y, LR.z);
    glVertex3f(LR.x, LR.y, LR.z);

    glVertex3f(LR.x, LR.y, LR.z);
    glVertex3f(LR.x, LR.y, UL.z);

    glVertex3f(LR.x, LR.y, UL.z);
    glVertex3f(LR.x, UL.y, UL.z);
    glEnd();

    //posY
    glBegin(GL_LINES);
    glVertex3f(UL.x, UL.y, UL.z);
    glVertex3f(LR.x, UL.y, UL.z);

    glVertex3f(LR.x, UL.y, UL.z);
    glVertex3f(LR.x, UL.y, LR.z);

    glVertex3f(LR.x, UL.y, LR.z);
    glVertex3f(UL.x, UL.y, LR.z);

    glVertex3f(UL.x, UL.y, LR.z);
    glVertex3f(UL.x, UL.y, UL.z);
    glEnd();

    //negY
    glBegin(GL_LINES);
    glVertex3f(UL.x, LR.y, UL.z);
    glVertex3f(LR.x, LR.y, UL.z);

    glVertex3f(LR.x, LR.y, UL.z);
    glVertex3f(LR.x, LR.y, LR.z);

    glVertex3f(LR.x, LR.y, LR.z);
    glVertex3f(UL.x, LR.y, LR.z);

    glVertex3f(UL.x, LR.y, LR.z);
    glVertex3f(UL.x, LR.y, UL.z);
    glEnd();

	/*
    Vect3 pos = 0.5 * (UL + LR);
    float length = abs(UL.x - LR.x);
    glTranslatef(pos.x,pos.y,pos.z);
    glutSolidSphere(.01,sphAcc,sphAcc);
    glutWireCube( (GLdouble) length);
    glTranslatef(-pos.x,-pos.y,-pos.z);
	*/
}
void KDtree::fullRender() {
    Vect3 point = 0.5 * (UL + LR);
    KDtree * local = mainTree;
    char localAxis = local->axisSplit;
    while (local->isLeaf == false) {
	local->render();
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
    
}
void KDtree::clearNodes() {
    localSpheres.clear();
    if (! isLeaf && ( leftChild != NULL )) {
        leftChild->clearNodes();
        rightChild->clearNodes();
    }
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

void renderNode(KDtree * node, sphere * throwAway) {
    node->render();
}
/*
struct Navigate {
    char sign;
    char axis;
};
*/
KDtree * visitEdge(KDtree * current, char sign1, char axis1, char sign2, char axis2, void (*interfaceNode)(KDtree *, sphere *), bool moveCall, sphere * currentSphere) {
    // TODO: How to capture all 4 of the Nodes.
    KDtree * nextVisit = current;
    (*interfaceNode)(nextVisit, currentSphere);
    // First Step
    nextVisit = turnHandle (nextVisit, sign1, axis1);
    (*interfaceNode)(nextVisit, currentSphere);
    // Second Step
    nextVisit = turnHandle (nextVisit, sign2, axis2);
    if (moveCall) return nextVisit;
    (*interfaceNode)(nextVisit, currentSphere);
    // Third Step
    nextVisit = turnHandle (nextVisit, flipSign(sign1), axis1);
    (*interfaceNode)(nextVisit, currentSphere);
}

void visitCorner(KDtree * current, char sign1, char axis1, char sign2, char axis2, char sign3, char axis3, void (*interfaceNode)(KDtree *, sphere *), sphere * currentSphere) {
    KDtree * edge1 = current;
    KDtree * edge2 = turnHandle(current, sign3, axis3);

    visitEdge(edge1, sign1, axis1, sign2, axis2, interfaceNode, false, currentSphere);
    visitEdge(edge2, sign1, axis1, sign2, axis2, interfaceNode, false, currentSphere);
}

void nodeNeighborTest(sphere * sph, void (* interfaceNode)(KDtree *, sphere *) ) {
    char cnt = 0;
    KDtree * myNode = sph->KDnode;
    Vect3 & center = sph->pos;
    float radius = sph->r;
    bool wallIntersect [] = {false, false, false, false, false, false}; // [ negX, posX, negY, posY, negZ, posZ ];

    float margin = center.x - radius;
    bool nextWall = margin < myNode->UL.x;

    if (myNode->prevX != NULL) {
	wallIntersect[0] = nextWall;
	cnt += nextWall;
    }

    if (myNode->nextX != NULL) {    
	margin = center.x + radius;
	nextWall = margin > myNode->LR.x;
	wallIntersect[1] = nextWall;
	cnt += nextWall;
    }

    if (myNode->prevY != NULL) {
	margin = center.y - radius;
	nextWall = margin < myNode->LR.y;
	wallIntersect[2] = nextWall;
	cnt += nextWall;
    }

    if (myNode->nextY != NULL) {
	margin = center.y + radius;
	nextWall = margin > myNode->UL.y;
	wallIntersect[3] = nextWall;
	cnt += nextWall;
    }

    if (myNode->prevZ != NULL) {
	margin = center.z - radius;
	nextWall = margin < myNode->UL.z;
	wallIntersect[4] = nextWall;
	cnt += nextWall;
    }

    if (myNode->nextZ != NULL) {
	margin = center.z + radius;
	nextWall = margin > myNode->LR.z;
	wallIntersect[5] = nextWall;
	cnt += nextWall;
    }
    
    if (cnt == 0) interfaceNode(myNode, sph);
    else {
	char sign [3];
	char axis [3];
	recoverNav(&wallIntersect[0], &sign[0], &axis[0]);
	if (cnt == 1) {
	    interfaceNode(myNode, sph);
	    interfaceNode(turnHandle(myNode, sign[0], axis[0]), sph);
	} else if (cnt == 2) {
	    visitEdge(myNode, sign[0], axis[0], sign[1], axis[1], interfaceNode, false, sph);
	} else if (cnt == 3) {
	    visitCorner(myNode, sign[0], axis[0], sign[1], axis[1], sign[2], axis[2], interfaceNode, sph);
	} else cout << "sound ze alarm!" << endl;
    }
}

void intersectNode(KDtree * node, sphere * currentSphere) {
    //TODO: apply collide over all the spheres in the node. default is to render for debugging.
    sphere * localSphere = NULL;
    for (int i = 0; i < node->localSpheres.size(); i++) {
	localSphere = node->localSpheres[i];
	if (currentSphere == localSphere) continue;
	if (currentSphere->intersect( *localSphere )) {
	    sphere & s1 = *currentSphere;
	    sphere & s2 = *localSphere;
	    collide(s1, s2);//collide them
	}
    }
    //node->render();
}
void doNothing(KDtree * dummyNode, sphere * dummySphere) {
    ;
}
void recoverNav(bool * wallIntersect, char * sign, char * axis) {
    int index = 0;
    for (int i = 0; i < 6; i++ ) {
	if (wallIntersect[i]) {
	    sign[index] = (i % 2 == 0) ? '-' : '+';

	    int j = i / 2;
	    switch(j) {
	    case 0:
		axis[index] = 'x';
		break;
	    case 1:
		axis[index] = 'y';
		break;
	    case 2:
		axis[index] = 'z';
		break;
	    }
	    index++;
	}
    }
}

sphere * removeSpherePtr(sphere * sph) {
    sphere * compareMe = NULL;
    KDtree * oldNode = sph->KDnode;
    sphere * dummy = NULL;

    vector<sphere *> & nodeSpheres = oldNode->localSpheres;

    if (nodeSpheres.size() > 0) dummy = nodeSpheres[0];

    for (int i = 0; i < nodeSpheres.size(); i++) {


	compareMe = nodeSpheres[i];

	if (compareMe == sph) {
	    nodeSpheres.erase(nodeSpheres.begin() + i);
	    break;
	}
    }
    return compareMe;
}
