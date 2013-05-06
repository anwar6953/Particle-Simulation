#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>
#include <string>
#include <sstream>
#include <cstdlib>
#include "globals.h"

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include "Geometry.h"
#include "ColorAndVector.h"
// #include "as4.cpp"

#define sphAcc 20 //higher number => prettier spheres.
#define thresholdForBounce 0 //higher number => bounces happen sooner.

using namespace std;


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
	float tP = 0;
	float step = 0.1;
	
while(tP <=1){
				
		Vect3 edge1 = p.apt1 + (p.apt2-p.apt1)*tP;
		Vect3 edge2 = p.apt2 + (p.apt3-p.apt2)*tP;
		Vect3 edge3 = p.apt3 + (p.apt4-p.apt3)*tP;
		Vect3 edge4 = p.apt4 + (p.apt1-p.apt4)*tP;
		bool d1 = (pos-edge1).getNorm()<r;
		bool d2 = (pos-edge2).getNorm()<r;
		bool d3 = (pos-edge3).getNorm()<r;
		bool d4 = (pos-edge4).getNorm()<r;
		if (d1 || d2 || d3 || d4){
			cout << "intersecting" << endl;
			sphere tmp = sphere();
			sphere tmp2 = sphere();
			if(d1){
				Vect3 sP = edge1-(edge1-p.center)*0.99;
				tmp = sphere(sP, Vect3(), (sP-edge1).getNorm(), 100000000);
				tmp2 = sphere(pos,vel,r,m);
				collide(tmp, tmp2);
				//pos = tmp2.pos;
				vel = tmp2.vel;
				//vel = -1*vel;
			}else if(d2){
				Vect3 sP = edge2-(edge2-p.center)*0.99;
				tmp = sphere(sP, Vect3(), (sP-edge2).getNorm(), 100000000);
				tmp2 = sphere(pos,vel,r,m);
				collide(tmp, tmp2);
				//pos = tmp2.pos;
				vel = tmp2.vel;
				//vel = -1*vel;
			}else if(d3){
				Vect3 sP = edge3-(edge3-p.center)*0.99;
				tmp = sphere(sP, Vect3(), (sP-edge3).getNorm(), 100000000);
				tmp2 = sphere(pos,vel,r,m);
				collide(tmp, tmp2);
				//pos = tmp2.pos;
				vel = tmp2.vel;
				//vel = -1*vel;
			}else if(d4){
				Vect3 sP = edge4-(edge4-p.center)*0.99;
				tmp = sphere(sP, Vect3(), (sP-edge4).getNorm(), 100000000);
				tmp2 = sphere(pos,vel,r,m);
				collide(tmp, tmp2);
				//pos = tmp2.pos;
				vel = tmp2.vel;
				//vel = -1*vel;
			}

			//vel = -1*vel;
			return true;
		}
		tP += step;
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
		
		if (!b1 || !b2 || !b3 || !b4){
			return false;
		}
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
    vel = (max(1.0*R,1.0))*mag*normalize(normalize(vel) + 2*d2*(normal));


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
void collide(sphere& s1, plane& p1){
    // float mag = s1.vel.getNorm();
    // float d = normalize(-1*s1.vel) * (p1.n);
    // Vect3 normal = p1.n;
    // if (d < 0){
	// normal = normal * -1;
	// d = normalize(-1*s1.vel)*(p1.n*-1);
    // }
    // s1.vel = mag*normalize(normalize(s1.vel) + 2*d*(normal));
    //move sphere OUT of plane, if necessary.
}
void collide(sphere& s1, sphere& s2){
    //**************************************************************
    // Exact handler of 3D collisions. Based off of Thomas Smid's implementation from http://www.plasmaphysics.org.uk/programs/coll3d_cpp.htm
    //**************************************************************
	
    Vect3& p1(s1.pos), p2(s2.pos), v1(s1.vel) , v2(s2.vel);
    float m1(s1.m), m2(s2.m);
    float r1(s1.r), r2(s2.r);
    Vect3 displacement(p2 - p1);
    Vect3 velDiff(v2 - v1);

    float phi2, fvz1r, thetav, phiv, dr, alpha, beta, sbeta, cbeta, dvz2;

    float totalRadius(r1 + r2);
    float massRatio(m2 / m1);

    Vect3 v_cm = (m1 * v1 + m2 * v2) * (1.0f / (m1 + m2));

    //Vect3 centerVect = ((m1 * s1.vel) + (m2 * s2.vel)) * (1 / (m1 + m2));
    float d = displacement.getNorm();
    float v = velDiff.getNorm();

    velDiff = -1 * velDiff;

    float theta2 = (!d) ? 0 : acos( displacement.z / d );
    if (displacement.x == 0 && displacement.y == 0) phi2 = 0;
    else phi2 = atan2( displacement.y, displacement.x );

    float st(sin(theta2));
    float ct(cos(theta2));
    float sp(sin(phi2));
    float cp(cos(phi2));

    Vect3 vel1r;
    vel1r.x = ct * cp * velDiff.x + ct * sp * velDiff.y - st * velDiff.z;
    vel1r.y =    - sp * velDiff.x      + cp * velDiff.y;
    vel1r.z = st * cp * velDiff.x + st * sp * velDiff.y + ct * velDiff.z;

    fvz1r = (!v) ? 0 : vel1r.z / v;
    if (fvz1r > 1) fvz1r = 1; 
    else if (fvz1r < -1) fvz1r = -1;

    thetav = acos(fvz1r);

    if (vel1r.x == 0 && vel1r.y == 0) phiv = 0;
    else phiv = atan2( vel1r.y, vel1r.x );

    dr = d * sin(thetav) / totalRadius;

    alpha = asin(- dr);
    beta = phiv;
    sbeta = sin(beta);
    cbeta = cos(beta);


    float t = (d * cos(thetav) - totalRadius * sqrt(1 - dr * dr) ) / v;

    float a = tan(thetav + alpha);
    dvz2 = 2 * (vel1r.z + a * (cbeta * vel1r.x + sbeta * vel1r.y)) /
	((1 + a * a) * (1 + massRatio));

    Vect3 vel2r = dvz2 * Vect3(a * cbeta, a * sbeta, 1.0);

    vel1r = vel1r - massRatio * vel2r;

    v1.x = ct * cp * vel1r.x - sp * vel1r.y + st * cp * vel1r.z + v2.x;
    v1.y = ct * sp * vel1r.x + cp * vel1r.y + st * sp * vel1r.z + v2.y;
    v1.z =    - st * vel1r.x                     + ct * vel1r.z + v2.z;

    float tmpx(v2.x), tmpy(v2.y), tmpz(v2.z);
    v2.x = ct * cp * vel2r.x - sp * vel2r.y + st * cp * vel2r.z + tmpx;
    v2.y = ct * sp * vel2r.x + cp * vel2r.y + st * sp * vel2r.z + tmpy;
    v2.z =    - st * vel2r.x                     + ct * vel2r.z + tmpz;

    v1 = (v1 - v_cm) * R + v_cm;
    v2 = (v2 - v_cm) * R + v_cm;
	
	if (pool){
		float thr = 0.01;
		if (v1.y < thr && v1.y > -thr) v1.y = 0;
		if (v2.y < thr && v2.y > -thr) v2.y = 0;
	}
	
	//debugging tool:
	// appendToFile("this","collision takign place\n");
	// ostringstream ss2;
	// if (v1.y != 0){
        // ss2 << "x velocity component was " << s1.vel.x
		// << " y velocity component was " << s1.vel.y 
		// << " z velocity component was " << s1.vel.z 
		// << " other sphere's x velocity component was " << s2.vel.x
		// << " other sphere's y velocity component was " << s2.vel.y 
		// << " other sphere's z velocity component was " << s2.vel.z 
		// << " x position component was " << s1.pos.x
		// << " y position component was " << s1.pos.y 
		// << " z position component was " << s1.pos.z 
		// << " other sphere's x position component was " << s2.pos.x
		// << " other sphere's y position component was " << s2.pos.y 
		// << " other sphere's z position component was " << s2.pos.z 
		// << " new y vel is " << v1.y
		// << "\n";
		// appendToFile("this",ss2.str());
		// exit(0);
	// }
	
	s1.vel = v1;
    s2.vel = v2;

    //**************************************************************

    /*Vect3 vel1 = s1.vel;
     Vect3 vel2 = s2.vel;

     s1.vel = (((s2.m)*(vel2 - vel1)) + (s2.m * vel2) + (s1.m * vel1))*(1/(s1.m+s2.m));
     s2.vel = (((s1.m)*(vel1 - vel2)) + (s2.m * vel2) + (s1.m * vel1))*(1/(s1.m+s2.m));

     //s1.vel = ((s1.m-s2.m)*vel1 + 2 * s2.m * vel2)*(1/(s1.m+s2.m));
     //s2.vel = (2*s1.m*vel1 - (s1.m-s2.m)*vel2)*(1/(s1.m+s2.m));

     // s1.vel = (2 * s2.vel)*(1/(s1.m+s2.m));
     // s2.vel = (2 * s1.vel)*(1/(s1.m+s2.m));*/

	 /////////////////////
	 if (s2.r < 1 && s1.r < 1){
	 // if (0){
		Vect3 pos1 = s1.pos;
		Vect3 pos2 = s2.pos;
		float diff = (s1.r+s2.r) - (s2.pos-s1.pos).getNorm();
		float delta = 0;
		if ((s2.pos-s1.pos).getNorm() < ((s1.r+s2.r)))
			delta = 0.0001;
		Vect3 deltaVector = Vect3(delta,delta,delta);
		s1.pos = pos1+(normalize(pos1-pos2+deltaVector))*((diff)/1);
		s2.pos = pos2+(normalize(pos2-pos1-deltaVector))*((diff)/1);
	}

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
