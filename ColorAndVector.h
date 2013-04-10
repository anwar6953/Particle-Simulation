//
//  ColorAndVector.h
//  
//
//  Created by Chrysolophus on 2/24/13.
//
//

#ifndef _ColorAndVector_h
#define _ColorAndVector_h
//{
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

using namespace std;
//}

class Vect3 {
public:
    float x, y, z, h;
    //constructors
    Vect3( void );
    Vect3( float x, float y, float z );
	Vect3( float x, float y, float z, float );
    Vect3( float arr [] );
    string printMe() const;
    //overload operators
    //void print( const std::ostream& os ); //cout << Vect3
    Vect3 operator+( const Vect3& addMe ) const; //vector addition
    Vect3 operator-( const Vect3& subMe )const ; //vector addition
    float operator*( const Vect3& dotMe ) const;//vector dot product
    Vect3 operator*( float scalar ) const; //vector scalar multiplication
    // Never going to need XOR for Vect3 anyways.
    Vect3 operator^( const Vect3& crossMe) const; //Cross Product
	float getNorm(void) const;
	void homogenize(float);

private:
    static const int dim = 3; //dimension is 3
    void init( float x, float y, float z );//This is used by constructors
};

//define Vect3 Constructors
Vect3::Vect3( void ) {
    init( 0.0f, 0.0f, 0.0f );
}
Vect3::Vect3( float x, float y, float z ) {
	this->h = 0;
    init( x, y, z );
}

Vect3::Vect3( float arr [] ) {
    init( arr[0], arr[1], arr[2] );
}
void Vect3::init( float x, float y, float z ) {
    this->x = x;
    this->y = y;
    this->z = z;
}

//define Vect3 member operator overloading
Vect3 Vect3::operator+( const Vect3& addMe ) const {
    return Vect3( this->x + addMe.x, this->y + addMe.y, this->z + addMe.z );
}

Vect3 Vect3::operator-( const Vect3& subMe ) const {
    return Vect3( this->x - subMe.x, this->y - subMe.y, this->z - subMe.z );
}

float Vect3::operator*( const Vect3& dotMe ) const {
    return this->x * dotMe.x + this->y * dotMe.y + this->z * dotMe.z;
}

Vect3 Vect3::operator*( float scalar ) const {
    return Vect3( this->x * scalar, this->y * scalar, this->z * scalar);
}

// Cross Product
Vect3 Vect3::operator^( const Vect3& crossMe ) const {
    return Vect3(this->y*crossMe.z - this->z*crossMe.y,
                 -(this->x*crossMe.z - this->z*crossMe.x),
                 this->x*crossMe.y - this->y*crossMe.x);
}

string Vect3::printMe() const {
    std::ostringstream sx,sy,sz;
    sx << this->x;
    sy << this->y;
    sz << this->z;
    string s = "Vect3( ";
    s = s + "x: " + sx.str() + " , y: " + sy.str() + " , z: " + sz.str() + " ) ";
    return s;
}

float Vect3::getNorm( void) const {
	return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}
void Vect3::homogenize(float homogeneity) {
	this->h = homogeneity;
}

//define Vect3 non-member operator overloading

/* ostream& operator<<( ostream& os, const Vect3& vect3 ) {
     os << "Vect3( " << vect3.x << ", " << vect3.y << ", " << vect3.z <<" ) ";
     return os;
 
 }*/
 
Vect3 operator*( float scalar, const Vect3& scaleVec ) {
    //return scaleVec * scalar;
    return Vect3( scalar * scaleVec.x, scalar * scaleVec.y, scalar * scaleVec.z );
}




//*************************************
// User-Defined Functions
//*************************************

Vect3 normalize( const Vect3 & vec ) {
    float norm = vec.getNorm();

    if ( norm == 0.0f ) {
        return Vect3( 0.0f, 0.0f, 0.0f );
    } else {
        return vec * ( 1.0f / norm );
    }
}

#endif
