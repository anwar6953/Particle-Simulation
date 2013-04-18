#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>
#include <string>
#include <sstream>

#include "ColorAndVector.h"

using namespace std;

// *****************************
// Vect3 Implementation
// *****************************
Vect3::Vect3( void ) {
    init( 0.0f, 0.0f, 0.0f );
}
Vect3::Vect3( float x, float y, float z ) {
    this->h = 0;
    init( x, y, z );
}
Vect3::Vect3( float x, float y, float z, float homogeneousParam ) {
    this->h = homogeneousParam;
    init( x, y, z);
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
    string s = "Vect3[ ";
    s = s + sx.str() + " , " + sy.str() + " , " + sz.str() + " ] ";
    return s;
}
float Vect3::getNorm( void) const {
    return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}
void Vect3::homogenize(float homogeneity) {
    this->h = homogeneity;
}

// *****************************
// VectColor Implementation
// *****************************
VectColor::VectColor( void ) {
    init( 0.0f, 0.0f, 0.0f ); //empty color instantiation vector gets black
}
VectColor::VectColor( float r, float g, float b ) {
    init( r, g, b );
}
VectColor::VectColor( float arrColor [] ) {
    init( arrColor[0], arrColor[1], arrColor[2] );
}
void VectColor::init( float r, float g, float b ) {
    this->r = r;
    this->g = g;
    this->b = b;
}
//define VectColor member overloading
VectColor VectColor::operator+( const VectColor& addMe ) const {
    return VectColor( this->r + addMe.r, this->g + addMe.g, this->b + addMe.b );
}
VectColor VectColor::operator-( const VectColor& subMe ) const {
    return VectColor( this->r - subMe.r, this->g - subMe.g, this->b - subMe.b );
}
VectColor VectColor::operator*( float scalar ) const {
    return VectColor( this->r * scalar, this->g * scalar, this->b * scalar );
}
VectColor VectColor::operator*( const VectColor& multMe ) const{
    return VectColor( this->r * multMe.r, this->g * multMe.g, this->b * multMe.b );
}
string VectColor::printMe() const {
    std::ostringstream sr,sg,sb;
    sr << this->r;
    sg << this->g;
    sb << this->b;
    string s = "VectColor( ";
    s = s + "r: " + sr.str() + " , g: " + sg.str() + " , b: " + sb.str() + " ) ";
    return s;
}

// *****************************
// Non-member operator overload Implementations
// *****************************
Vect3 operator*( float scalar, const Vect3& scaleVec ) {
    //return scaleVec * scalar;
    return Vect3( scalar * scaleVec.x, scalar * scaleVec.y, scalar * scaleVec.z );
}
VectColor operator*( float scalar, const VectColor& scaleVec ) {
    return VectColor( scalar * scaleVec.r, scalar * scaleVec.g, scalar * scaleVec.b );
}

// *****************************
// Function Implementations
// *****************************
Vect3 normalize( const Vect3 & vec ) {
    float norm = vec.getNorm();
    if ( norm == 0.0f ) {
	//cout << "we see a norm of 0. this is bad" << endl;
	return Vect3( 0.0f, 0.0f, 0.0f );
    } else {
	return vec * ( 1.0f / norm );
    }
}

