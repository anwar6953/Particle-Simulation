//
//  ColorAndVector.h

#ifndef _ColorAndVector_h
#define _ColorAndVector_h
//{
#include <string>

using namespace std;
//}

// *****************************
// forward Declaration
// *****************************
class Vect3;
class VectColor; 

// *****************************
// Vect3 prototype
// *****************************
class Vect3 {
 public:
    float x, y, z, h;
    //constructors
    Vect3( void );
    Vect3( float x, float y, float z );
    Vect3( float x, float y, float z, float );
    Vect3( float arr [] );
    string printMe() const;
    //overload o	perators
    Vect3 operator+( const Vect3& addMe ) const; //vector addition
    Vect3 operator-( const Vect3& subMe )const ; //vector addition
    float operator*( const Vect3& dotMe ) const;//vector dot product
    Vect3 operator*( float scalar ) const; //vector scalar multiplication
    // Never going to need XOR for	 Vect3 anyways.
    Vect3 operator^( const Vect3& crossMe) const; //Cross Product
    float getNorm(void) const;
    void homogenize(float);
 private:
    static const int dim = 3; //dimension is 3
    void init( float x, float y, float z );//This is used by constructors
};

// *****************************
// VectColor prototype
// *****************************
class VectColor {
 public:
    float r, g, b;
    //constructors
    VectColor( void );
    VectColor( float r, float g, float b );
    VectColor( float arrColor [] );
    //overload operators
    VectColor operator+( const VectColor& addMe ) const; //vector addition
    VectColor operator-( const VectColor& subMe ) const; //vector addition
    VectColor operator*( const VectColor& MultMe ) const;//Color vector product	
    VectColor operator*( float scalar ) const; //Color scalar multiplication
    string printMe() const;  
 private:
    static const int length = 3; //length is 3
    void init( float x, float y, float z );//This is used by constructors
};

// *****************************
// Non-member operator overload prototypes
// *****************************
Vect3 operator*( float scalar, const Vect3& scaleVec );
VectColor operator*( float scalar, const VectColor& scaleVec );

//*************************************
// User-Defined Functions
//*************************************
Vect3 normalize( const Vect3 & vec );

#endif
