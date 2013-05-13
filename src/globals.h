//
//	Geometry.h
#ifndef _globals_h
#define _globals_h

#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <stdlib.h>
#include <cstdlib>
#include "ColorAndVector.h"
#include "Geometry.h"
using namespace std;
extern bool lock;
extern int posInFile;
extern int lastTimeElapsed;
extern FILE * pFile;

class Viewport;

extern int firstTime ;
extern int prevCounter ;
extern int counter ;
extern int fDataCounter ;
extern int lineNo ;
extern ifstream infile;

extern string fname ;
extern bool loadFromFile ;
extern bool saveToFile ;
extern bool openGLrender ;
extern float R;
extern bool dragOn ;
extern bool gravityOn ;
extern float gConst ;
extern bool downwardGravity ;
extern float downwardC ;

extern int dragAlt ;
extern int gravAlt ;
extern int downGrav ;
extern int totalSpheres ;

extern int tmpID;
//extern int convertInt ;

extern bool pool ;
extern bool box ;
extern bool jeromiesScene ;
extern bool jeromieScene2 ;
extern float incrAmtJeromie ;
extern float jx ;
extern float jy ;
extern float jz ;

extern float defZ ;
extern bool removeSpheres ;
extern float bound ;
extern int numSpheresPerClick ;
extern float timeStp ;
extern float defMass ;
extern float originalRadius ;
extern float defRadius ;

extern float dragCoef ;

extern bool specialCase ;
extern float rSqrd ;
extern int numCubed ;
extern float numDivs ;
extern float cameraTracking ;
extern Viewport viewport;

extern KDtree * mainTree;
extern vector<sphere> listOfSpheres;
extern vector<sphere> listOfLargeSpheres;
extern vector<plane> listOfPlanes;

extern int prevX, prevY;
extern bool paused ;

extern time_t initTime;  //for performance uses.
extern vector<vector<Vect3> > fData; //The pre-rendered data from file.
extern string globalToAppend ;
extern string globalTyping ;

extern float scaleAmt ;
extern float rotX ;
extern float rotY ;
extern float transY ;
extern float transX ;
extern float transZ ;
extern float rotAmount ;
extern float transAmt ;
extern float zoomAmt ;

extern float dxL ;
extern float dyL ;
extern float dzL ;
extern float xLookAt ;
extern float yLookAt ;
extern float zLookAt ;
extern bool alisCrack ;


#endif
