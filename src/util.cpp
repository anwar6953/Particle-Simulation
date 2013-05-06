#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include "Geometry.h"
#include "ColorAndVector.h"
#include "globals.h"

#define sphAcc 20 //higher number => prettier spheres.
#define thresholdForBounce 0 //higher number => bounces happen sooner.

using namespace std;



void myMsHndlr(int button, int state, int x, int y){
switch (button)
  {
  case GLUT_LEFT_BUTTON:
    {
      if (state == 0){prevX = x; prevY = y; return;}
      float r = 1 - 2 *(float) rand() / RAND_MAX;
      float lenOfDrag = sqrt(pow(x-prevX,2.0f)+pow(y-prevY,2.0f)) / 40;
	//43.4
	//65.4
	float amt = 54;
      float xx = (float)(-viewport.w*0.5 + prevX)/amt;
      float yy = (float)(-viewport.h*0.5 + prevY)/-amt;
      Vect3 vel = lenOfDrag * 0.02 * normalize(Vect3(x-prevX,-y+prevY,0)); //MOUSE DRAG decides direction of vel.
	  
	  // cout << xx << " " << yy << endl;
	  for (int i = 0; i < numSpheresPerClick; i++){
	  if (alisCrack && defRadius > originalRadius){
		  if (pool){
		listOfLargeSpheres.push_back(sphere(Vect3(xx,defRadius,-yy),Vect3(vel.x,0,vel.z),defRadius,defMass));}
		else{
		listOfLargeSpheres.push_back(sphere(Vect3(xx,yy,defZ),vel,defRadius,defMass));
	}
		}
	  else{
		  if (pool){
			  listOfSpheres.push_back(sphere(Vect3(xx,defRadius,-yy),Vect3(vel.x,0,vel.z),defRadius,defMass));
		  }
		else{
		listOfSpheres.push_back(sphere(Vect3(xx,yy,defZ),vel,defRadius,defMass));
	}
		}
	  }
	  
	
	
    } break;
  default:
    return;
  }
  glutPostRedisplay ();

}


