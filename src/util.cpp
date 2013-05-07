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


void myReshape(int w, int h) {
  viewport.w = w;
  viewport.h = h;

  glViewport (0,0,viewport.w,viewport.h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0f,(GLfloat)viewport.w/(GLfloat)viewport.h,0.1f,100.0f);
}

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

void myKybdHndlr(int key, int x, int y){
    if (key == GLUT_KEY_UP)
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
            transY += transAmt;
        else
            rotX -= rotAmount;

    else if (key == GLUT_KEY_DOWN)
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
            transY -= transAmt;
        else
            rotX += rotAmount;

    else if (key == GLUT_KEY_LEFT)
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
            transX -= transAmt;
        else
            rotY -= rotAmount;

    else if (key == GLUT_KEY_RIGHT)
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
            transX += transAmt;
        else
            rotY += rotAmount;
    else if (key == GLUT_KEY_PAGE_UP){
        // if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
            // transX += transAmt;
        // elsedefRadius,defMass
			if (!specialCase){
            defMass *= 2;
            defRadius *= 1.5;}}
    else if (key == GLUT_KEY_PAGE_DOWN){
        // if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
            // transX += transAmt;
        // else
			if (!specialCase){
            defMass /= 2;
            defRadius /= 1.5;}}
    else
        return;

    glutPostRedisplay ();
}

void dump(){
	appendToFile("this","DUMP:\n");
	ostringstream ss2;
	ss2 << "listOfSpheres size: " << listOfSpheres.size() << endl;
	for (int i = 0; i < listOfSpheres.size(); i++){
		sphere s1 = listOfSpheres[i];
        ss2 << "Velocity: " << s1.vel.x << " " <<  s1.vel.y << " " <<  s1.vel.z;
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
	}
	ss2 << "listOfLargeSpheres size: " << listOfLargeSpheres.size() << endl;
		appendToFile("this",ss2.str());
		exit(0);
}
void sparse(string s){
	if (s == "clear")
		listOfSpheres.clear();	
	if (s == "numspheres")
		cout << "There are " << listOfSpheres.size() + listOfLargeSpheres.size() << " currently." << endl;
	if (s == "pause")
		paused = !paused;
	if (s == "dump")
		dump();

	if (s == "gravity")
		downwardGravity = !downwardGravity;
}

void myKybdHndlr(unsigned char key, int x, int y){
	//cout << (int) key << endl;

	if (key == '1'){
		//cout << alisCrack << endl;
		 timeStp*=0.5;
		 cout << "timestep is now " << timeStp << endl;
		//cout << running1 << endl;
		//cout << running2 << endl;
		// cout << counter << endl;
	}
 
	if (key == '2'){
		timeStp*=2;
		cout << "timestep is now " << timeStp << endl;
	}

	if (key == '3'){
		numSpheresPerClick /= 2;
		cout << "numspheresperclick is now " << numSpheresPerClick << endl;
	}
	if (key == '4'){
		
		numSpheresPerClick = max(numSpheresPerClick*2,1);
		cout << "numspheresperclick is now " << numSpheresPerClick << endl;
	}
	if (key == '5'){
		defZ += 0.2;
		cout << "defZ is now " << defZ << endl;
	}
	if (key == '6'){
		defZ -= 0.2;
		cout << "defZ is now " << defZ << endl;
	}
	if ((int) key == 27){  // ESC key.
		infile.close();
		if (saveToFile)
			appendToFile(fname,globalToAppend);
        exit(0);
	}
    else if (key == '+')
        transZ += zoomAmt;

    else if (key == '-')
        transZ -= zoomAmt;
        
    else if ((int) key == 13){
		//cout << "enter pressed. The globalTyping var was: "<< globalTyping << endl;
		sparse(globalTyping);
		globalTyping = "";
	}
    else{
		globalTyping += key;
        return;
        
	}

    glutPostRedisplay ();
}

void myParse(std::string file) {
  std::ifstream inpfile(file.c_str());
  vector<Vect3> tmpFrame;
  if(!inpfile.is_open()) {
    std::cout << "Unable to open file" << std::endl;
  } else {
    std::string line;
    while(inpfile.good()) {
      std::vector<std::string> splitline;
      std::string buf;

      std::getline(inpfile,line);
      std::stringstream ss(line);

      while (ss >> buf) {
        splitline.push_back(buf);
      }
      //Ignore blank lines
      if(splitline.size() == 0) {
        continue;
      }
      //Ignore comments
      if(splitline[0][0] == '#') {
        continue;
      } 
      if(splitline[0]=="subdivs") {
		alisCrack = 1;
        numDivs = atof(splitline[1].c_str());
      } 
      else if(splitline[0]=="numSpheresPerClick") {
        numSpheresPerClick = atof(splitline[1].c_str());
      } 
      else if(splitline[0]=="scene1") {
        pool = 1;
        alisCrack = 0;
      } 
      else if(splitline[0]=="scene4") {
		box = 1;
	  }
      else if(splitline[0]=="scene2") {
        jeromiesScene = 1;
      } 
      else if(splitline[0]=="scene3") {
        jeromieScene2 = 1;
      } 
      else if(splitline[0]=="scene5") {
        gravityOn=1;
        gConst=0.000005;
        alisCrack=0;
        defMass=1;
        bound = 30;
      } 
      else if(splitline[0]=="numcubed") {
        numCubed = atof(splitline[1].c_str());
      } 
      else if (splitline[0]=="dragOn"){
		dragOn = 1;
	  }
	  else if (splitline[0]=="dragoff"){
		dragOn = 0;
	  }
	  else if (splitline[0]=="gravity"){
		gravityOn = atof(splitline[1].c_str());
	  }
	  else if (splitline[0]=="gravityC"){
		gConst = atof(splitline[1].c_str());
	  }
	  else if (splitline[0]=="R"){
		R = atof(splitline[1].c_str());
	  }
	  else if (splitline[0]=="defRadius"){
		  originalRadius = atof(splitline[1].c_str());
		defRadius = originalRadius;
		  
		  }else {
        std::cerr << "Unknown command: " << splitline[0] << std::endl;
      }
    }
    inpfile.close();
  }

}
void myParse2(std::string file) {
  std::ifstream inpfile(file.c_str());
  vector<Vect3> tmpFrame;
  if(!inpfile.is_open()) {
    std::cout << "Unable to open file" << std::endl;
  } else {
    std::string line;
    while(inpfile.good()) {
      std::vector<std::string> splitline;
      std::string buf;

      std::getline(inpfile,line);
      std::stringstream ss(line);

      while (ss >> buf) {
        splitline.push_back(buf);
      }
      //Ignore blank lines
      if(splitline.size() == 0) {
        continue;
      }
      //Ignore comments
      if(splitline[0][0] == '#') {
        continue;
      } 
      if(splitline[0][0]=='v') {
		//float c = 20;  //bunny
		//float c = 0.006;  //elephant
		//float c = 0.05;  //harley radius = 0.005
		//float c = 0.05;  //vase
		//float c = 0.1;  //design1 - set z = 0
		//float c = 0.04; //design2 tradius = 0.01
		//float c = 0.1; //design3 design4, 0.02
		//float c = 0.2; //designs 4-8, 0.02
		//float c = 0.5;  //design 9, 0.2
		//float c = 0.1;  //design 10, 0.02
		float c = 0.06;
		float tRadius = 0.01;
		float tMass = 1.0f;
		listOfSpheres.push_back(sphere(c*Vect3(atof(splitline[1].c_str()),atof(splitline[2].c_str()),atof(splitline[3].c_str())),Vect3(),tRadius,tMass,Vect3(1,0,0)));
      } 
      else if(splitline[0]=="numSpheresPerClick") {
        numSpheresPerClick = atof(splitline[1].c_str());
      } 
      else if(splitline[0]=="scene1") {
        pool = 1;
        alisCrack = 0;
      } 
      else if(splitline[0]=="scene4") {
		box = 1;
	  }
      else if(splitline[0]=="scene2") {
        jeromiesScene = 1;
      } 
      else if(splitline[0]=="scene3") {
        jeromieScene2 = 1;
      } 
      else if(splitline[0]=="scene5") {
        gravityOn=1;
        gConst=0.000005;
        alisCrack=0;
        defMass=1;
        bound = 30;
      } 
      else if(splitline[0]=="numcubed") {
        numCubed = atof(splitline[1].c_str());
      } 
      else if (splitline[0]=="dragOn"){
		dragOn = 1;
	  }
	  else if (splitline[0]=="dragoff"){
		dragOn = 0;
	  }
	  else if (splitline[0]=="gravity"){
		gravityOn = atof(splitline[1].c_str());
	  }
	  else if (splitline[0]=="gravityC"){
		gConst = atof(splitline[1].c_str());
	  }
	  else if (splitline[0]=="R"){
		R = atof(splitline[1].c_str());
	  }
	  else if (splitline[0]=="defRadius"){
		  originalRadius = atof(splitline[1].c_str());
		defRadius = originalRadius;
		  
		  }else {
        std::cerr << "Unknown command: " << splitline[0] << std::endl;
      }
    }
    inpfile.close();
  }

}
