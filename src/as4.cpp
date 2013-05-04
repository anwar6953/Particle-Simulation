//{ Includes, Declarations of functions/classes:
#include "ColorAndVector.h"
#include "Geometry.h"
#include <omp.h>
#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <time.h>
#include <math.h>
#include <list>
#include <stdlib.h>
#include <cstdlib>
#include <map>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif


#define PI 3.14159265  // Should be used from mathlib
#define sphAcc 7

inline float sqr(float x) { return x*x; }
using namespace std;
void myDisplay();

int firstTime = 0;
int prevCounter = 0;
int counter = 0;
int fDataCounter = 0;
int lineNo = 0;
		ifstream infile;


//{ SETTINGS:
string fname = "scenes/test1";
bool loadFromFile = 0;
bool saveToFile = 0;
bool openGLrender = 1;
float R;
bool dragOn = 0;
bool gravityOn = 0;
float gConst = 0.00000005;
bool downwardGravity = 0;
float downwardC = 0.0008;

//SCENES
bool pool = 0;
bool box = 0;
bool jeromiesScene = 0;
bool jeromieScene2 = 0;
	float incrAmtJeromie = 0.3;
		float jx = -3.5;
		float jy = 4.5;
		float jz = 1;

float defZ = 0;
//turn on to remove spheres that go out of bounds.
bool removeSpheres = 1;
float bound = 6;


int numSpheresPerClick = 1;
float timeStp = 1;
float defMass = 0.00001;
// float originalRadius = 0.04;
float originalRadius = 0.05;
float defRadius = originalRadius;

float dragCoef = 0.01;
//special case means all the spheres have same radii value.
bool specialCase = 0;
float rSqrd = (2*defRadius)*(2*defRadius);
int numCubed = 0;
	float numDivs = 100;
	float cameraTracking = 0;
//}

// *********** forward declarations
void initSphere(sphere newSphere);

void appendToFile(string fnameParam, string toAppend){
    ofstream outfile;
    const char * fname = fnameParam.c_str();
    outfile.open(fname, ios_base::app);
    outfile << toAppend;
    outfile.close();
}


//{PERFORMANCE STUFF FOLLOWS (works only on Windows).:
// struct timezone2 
// {
  // __int32  tz_minuteswest; /* minutes W of Greenwich */
  // bool  tz_dsttime;     /* type of dst correction */
// };

// struct timeval2 {
// __int32    tv_sec;         /* seconds */
// __int32    tv_usec;        /* microseconds */
// };

// const __int64 DELTA_EPOCH_IN_MICROSECS= 11644473600000000;
// struct timeval timeStart,timeEnd;
// struct timeval2 tvi;
// struct timeval2 tvf;
// struct timezone2 tz;
// struct tm *tm1; 
// time_t time1;

// int gettimeofday(struct timeval2 *tv/*in*/, struct timezone2 *tz/*in*/)
// {
  // FILETIME ft;
  // __int64 tmpres = 0;
  // TIME_ZONE_INFORMATION tz_winapi;
  // int rez=0;

   // ZeroMemory(&ft,sizeof(ft));
   // ZeroMemory(&tz_winapi,sizeof(tz_winapi));

    // GetSystemTimeAsFileTime(&ft);

    // tmpres = ft.dwHighDateTime;
    // tmpres <<= 32;
    // tmpres |= ft.dwLowDateTime;

    // /*converting file time to unix epoch*/
    // tmpres /= 10;  /*convert into microseconds*/
    // tmpres -= DELTA_EPOCH_IN_MICROSECS; 
    // tv->tv_sec = (__int32)(tmpres*0.000001);
    // tv->tv_usec =(tmpres%1000000);


    //_tzset(),don't work properly, so we use GetTimeZoneInformation
    // rez=GetTimeZoneInformation(&tz_winapi);
    // tz->tz_dsttime=(rez==2)?true:false;
    // tz->tz_minuteswest = tz_winapi.Bias + ((rez==2)?tz_winapi.DaylightBias:0);

  // return 0;
// }

// void resetTimer(){
	// gettimeofday(&tvi, &tz);
// }
// int getTimeSince(int i){
	// gettimeofday(&tvf, &tz); // call gettimeofday()
	// int val = ((tvf.tv_sec - tvi.tv_sec) * 1000000 + tvf.tv_usec - tvi.tv_usec);
	// if (i){
		// std::cout << "This effing slow piece of code took "
        // << val
        // << " us to execute."
        // << std::endl;
	// }
	// return val;
// }

/*
#include <sys/timeb.h>
int startTime;
void resetTimer(){
	timeb tb;
	ftime(&tb);
	startTime = tb.millitm + (tb.time & 0xfffff) * 1000;
}
int getMilliCount2(){
	timeb tb;
	ftime(&tb);
	return tb.millitm + (tb.time & 0xfffff) * 1000;
}
void getTimer(){
	int nSpan = getMilliCount2()-startTime;
	if(nSpan < 0)
		nSpan += 0x100000 * 1000;
	printf("Elapsed time = %u milliseconds\n", nSpan);
}
struct timeval tvi1;
struct timeval tvi2;
struct timeval tvi3;
struct timeval tvi4;
struct timeval tvi5;
struct timeval tvf;
struct timezone tz;
int running1 = 0;
int running2 = 0;
int running3 = 0;
int running4 = 0;
int running5 = 0;
void rT(timeval & t){
    gettimeofday(&t, &tz);
    localtime(&t.tv_sec);
}	
int gT(timeval & t){
    gettimeofday(&tvf, &tz);
    localtime(&tvf.tv_sec);
	return (int)(1000000*(tvf.tv_sec - t.tv_sec)+(tvf.tv_usec-t.tv_usec));
}

* */
//} /////PERFORMANCE STUFF ENDS HERE.
//}

//{ Global Variables
//{ Other


Viewport viewport;
GLfloat light_diffuse[] = {1.0, 1.0, 1.0};  /* white diffuse light. */
GLfloat light_ambient[] = {0.1, 0.1, 0.1};  /* white ambient light. */
GLfloat light_specular[] = {1.0, 1.0, 1.0};  /* white specular light. */
GLfloat light_position[] = {3.0, -23.0, 3.0, 0.0};  /* Infinite light location. */
GLfloat light_position2[] = {3.0, 23.0, 3.0, 0.0};  /* Infinite light location. */

KDtree * mainTree = NULL;
sphere * sceneSpheres = NULL;
bool * activeSpheres = NULL;
int numSpheres = 500; //default spheres in a scene
int sceneSpheresIndex = 0; //used to create new spheres (in the place of the sceneSpheres)
vector<sphere> listOfSpheres;
vector<sphere> listOfLargeSpheres;
vector<plane> listOfPlanes;
int prevX, prevY;



bool paused = 0;

time_t initTime;  //for performance uses.
vector<vector<Vect3> > fData; //The pre-rendered data from file.
string globalToAppend = "";
string globalTyping = "";
//}

//{ Defaults for rotations, translations, zooms:
float scaleAmt = 1;
float rotX = 0;
float rotY = 0;
float transY = 0;
float transX = 0;
float transZ = -12;
float rotAmount = 6;  // up, down, left, right
float transAmt = 0.4;  // shift + [the above]
float zoomAmt = 1.6;     // +, -
//Defaults for Lookat:
float dxL = 0;
float dyL = 0;
float dzL = 1;
float xLookAt = 0;
float yLookAt = 0;
float zLookAt = -1;


//}

bool alisCrack = 0;
//}

//{ Functions.
//{ Keyboard/Mouse Handlers,reShape, myParse:
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
      float xx = (float)(-viewport.w*0.5 + prevX)/65.4;
      float yy = (float)(-viewport.h*0.5 + prevY)/-65.4;
      Vect3 vel = lenOfDrag * 0.02 * normalize(Vect3(x-prevX,-y+prevY,0)); //MOUSE DRAG decides direction of vel.
	  
	  // cout << xx << " " << yy << endl;
      //      cout << pool << endl; // its false in default
      for (int i = 0; i < numSpheresPerClick; i++){
	  if (alisCrack && defRadius > originalRadius){
	      if (pool){
		  initSphere(sphere(Vect3(xx,defRadius,-yy),Vect3(vel.x,0,vel.z),defRadius,defMass)); //new
		  listOfLargeSpheres.push_back(sphere(Vect3(xx,defRadius,-yy),Vect3(vel.x,0,vel.z),defRadius,defMass));}
	      else{
		  initSphere(sphere(Vect3(xx,yy,defZ),vel,defRadius,defMass)); //new
		  listOfLargeSpheres.push_back(sphere(Vect3(xx,yy,defZ),vel,defRadius,defMass));
	      }
	  }
	  else{
	      if (pool){
		  initSphere(sphere(Vect3(xx,defRadius,-yy),Vect3(vel.x,0,vel.z),defRadius,defMass)); //new
		  listOfSpheres.push_back(sphere(Vect3(xx,defRadius,-yy),Vect3(vel.x,0,vel.z),defRadius,defMass));
	      }
	      else{
		  initSphere(sphere(Vect3(xx,yy,defZ),vel,defRadius,defMass)); //new //default case
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

void initSphere(sphere newSphere) {
    sceneSpheres[sceneSpheresIndex].copy(newSphere);
    sceneSpheresIndex ++;
    if (sceneSpheresIndex == numSpheres) {
	cout << "allocate more numSpheres, wrapping around to beginning" << endl;
	sceneSpheresIndex = sceneSpheresIndex % numSpheres;
    }
}

void sparse(string s){
	if (s == "clear")
		listOfSpheres.clear();	
	if (s == "numspheres")
		cout << "There are " << listOfSpheres.size() + listOfLargeSpheres.size() << " currently." << endl;
	if (s == "pause")
		paused = !paused;

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
//}

void applyVectorField(sphere & thisSph) {
    //apply vortex field
    // Vect3 fieldVector = thisSph.pos;
    //float denom = pow(fieldVector.x, 2) + pow(fieldVector.y, 2);
    //    fieldVector = Vect3(-1.0f * (1.0f / denom ) * fieldVector.y, (1.0f / denom ) * fieldVector.x, 0.0);
    //    thisSph.vel = thisSph.vel + 0.0001 * fieldVector;

    //fieldVector = Vect3(0,sin(thisSph.pos.x),0);
    //float mag = thisSph.vel.getNorm();
    //    thisSph.vel = mag * normalize(thisSph.vel + 0.01 * fieldVector);
    //    thisSph.vel.y = 0.01 *  (cos(thisSph.pos.y));
    }

void initScene() {
    sceneSpheres = new sphere [numSpheres];
    activeSpheres = new bool [numSpheres];
    for (int i = 0; i < numSpheres; i++) activeSpheres[i] = false; //initialize array to false
    Vect3 UL(-1, 1, -1), LR(1, -1, 1);
    UL = 16 * UL;
    LR = 16 * LR;
    mainTree = new KDtree(UL, LR);
    mainTree->constructTree(2 * sqrt(3) + 0.1, 'x', mainTree);
    mainTree->constructWeb();

    glLineWidth(0.5);
    //glColor3f(1,1,1);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
	
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glShadeModel(GL_SMOOTH);


	float width = 2;
	if (jeromieScene2){
        bool apr2 = 0;
	 dragOn = 1;
 gravityOn = 0;
 downwardGravity = 1;
 downwardC = 0.002;
 originalRadius = 0.05;
 defRadius = originalRadius;
 removeSpheres = 1;
 bound = 6; numSpheresPerClick = 100;
 R=0.6;
        
		float width = 1;
		//plane 3
        listOfPlanes.push_back(plane(Vect3(-1,3,-width), Vect3(1,2,-width), Vect3(1,2,width), Vect3(-1,3,width), apr2));
        
        listOfPlanes.push_back(plane(Vect3(0.1,4,-width), Vect3(0.1,4,width), Vect3(4,5,width), Vect3(4,5,-width), apr2));
        
        listOfPlanes.push_back(plane(Vect3(-0.1,4,-width), Vect3(-0.1,4,width), Vect3(-4,5,width), Vect3(-4,5,-width), apr2));
        
		//5th plane
        listOfPlanes.push_back(plane(Vect3(0.1,0,-width), Vect3(0.1,0,width), Vect3(4,1.5,width), Vect3(4,1.5,-width), apr2));
        
        listOfPlanes.push_back(plane(Vect3(-0.1,-1,-width), Vect3(-0.1,-1,width), Vect3(-4,0.5,width), Vect3(-4,0.5,-width), apr2));
        
        listOfPlanes.push_back(plane(Vect3(-1,-2,-width), Vect3(1,-2,-width), Vect3(1,-3,width), Vect3(-1,-3,width), apr2));
        
   	}
	if(pool){
		float hwidth = 2;
		float hlength = 4;
		float hrails = 0.3;
		float r = 0.2;
		defRadius = r;
		Vect3 sC = Vect3(1,0,0);
		bool apr = 0;
		listOfPlanes.push_back(plane(Vect3(-hlength,0,-hwidth),Vect3(-hlength,0,hwidth),Vect3(hlength,0,hwidth),Vect3(hlength,0,-hwidth),Vect3(0.059,0.330,0.157),apr));
		listOfPlanes.push_back(plane(Vect3(-hlength,0,-hwidth),Vect3(-hlength,0,hwidth),Vect3(-hlength,hrails,hwidth),Vect3(-hlength,hrails,-hwidth),Vect3(0.173,0.094,0.0588),apr));
		listOfPlanes.push_back(plane(Vect3(hlength,0,-hwidth),Vect3(hlength,0,hwidth),Vect3(hlength,hrails,hwidth),Vect3(hlength,hrails,-hwidth),Vect3(0.173,0.094,0.0588),apr));
		listOfPlanes.push_back(plane(Vect3(-hlength,0,-hwidth),Vect3(-hlength,hrails,-hwidth),Vect3(hlength,hrails,-hwidth),Vect3(hlength,0,-hwidth),Vect3(0.173,0.094,0.0588),apr));
		listOfPlanes.push_back(plane(Vect3(-hlength,0,hwidth),Vect3(-hlength,hrails,hwidth),Vect3(hlength,hrails,hwidth),Vect3(hlength,0,hwidth),Vect3(0.173,0.094,0.0588),apr));

		listOfSpheres.push_back(sphere(Vect3(0,r,0),Vect3(.01,0,0),r,defMass,sC));
		listOfSpheres.push_back(sphere(Vect3(1,r,0),Vect3(.01,0,0.01),r,defMass,sC));
		listOfSpheres.push_back(sphere(Vect3(-1,r,0),Vect3(.01,0,-0.01),r,defMass,sC));
		listOfSpheres.push_back(sphere(Vect3(1,r,1),Vect3(.01,0,0),r,defMass,sC));
		// listOfSpheres.push_back(sphere(Vect3(-1,r,0),Vect3(.01,0,-0.01),r,defMass,sC));
		listOfSpheres.push_back(sphere(Vect3(-2,r,1),Vect3(.01,0,0),r,defMass,sC));
		listOfSpheres.push_back(sphere(Vect3(-3,r,0),Vect3(.01,0,-0.01),r,defMass,sC));
		listOfSpheres.push_back(sphere(Vect3(2,r,-1),Vect3(.01,0,0),r,defMass,sC));
	}


	if (jeromiesScene) {
	dragOn = 0;
gravityOn = 0;
downwardGravity = 1;
 downwardC = 0.0008;
 originalRadius = 0.05;
 defRadius = originalRadius;
removeSpheres = 1;
 bound = 6;
numSpheresPerClick = 100;
 R=0.6;
        bool apr2 = 0;
		//lower box depth
		float lBD = -1;
        listOfPlanes.push_back(plane(Vect3(0.2,lBD,0.2),Vect3(0.2,lBD,-0.2),Vect3(2,2,-2), Vect3(2,2,2), apr2));
        
        listOfPlanes.push_back(plane(Vect3(0.2,lBD,-0.2),Vect3(-0.2,lBD,-0.2),Vect3(-2,2,-2), Vect3(2,2,-2), apr2));
        listOfPlanes.push_back(plane(Vect3(-0.2,lBD,-0.2),Vect3(-0.2,lBD,0.2),Vect3(-2,2,2), Vect3(-2,2,-2), apr2));
        
        listOfPlanes.push_back(plane(Vect3(-0.2,lBD,0.2),Vect3(0.2,lBD,0.2),Vect3(2,2,2), Vect3(-2,2,2), apr2));
        
        
       listOfPlanes.push_back(plane(Vect3(2,2.8,-2.8), Vect3(2,2.8,2.8), Vect3(4,4,4), Vect3(4,4,-4), apr2));
        
        listOfPlanes.push_back(plane(Vect3(-2,2.8,-2.8), Vect3(2,2.8,-2.8), Vect3(4,4,-4), Vect3(-4,4,-4), apr2));
        
        listOfPlanes.push_back(plane(Vect3(-2,2.8,2.8), Vect3(-2,2.8,-2.8),  Vect3(-4,4,-4), Vect3(-4,4,4), apr2));
        
        listOfPlanes.push_back(plane(Vect3(2,2.8,2.8), Vect3(-2,2.8,2.8),  Vect3(-4,4,4), Vect3(4,4,4), apr2));
    }


	if (box){
		float overlap = 0.1;
		bool apr2 = 1;
		//left box.
		listOfPlanes.push_back(plane(Vect3(0,1,1 + overlap),  Vect3(0,-1-overlap,1+overlap),  Vect3(0,-1-overlap,-1-overlap),  Vect3(0,1,-1-overlap)       ,apr2));
		//right box.
		listOfPlanes.push_back(plane(Vect3(width,1,1+overlap),Vect3(width,-1-overlap,1+overlap),Vect3(width,-1-overlap,-1-overlap),Vect3(width,1,-1-overlap),apr2));
		//bottom box.
		listOfPlanes.push_back(plane(Vect3(-overlap,-1,1+overlap),Vect3(-overlap,-1,-1-overlap),Vect3(width+overlap,-1,-1-overlap),Vect3(width+overlap,-1,1+overlap),apr2));
		//far box.
		listOfPlanes.push_back(plane(Vect3(-overlap,1,-1),Vect3(-overlap,-1,-1),Vect3(width+overlap,-1,-1),Vect3(width+overlap,1,-1),apr2));
		//near box.
		listOfPlanes.push_back(plane(Vect3(-overlap,1,1),Vect3(-overlap,-1,1),Vect3(width+overlap,-1,1),Vect3(width+overlap,1,1),apr2));
    }
	
    float dist = 0.16;
	if (loadFromFile) numCubed = 0;
    for (int i = 0; i < numCubed; i++) {
        for (int j = 0; j < numCubed; j++){
            for (int k = 0; k < numCubed; k++){
                listOfSpheres.push_back(sphere(Vect3(i*dist+0.2,j*dist,k*dist-0.8),Vect3(0,0,0),0.05,0.00001));

            }
        }
    }
    // listOfSpheres.push_back(sphere(Vect3(-1,0,0),Vect3(0.01,0.005,0),0.2));
    // listOfSpheres.push_back(sphere(Vect3(0,0,0),Vect3(-0.01,0,0),0.2));


}
void removeSphere(int x){
sphere& s = listOfSpheres[x];

	listOfSpheres.erase(listOfSpheres.begin()+x);
}
void removeSphere2(int x){
sphere& s = listOfLargeSpheres[x];


	listOfLargeSpheres.erase(listOfLargeSpheres.begin()+x);
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

void preRender(){
    if (firstTime==1){  time(&initTime); firstTime = 0;}
    if ((counter % 50)==49){
		time_t finalTime;
		time(&finalTime);
		// cout << (float) counter/(finalTime - initTime) << endl;

    }
    counter++;
    fDataCounter++;

    for (int k = 0; k < listOfSpheres.size(); k++) {
        sphere& s1 = listOfSpheres[k];

        s1.move();
        if (dragOn)
            s1.drag();


        // gravity loop
		if (downwardGravity)
			s1.vel.y -= downwardC;
        if (gravityOn)
            for (int j = 0; j < listOfSpheres.size(); j++){
                if (j == k){ continue; }
                sphere& s2 = listOfSpheres[j];
                if ((s2.pos-s1.pos).getNorm() < 0.001) continue;
                // else s1.vel = s1.vel + 0.00000005*(s2.pos-s1.pos)*(s1.m+s2.m)*(1/(s1.m*(s2.pos-s1.pos).getNorm()));
                else s1.vel = s1.vel + 0.00000005*(s2.pos-s1.pos)*(s1.m+s2.m)*(1/(s1.m*(s2.pos-s1.pos).getNorm()));
            }

        //intersection loop
        for (int j = 0; j < listOfSpheres.size(); j++) {
            if (j == k) { continue; }
            sphere& s2 = listOfSpheres[j];
            if (s1.intersect(s2)) {
				collide(s1,s2);
				// s1.move();
			}
        }

        for (int j = 0; j < listOfPlanes.size(); j++) {
            plane p = listOfPlanes[j];
            if (s1.intersect(p)) {
				collide(s1,p);  
				s1.move();
			}
        }
		
		if(saveToFile){
			Vect3 pos = s1.pos;
			if (counter != prevCounter){
				globalToAppend += "EOF\n";
				prevCounter = counter;
			}
			ostringstream ss;
			ss << pos.x << " " << pos.y << " " << pos.z << "\n";
			globalToAppend += ss.str();
		}
    }
	
	
	
}
	struct mykey{
		int a;
		int b;
		int c;
	};		
	bool operator < (const mykey &l, const mykey &r) { 
		return l.a < r.a || (l.a == r.a && l.b < r.b) || (l.a == r.a && l.b == r.b && l.c < r.c ); 
	}
	void jeromiesSphereInit(){
		float r = (float)rand()/RAND_MAX;
		float r2 = (float)rand()/RAND_MAX;
		r -= 0.5;
		r *= 0.5;	
		r2 -= 0.5;
		r2 *= 0.5;
		listOfSpheres.push_back(sphere(Vect3(jx+r,jy+r,jz+r),Vect3(0,-0.02,0),originalRadius,defMass));
		listOfSpheres.push_back(sphere(Vect3(jx+r2,jy+r2,jz+r2),Vect3(0,-0.02,0),originalRadius,defMass));
		if (jx <= -3.5){
			if (jz <= 3.5)
				jz += incrAmtJeromie;
		}
		else if (jx >= 3.5){
			// if (jz >= -3.5)
				jz -= incrAmtJeromie;
		}
		if (jz >= 3.5){
			jz = 3.5;
			jx += incrAmtJeromie;
		} else if (jz <= -3.5){
			jz = -3.5;
			jx-=incrAmtJeromie;
		}
	}
void jeromiesSphereInit2(){

		float r = (float)rand()/RAND_MAX;
		if (r < timeStp){
		listOfSpheres.push_back(sphere(Vect3(-2.5,5.5,0),Vect3(0,0,0),originalRadius,defMass,Vect3(1,0,0)));
		listOfSpheres.push_back(sphere(Vect3(2.5,5.5,0),Vect3(0,0,0),originalRadius,defMass,Vect3(0,0,1)));
	}
	}
void myDisplay() {
	//rT(tvi2);
	

	if (paused)
		return;
    //{ Buffers and Matrices:
    glClear(GL_COLOR_BUFFER_BIT);		    // clear the color buffer
    glClear(GL_DEPTH_BUFFER_BIT);           // clear the depth buffer

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();				        // make sure transformation is "zero'd"
    gluPerspective(60.0f,(GLfloat)viewport.w/(GLfloat)viewport.h,0.1f,100.0f);

	float rradius;
	//following a sphere:
	if (listOfSpheres.size()>1000000){
		xLookAt = listOfSpheres[0].pos.x;
		yLookAt = listOfSpheres[0].pos.y;
		zLookAt = listOfSpheres[0].pos.z;
		// rradius = listOfSpheres[0].r;

		gluLookAt( xLookAt+dxL, yLookAt+dyL, zLookAt+dzL, /* look from camera XYZ */
				xLookAt, yLookAt, zLookAt, /* look at the origin */
				// 0,0, -1, /* look at the origin */
				0, 1, 0); /* positive Y up vector */
	}

    glMatrixMode(GL_MODELVIEW);			    // indicate we are specifying camera transformations
    glLoadIdentity();				        // make sure transformation is "zero'd"

    glTranslatef(transX, transY, transZ);
    glRotatef(rotY,0,1,0);  //rotate about y axis
    glRotatef(rotX,1,0,0);  //} rotate about x



    if (firstTime==1){  time(&initTime); firstTime = 0;}
    if ((counter % 50)==49){
		time_t finalTime;
		time(&finalTime);
		 //cout << (float) counter/(finalTime - initTime) << endl;

    }
    counter++;
    fDataCounter++;

		
	if (jeromiesScene)
		jeromiesSphereInit();
	if (jeromieScene2)
		jeromiesSphereInit2();
	//move all the spheres before doing any calculations.
	for (int k = 0; k < listOfSpheres.size(); k++) {
		//rT(tvi1);
        sphere& s1 = listOfSpheres[k];
        s1.move();
		//running1 += gT(tvi1);
	}
	// vector<vector<vector<vector<int> > > > xMap;

	// vector<vector< vector< vector<int> > > > xMap ( numDivs+1, vector<vector<vector<int> > >(numDivs+1, vector<vector<int> >(numDivs+1, vector<int>(0, 0))));


	map<mykey, vector<int> > zMap;
	map<mykey, vector<int> >::iterator it;
	
	float range = 7;

	if (alisCrack){
		for (int k = 0; k < listOfSpheres.size(); k++) {
			float valx = listOfSpheres[k].pos.x;
			float valy = listOfSpheres[k].pos.y;
			float valz = listOfSpheres[k].pos.z;
			valx *= (numDivs) / (range*2);
			valy *= (numDivs) / (range*2);
			valz *= (numDivs) / (range*2);
			int valfx = floor(valx)+numDivs/2;
			int valfy = floor(valy)+numDivs/2;
			int valfz = floor(valz)+numDivs/2;
			//cout << valfx << " " << valfy << " " << valfz << endl;
			mykey tmp;
			tmp.a = valfx;
			tmp.b = valfy;
			tmp.c = valfz;
			// xMap.at(valfx).at(valfy).at(valfz).push_back(k);
			zMap[tmp].push_back(k);
		}
	}
	//#pragma omp parallel for
    for (int k = 0; k < listOfSpheres.size(); k++) {
        sphere& s1 = listOfSpheres[k];


        if (dragOn)
            s1.drag();


        // gravity loop
		if (downwardGravity)
			s1.vel.y -= downwardC*timeStp;

        if (gravityOn)
            for (int j = 0; j < listOfSpheres.size(); j++){
                if (j == k){ continue; }
                sphere& s2 = listOfSpheres[j];
                if ((s2.pos-s1.pos).getNorm() < 0.001) continue;
                // else s1.vel = s1.vel + 0.00000005*(s2.pos-s1.pos)*(s1.m+s2.m)*(1/(s1.m*(s2.pos-s1.pos).getNorm()));
                else s1.vel = s1.vel + gConst*timeStp*(s2.pos-s1.pos)*(s1.m+s2.m)*(1/(s1.m*(s2.pos-s1.pos).getNorm()));
            }

        //intersection loop
		if (alisCrack){
			float valx = s1.pos.x;
			float valy = s1.pos.y;
			float valz = s1.pos.z;
			valx *= (numDivs) / (range*2);
			valy *= (numDivs) / (range*2);
			valz *= (numDivs) / (range*2);
			int valfx = floor(valx)+numDivs/2;
			int valfy = floor(valy)+numDivs/2;
			int valfz = floor(valz)+numDivs/2;
			
			
			vector<int> thisVector;
			// thisVector = xMap.at(valfx).at(valfy).at(valfz);
			mykey tmp;
			tmp.a = valfx;
			tmp.b = valfy;
			tmp.c = valfz;
			it = zMap.find(tmp);
			if (it != zMap.end()){
				thisVector = it->second;
			}
			for (int i = 0; i < thisVector.size(); i++){
				int thisSph = thisVector[i];
				if (k == thisSph){ continue; }
				sphere& s2 = listOfSpheres[thisSph];
				if (s1.intersect(s2)) {
					collide(s1,s2);
				}
			}
			
		}
		//end of ali's crack.
		else{
			for (int j = 0; j < listOfSpheres.size(); j++) {
				if (j == k) { continue; }
				sphere& s2 = listOfSpheres[j];
				if (specialCase){
					if (s1.intersect(s2,rSqrd)) {
						collide(s1,s2);
						// s1.move();
					}
				}
				else{
					if (s1.intersect(s2)) {
						collide(s1,s2);
						// s1.move();
					}
				}
			}
		}
        for (int j = 0; j < listOfPlanes.size(); j++) {
            plane p = listOfPlanes[j];
            if (s1.intersect(p)) {
				collide(s1,p);  
				s1.move();
			}
        }

		//GLfloat light_spec[] = {s1.color.x, s1.color.y, s1.color.z};
		GLfloat light_spec[] = {1,1,1,1};
		//GLfloat light_amb[] = {s1.color.x, s1.color.y, s1.color.z, 0.01};
		GLfloat light_amb[] = {s1.color.x*0.3,s1.color.y*0.3,s1.color.z*0.3};
		GLfloat light_diff[] = {s1.color.x*0.7, s1.color.y*0.7, s1.color.z*0.7};
		// applyVectorField(s1);
		//glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, light_diff);
		//glMaterialf(GL_FRONT_AND_BACK, GL_DIFFUSE, 1.0);
		//glMaterialf(GL_FRONT_AND_BACK, GL_AMBIENT, 1.0);
		//glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, light_spec);
		//glColor3f(1,1,0);
		//glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, light_amb);
		//glMaterialf(GL_FRONT_AND_BACK, GL_SPECULAR, 10.0);
		//glMaterialf(GL_FRONT_AND_BACK, GL_SPECULAR, 1.0);
		//glColorMaterial(GL_FRONT_AND_BACK, GL_SPECULAR);
		
		//glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
		//glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
		
		//glEnable(GL_COLOR_MATERIAL);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, light_amb);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, light_diff);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, light_spec);
		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 15.0);
		glColor3f(s1.color.x, s1.color.y, s1.color.z);

		//glColor3f(1,1,1);
        s1.render();
		//glDisable(GL_COLOR_MATERIAL);
    }
	for (int i = 0; i < listOfLargeSpheres.size(); i++){  //TESTING large spherse vs small.
		sphere & s1 = listOfLargeSpheres[i];
		s1.move();
		for (int j = 0; j < listOfSpheres.size(); j++){
				sphere & s2 = listOfSpheres[j];
				if (s1.intersect(s2)) {
					collide(s1,s2);
				}
		}
		for (int j = 0; j < listOfLargeSpheres.size(); j++){
				sphere & s2 = listOfLargeSpheres[j];
				if (i == j){ continue; }
				if (s1.intersect(s2)) {
					collide(s1,s2);
				}
		}
		for (int j = 0; j < listOfPlanes.size(); j++) {
            plane p = listOfPlanes[j];
            if (s1.intersect(p)) {
				collide(s1,p);  
				s1.move();
			}
        }
		GLfloat light_spec[] = {1,1,1,1};
		GLfloat light_amb[] = {s1.color.x*0.3,s1.color.y*0.3,s1.color.z*0.3};
		GLfloat light_diff[] = {s1.color.x*0.7, s1.color.y*0.7, s1.color.z*0.7};
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, light_amb);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, light_diff);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, light_spec);
		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 15.0);
		glColor3f(s1.color.x, s1.color.y, s1.color.z);
		s1.render();
	}
	for (int j = 0; j < listOfPlanes.size(); j++) {
        plane p = listOfPlanes[j];
		if (p.isRect){
			GLfloat light_spec[] = {0,0,0,0};
			GLfloat light_amb[] = {p.color.x*0.3, p.color.y*0.3, p.color.z*0.3};
			GLfloat light_diff[] = {p.color.x*0.7, p.color.y*0.7, p.color.z*0.7};
			//GLfloat light_spec[] = {p.color.x, p.color.y, p.color.z};
			//glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
			//glMaterialf(GL_FRONT_AND_BACK, GL_DIFFUSE, 1.0);
			//glMaterialf(GL_FRONT_AND_BACK, GL_AMBIENT, 1.0);
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, light_diff);
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, light_spec);
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, light_amb);
			//glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 1.0);
			//glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);
			//glEnable(GL_COLOR_MATERIAL);
			glColor3f(p.color.x, p.color.y, p.color.z);
			p.render();
			glDisable(GL_COLOR_MATERIAL);
		}
    }

    if(loadFromFile) {
		std::string line;
		while(infile.good()) {
			std::vector<std::string> splitline;
			std::string buf;
		
			std::getline(infile,line);
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
			if(splitline[0][0] == 'E') {
				break;
			}
			if(splitline.size() == 7){
				glColor3f(atof(splitline[4].c_str()),atof(splitline[5].c_str()),atof(splitline[6].c_str()));
				glPushMatrix();
				glTranslatef(atof(splitline[0].c_str()),atof(splitline[1].c_str()),atof(splitline[2].c_str()));
				float r = atof(splitline[3].c_str());
				if (r < 0.3)
					glutSolidSphere(r,7,7);
				else
					glutSolidSphere(r,20,20);
				

				glPopMatrix();
			} else {
				std::cerr << "Unknown command: " << splitline[0] << std::endl;
			}
		}
    }

	if (removeSpheres){
		bool erase;
		for (int j = 0; j < listOfLargeSpheres.size(); j++){
			Vect3 pos = listOfLargeSpheres[j].pos;
			erase = 0;
			if ((pos.x > bound)||(pos.x < -bound)||(pos.y > bound)
			||(pos.y < -bound)||(pos.z > bound)||(pos.z < -bound))
				erase = 1;
			else
				j++;
			if (erase)
				removeSphere2(j);
		}
		for (int j = 0; j < listOfSpheres.size();) {
			Vect3 pos = listOfSpheres[j].pos;
			erase = 0;
			if ((pos.x > bound)||(pos.x < -bound)||(pos.y > bound)
			||(pos.y < -bound)||(pos.z > bound)||(pos.z < -bound))
				erase = 1;
			else
				j++;
			if (erase)
				removeSphere(j);
		}
	}
    glFlush();
    glutSwapBuffers();					// swap buffers (we earlier set double buffer)

	//running2 += gT(tvi2);
}
//}

int main(int argc, char *argv[]) {

	// map<int,vector<int> > xMap;
	
	
    // if (it != xMap.end())
        // it->second.push_back(5);
	// else
		// cout << "no" << endl;	
		
	// xMap[2].push_back(5);
	

	// else
		// cout << "no" << endl;
	// exit(0);
	
	// map<int,string> testMap;
	// string s1 = "first";
	// string s2 = "second";
	// string s3 = "third";
	// string s4 = "fourth";
	// testMap.insert(std::pair<int,string>(2,"first"));
	// testMap.insert(std::pair<int,string>(2,"second"));
	
	

	// exit(0);

	
//{ PARSING:
    // if (argc < 3) { cout << "Please provide the filename and subdivision_Parameter" << endl; exit(0); }
    //string fname = argv[1];
	// string fname = "test1";
    // stpSize = atof(argv[2]);
    // if (argc == 4){
        // if (string(argv[3]) == "-a"){ adaptive = true; }
        // else { adaptive = false; }
    // }
    // if (!adaptive){ //if uniform...
        // float numDivs = floor(1.0f / stpSize);
        // if (numDivs != 1.0f/stpSize){ numDivs++; }
        // stpSize = 1 / numDivs;
    // }
    /*
    if (argc < 2){
		printf("Dude, are you passing something at all? We'll just read input from an arbitrary file for now.\n");
    } else {
        int i = 1;

        while (i < argc) {
            if ((i == 2) && (*(argv[i]) != '-')) {
                fname = argv[2];
                loadFromFile = true;
            }
            if (*(argv[i]) == '-') {
                char *p = argv[i];

                if (*(p+1) == 'g' ) {
                    gravityOn = true;
                } else if (*(p+1) == 'd' ) {
                    dragOn = true;
                } else if (*(p+1) == 'd' ) {
                    saveToFile = true;
                }
            }
            i++;
        }
    }
    */
    myParse(fname);  //}
    if (loadFromFile) {
		infile.open (fname.c_str());
    }

	
//{ Initialization of glut and window:
	if (openGLrender){
    viewport.w = 900;
    viewport.h = 900;
    glutInit(&argc, argv);                        // This initializes glut
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );  // Use a double-buffered window with red, green, and blue channels
    glutInitWindowSize(viewport.w, viewport.h);   //The size and position of the window
    glutInitWindowPosition(0,0);      // x-,y- coords of the topLeft of new window.
    glutCreateWindow(argv[0]);        //} name of window.
	}


//{ initScene() and callBack function:
    initScene();
	if (openGLrender){
    glutDisplayFunc(myDisplay);	    // function to run when its time to draw something
    glutIdleFunc(myDisplay);	    // function to run when its time to draw something
    glutReshapeFunc(myReshape);	    // function to run when the window gets resized
    glutMouseFunc(myMsHndlr);
    glutKeyboardFunc(myKybdHndlr);
    glutSpecialFunc(myKybdHndlr);
    glutMainLoop();				    // infinite loop that will keep drawing and resizing
	}
    else {
		for (int i = 0; i < 10000; i++){
			preRender();
		}
		appendToFile(fname,globalToAppend);
        exit(0);
	}
	return 0;                       //} never reaches here?

}
