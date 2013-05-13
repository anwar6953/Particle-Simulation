//{ Includes, Declarations of functions/classes:
#include "globals.h"
#include "ColorAndVector.h"
#include "Geometry.h"
#include "util.h"
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

#include "../glui/include/GL/glui.h"

#define PI 3.14159265  // Should be used from mathlib
#define sphAcc 7

using namespace std;
void myDisplay();

GLUI *cmd_line_glui=0, *glui;

bool lock = 0;
int posInFile = 0;
int lastTimeElapsed = 0;
FILE * pFile;


int firstTime = 0;
int prevCounter = 0;
int counter = 0;
int fDataCounter = 0;
int lineNo = 0;
ifstream infile;
int main_window;


//{ SETTINGS:
string fname = "scenes/test1";
bool loadFromFile = 0;
bool saveToFile = 0;
bool openGLrender = 1;
float R = 1;
bool dragOn = 0;
int dragAlt = 0;
bool gravityOn = 0;
int gravAlt = 0;
float gConst = 0.00000005; 
bool downwardGravity = 0;
int downGrav = 0;
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


int totalSpheres;
int numSpheresPerClick = 1;
float timeStp = 1;
//float defMass = 0.00001;
 float defMass = 10;
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
float load = 0;

Viewport viewport;
GLfloat light_diffuse[] = {1.0, 1.0, 1.0};  /* white diffuse light. */
GLfloat light_ambient[] = {0.1, 0.1, 0.1};  /* white ambient light. */
// GLfloat light_specular[] = {1.0, 1.0, 1.0};  /* white specular light. */
GLfloat light_specular[] = {0.0,0.0,0.0};  /* white specular light. */
GLfloat light_position[] = {3.0, -23.0, 3.0, 0.0};  /* Infinite light location. */
GLfloat light_position2[] = {3.0, 23.0, 3.0, 0.0};  /* Infinite light location. */

KDtree * mainTree = NULL;
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

bool alisCrack = 1;
//}

//{ Functions.
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

//bool moveable(sphere& s){
	//sphere sCopy = sphere(s.pos, s.vel);
//}

void myGluiIdle( void )
{
    /* According to the GLUT specification, the current window is
     undefined during an idle callback.  So we need to explicitly change
     it if necessary */
    if ( glutGetWindow() != main_window )
        glutSetWindow(main_window);
    totalSpheres = listOfLargeSpheres.size() + listOfSpheres.size();
    glui->sync_live();
    //cout << totalSpheres << endl;
    
    glutPostRedisplay();
}

void initScene() {
// "C:/Users/Gateway/Documents/GitHub/Particle-Simulation/myfile.bin"
	// pFile = fopen ( "myfile.bin" , "ab" );
	pFile = fopen("myfile.bin","rb");
	
	glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
	glClearColor (0.0, 0.0, 0.0, 0.0);
	// if (load){ return; }
	// listOfPlanes.push_back(plane(Vect3(0,1,1 + 0),  Vect3(0,-1-0,1+0),  Vect3(0,-1-0,-1-0),  Vect3(0,1,-1-0)       ,0));
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
	 dragAlt = 1;
 gravAlt = 0;
 downGrav = 1;
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
	dragAlt = 0;
gravAlt = 0;
downGrav = 1;
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
	// if (glutGet(GLUT_ELAPSED_TIME) - lastTimeElapsed < 1000){
		// return;
	// }
	// lastTimeElapsed = glutGet(GLUT_ELAPSED_TIME);
    
	if (paused)
		return;
    //{ Buffers and Matrices:
    glClear(GL_COLOR_BUFFER_BIT);		    // clear the color buffer
    glClear(GL_DEPTH_BUFFER_BIT);           // clear the depth buffer
	if (!load){

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
			mykey tmp;
			tmp.a = valfx;
			tmp.b = valfy;
			tmp.c = valfz;
			zMap[tmp].push_back(k);
		}
	}
    for (int k = 0; k < listOfSpheres.size(); k++) {
        sphere& s1 = listOfSpheres[k];


        if (dragAlt == 1)
            s1.drag();


        // gravity loop
		if (downGrav == 1)
			s1.vel.y -= downwardC*timeStp;

        if (gravAlt == 1)
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
	
	}
	
	
	//****
	// unsigned char* pixels;
	// FILE * pFile;
	// pixels = new unsigned char [4 * viewport.w * viewport.h];
	//****
	//****
	// glReadPixels( 0, 0, viewport.w, viewport.h, GL_RGBA, GL_UNSIGNED_BYTE, &pixels[0] );
    // pFile = fopen ( "myfile.bin" , "wb" );
    // fwrite (pixels , 1 , 4 * viewport.w * viewport.h , pFile );
    // fclose (pFile);
	//****
	//****
	// pFile = fopen ( "myfile.bin" , "r" );
	// fread(pixels,1,4 * viewport.w * viewport.h,pFile);
	// fclose (pFile);
	//****
	//****
	// glDrawPixels(viewport.w, viewport.h, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
	// delete [] pixels;
	//****
	
	float* arr;
	int entries = 4;
	string st ;
	float numS;
	bool lfd = 0;
	bool std = 0;
	
	if(std){
		lock = 1;
		arr = new float [entries*listOfSpheres.size()];
		for (int i = 0; i < listOfSpheres.size(); i++){
			sphere s = listOfSpheres[i];
			arr[i*entries + 0] = s.pos.x;
			arr[i*entries + 1] = s.pos.y;
			arr[i*entries + 2] = s.pos.z;
			arr[i*entries + 3] = s.r;
			cout << s.r << endl;
		}
		float f = listOfSpheres.size();
		if (fwrite(&f, sizeof(float),1,pFile)==0){cout << "yes1" << endl;exit(0);}
		if (fwrite(arr, sizeof(float),entries*listOfSpheres.size(),pFile)==0 && listOfSpheres.size()!=0){cout << "yes2" << endl; exit(0);}
		delete [] arr;
		/************************************************************/
		/***only difference is listOfSpheres vs listOfLargeSpheres***/
		/************(between top and bottom of serialize)***********/
		/************************************************************/
		arr = new float [entries*listOfLargeSpheres.size()];
		for (int i = 0; i < listOfLargeSpheres.size(); i++){
			sphere s = listOfLargeSpheres[i];
			arr[i*entries + 0] = s.pos.x;
			arr[i*entries + 1] = s.pos.y;
			arr[i*entries + 2] = s.pos.z;
			arr[i*entries + 3] = s.r;
			cout << s.r << endl;
		}
		f = listOfLargeSpheres.size();
		if (fwrite(&f, sizeof(float),1,pFile)==0){cout << "yes1" << endl;exit(0);}
		if (fwrite(arr, sizeof(float),entries*listOfLargeSpheres.size(),pFile)==0 && listOfLargeSpheres.size()!=0){cout << "yes2" << endl; exit(0);}
		delete [] arr;
		lock = 0;
	}
	
	if(lfd){
		for (int repeat2x = 0; repeat2x < 2; repeat2x ++){
			// fseek(pFile,posInFile,0);
			int ret = fread(&numS,sizeof(float),1,pFile);
			if (ret){
                int tmp = entries*numS;
				arr = new float[tmp];
				for (int i = 0; i < entries*numS; i++){
					arr[i]=-99999;
				}
				int numRead = fread(arr,sizeof(float),numS*entries,pFile);
				if (numS != 0){
					// cout << "read: " << numRead << endl;
					// cout << "exist: " << numS*entries << endl;
					if (numS*entries != numRead){ cout << "err: " << feof(pFile) << " " << ferror(pFile) << endl; }
					for (int i = 0; i < numS; i++){
						// cout << arr[i*entries+0] << " " << arr[i*entries+1] << " " << arr[i*entries+2] << endl;
						// cout << arr[i*entries+3] << endl;
						sphere(Vect3(arr[i*entries+0],arr[i*entries+1],arr[i*entries+2]),Vect3(),arr[i*entries+3]).render();
					}
				}
				delete [] arr;
				// posInFile += sizeof(float)*(1);
				// posInFile += sizeof(float)*(entries*numS);
			} else { cout << "End of file." << feof(pFile) << endl; /*xit(ret);*/ }
		}
	}

    glFlush();
    glutSwapBuffers();					// swap buffers (we earlier set double buffer)

	//running2 += gT(tvi2);
}
//}

int main(int argc, char *argv[]) {

    // myParse(fname);  //}
    myParse2("scenes/design16");  //}
    if (loadFromFile) {
		infile.open (fname.c_str());
    }

	
//{ Initialization of glut and window:
	if (openGLrender){
    viewport.w = 400;
    viewport.h = 400;
    glutInit(&argc, argv);                        // This initializes glut
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );  // Use a double-buffered window with red, green, and blue channels
    glutInitWindowSize(viewport.w, viewport.h);   //The size and position of the window
    glutInitWindowPosition(0,0);      // x-,y- coords of the topLeft of new window.
    main_window = glutCreateWindow(argv[0]);        //} name of window.
        
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
   
        glui = GLUI_Master.create_glui( "GLUI" );
        (new GLUI_Spinner( glui, "Spheres per Click:", &numSpheresPerClick ))
        ->set_int_limits( 1, 100 );
        (new GLUI_Spinner( glui, "Sphere Radius:", & defRadius ))
        ->set_float_limits( .01f, 5 );
        (new GLUI_Spinner( glui, "Collision Elasticity:", & R ))
        ->set_float_limits( 0.0f, 1 );
        //new GLUI_Checkbox( glui, "Wireframe", &wireframe );
        new GLUI_Checkbox( glui, "Drag?", &dragAlt );
        new GLUI_Checkbox( glui, "Downward Gravity?", &downGrav );
        new GLUI_Checkbox( glui, "Interparticle Gravity?", &gravAlt );
        GLUI_EditText *counter_edittext =
        new GLUI_EditText( glui, "Spheres in Scene:", & totalSpheres);
        counter_edittext->disable();
        new GLUI_Button( glui, "Quit", 0,(GLUI_Update_CB)exit );
        
        glui->set_main_gfx_window( main_window );
        
        /* We register the idle callback with GLUI, *not* with GLUT */
        GLUI_Master.set_glutIdleFunc( myGluiIdle );
        
    glutMainLoop();				    // infinite loop that will keep drawing and resizing
	}
	return 0;                       //} never reaches here?

}
