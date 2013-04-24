//{ Includes, Declarations of functions/classes:
#include "ColorAndVector.h"
#include "Geometry.h"

#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
// #include <unistd.h>
#include <time.h>
#include <math.h>
#include <list>
#include <stdlib.h>
#include <cstdlib> 

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif


#define PI 3.14159265  // Should be used from mathlib
#define sphAcc 20

inline float sqr(float x) { return x*x; }
using namespace std;
void myDisplay();

bool debug = false;
int firstTime = 1;
int prevCounter = 0;
int counter = 0;
int fDataCounter = 0;

//{ SETTINGS:
bool loadFromFile = 0;
bool saveToFile = 0;
bool dragOn = 1;
bool gravityOn = 1;
string fname = "scenes/test1";

float defRadius = 0.2;
float defMass = 1;
//}

void appendToFile(string fnameParam, string toAppend){
    ofstream outfile;
    const char * fname = fnameParam.c_str();
    outfile.open(fname, ios_base::app);
    outfile << toAppend; 
    outfile.close();
}

//}

//{ Global Variables
//{ Other
Viewport viewport;
GLfloat light_diffuse[] = {0.2, 0.1, 1, 1.0};  /* Red diffuse light. */
GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};  /* Infinite light location. */


vector<sphere> listOfSpheres;
vector<plane> listOfPlanes;
int prevX, prevY;

time_t initTime;  //for performance uses.
vector<vector<Vect3> > fData; //The pre-rendered data from file.
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
//}


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
      float xx = (float)(-300 + prevX)/50;
      float yy = (float)(-300 + prevY)/-50;
      Vect3 vel = lenOfDrag * 0.02 * normalize(Vect3(x-prevX,-y+prevY,0)); //MOUSE DRAG decides direction of vel.
      // listOfSpheres.push_back(sphere(Vect3(xx,yy,0),Vect3(0.02*r,0.02*r,0),0.2)); //RANDOM vel dir.
      listOfSpheres.push_back(sphere(Vect3(xx,yy,0),vel,defRadius,defMass));
    } break;
  default:
    return;
  }
  glutPostRedisplay ();

}
void myKybdHndlr(int key, int x, int y){
    // if (key == 27)  //ESC key
		// exit(0);
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
            defMass += 10;
            defRadius += 0.2;}
    else if (key == GLUT_KEY_PAGE_DOWN){
        // if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
            // transX += transAmt;
        // else
            defMass -= 10;
            defRadius -= 0.2;}
    else 
        return;
    
    glutPostRedisplay ();
}
void myKybdHndlr(unsigned char key, int x, int y){ 
	if (key == ' ')  // SPACE key
        exit(0);
	
    else if (key == '+')
        transZ += zoomAmt;
	
    else if (key == '-')
        transZ -= zoomAmt;
    else 
        return;
    
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
      if(splitline[0][0] == 'E') {
        fData.push_back(tmpFrame);
        tmpFrame.clear();
        continue;
      }
      if(splitline.size() == 3){
        tmpFrame.push_back(Vect3(atof(splitline[0].c_str()),atof(splitline[1].c_str()),atof(splitline[2].c_str())));
      } else {
        std::cerr << "Unknown command: " << splitline[0] << std::endl;
      }
    }
    
    inpfile.close();
  }

} 
//}


void applyVectorField(sphere & thisSph) {
    //apply vortex field
    Vect3 fieldVector = thisSph.pos;
    //float denom = pow(fieldVector.x, 2) + pow(fieldVector.y, 2);
    //    fieldVector = Vect3(-1.0f * (1.0f / denom ) * fieldVector.y, (1.0f / denom ) * fieldVector.x, 0.0);
    //    thisSph.vel = thisSph.vel + 0.0001 * fieldVector;

    //fieldVector = Vect3(0,sin(thisSph.pos.x),0);
    //float mag = thisSph.vel.getNorm();
    //    thisSph.vel = mag * normalize(thisSph.vel + 0.01 * fieldVector);
    thisSph.vel.y = 0.01 *  (sin(thisSph.pos.x / 2.0f));
    //    thisSph.vel.y = 0.01 *  (cos(thisSph.pos.y));
    }

void initScene() {
    glLineWidth(0.5);
    glColor3f(1,1,1);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    
    // listOfPlanes.push_back(plane(1,1,1,0));
    
    int numCubed = 0;
    for (int i = 0; i < numCubed; i++) {
        for (int j = 0; j < numCubed; j++){
            for (int k = 0; k < numCubed; k++){
                listOfSpheres.push_back(sphere(Vect3(i,j,k),Vect3(0,0,0),0.2));
            }
        }
    }
    // listOfSpheres.push_back(sphere(Vect3(-1,0,0),Vect3(0.01,0.005,0),0.2));
    // listOfSpheres.push_back(sphere(Vect3(0,0,0),Vect3(-0.01,0,0),0.2));


}

void collide(sphere& s1, sphere& s2){    
    //**************************************************************
    // Exact handler of 3D collisions. Based off of Thomas Smid's implementation from http://www.plasmaphysics.org.uk/programs/coll3d_cpp.htm
    //**************************************************************
    
    /*void collision3D(double R, double m1, double m2, double r1, double r2,
     double& x1, double& y1,double& z1,
     double& x2, double& y2, double& z2,
     double& vx1, double& vy1, double& vz1,
     double& vx2, double& vy2, double& vz2,
     int& error)     {*/

    Vect3& p1(s1.pos), p2(s2.pos), v1(s1.vel) , v2(s2.vel);

    float x1,x2,y1,y2,z1,z2, vx1,vx2,vy1,vy2,vz1,vz2, m1,m2, r1, r2, R;
    
    R = 1.0f;
    
    x1 = p1.x;
    y1 = p1.y;
    z1 = p1.z;
    
    x2 = p2.x;
    y2 = p2.y;
    z2 = p2.z;
    
    m1 = s1.m;
    m2 = s2.m;
    
    r1 = s1.r;
    r2 = s2.r;
    
    vx1 = v1.x;
    vy1 = v1.y;
    vz1 = v1.z;
    
    vx2 = v2.x;
    vy2 = v2.y;
    vz2 = v2.z;

    Vect3 displacement(p2 - p1), velDiff(v2 - v1);

    float  r12,m21,d,v,theta2,phi2,st,ct,sp,cp,vx1r,vy1r,vz1r,fvz1r,
    thetav,phiv,dr,alpha,beta,sbeta,cbeta,t,a,dvz2,
	vx2r,vy2r,vz2r,vx_cm,vy_cm,vz_cm, totalMass;
    
    //     **** initialize some variables ****
    totalMass = m1 + m2;
    r12 = r1 + r2;
    m21 = m2 / m1;
    
    //Vect3 posDiff = s2.pos - p1.pos;
    //Vect3 velDiff = s2.vel - s1.vel;
    Vect3 v_cm = ((m1 * v1) + (m2 * v2)) * (1.0f / totalMass);

    vx_cm = v_cm.x;
    vy_cm = v_cm.y;
    vz_cm = v_cm.z;
    
    //Vect3 centerVect = ((m1 * s1.vel) + (m2 * s2.vel)) * (1 / (m1 + m2));
    //     **** calculate relative distance and relative speed ***
    d = displacement.getNorm();
    v = velDiff.getNorm();
    
    //     **** shift coordinate system so that ball 1 is at the origin ***
    x2 = displacement.x;
    y2 = displacement.y;
    z2 = displacement.z;
    
    //     **** boost coordinate system so that ball 2 is resting ***
    vx1 = - velDiff.x;
    vy1 = - velDiff.y;
    vz1 = - velDiff.z;
    
    //     **** find the polar coordinates of the location of ball 2 ***
    theta2 = (!d) ? 0 : acos( z2 / d );
    phi2 = (x2 == 0 && y2 == 0) ? 0 : atan2( y2, x2 );

    st = sin(theta2);
    ct = cos(theta2);
    sp = sin(phi2);
    cp = cos(phi2);
    
    //     **** express the velocity vector of ball 1 in a rotated coordinate
    //          system where ball 2 lies on the z-axis ******
    vx1r = ct * cp * vx1 + ct * sp * vy1 - st * vz1;
    vy1r = cp * vy1 - sp * vx1;
    vz1r = st * cp * vx1 + st * sp * vy1 + ct * vz1;

    fvz1r = (!v) ? 0 : vz1r / v;
    if (fvz1r > 1) fvz1r = 1;   // fix for possible rounding errors
    else if (fvz1r < -1) fvz1r = -1;

    thetav = acos(fvz1r);
    phiv = (vx1r==0 && vy1r==0) ? 0 : atan2( vy1r, vx1r );
    
    //     **** calculate the normalized impact parameter ***
    dr = d * sin(thetav) / r12;
    
    //     **** calculate impact angles if balls do collide ***
    alpha = asin(- dr);
    beta = phiv;
    sbeta = sin(beta);
    cbeta = cos(beta);
    
    
    //     **** calculate time to collision ***
    t=(d*cos(thetav) -r12*sqrt(1-dr*dr))/v;
    //cout << t << endl;
    //     **** update positions and reverse the coordinate shift ***
    x2=x2+vx2*t +x1;
    y2=y2+vy2*t +y1;
    z2=z2+vz2*t +z1;
    x1=(vx1+vx2)*t +x1;
    y1=(vy1+vy2)*t +y1;
    z1=(vz1+vz2)*t +z1;
    
    
    
    //  ***  update velocities ***
    
    a=tan(thetav+alpha);
    
    dvz2=2*(vz1r+a*(cbeta*vx1r+sbeta*vy1r))/((1+a*a)*(1+m21));
    
    vz2r=dvz2;
    vx2r=a*cbeta*dvz2;
    vy2r=a*sbeta*dvz2;
    vz1r=vz1r-m21*vz2r;
    vx1r=vx1r-m21*vx2r;
    vy1r=vy1r-m21*vy2r;
    
    
    //     **** rotate the velocity vectors back and add the initial velocity
    //           vector of ball 2 to retrieve the original coordinate system ****
    
    vx1=ct*cp*vx1r-sp*vy1r+st*cp*vz1r +vx2;
    vy1=ct*sp*vx1r+cp*vy1r+st*sp*vz1r +vy2;
    vz1=ct*vz1r-st*vx1r               +vz2;
    vx2=ct*cp*vx2r-sp*vy2r+st*cp*vz2r +vx2;
    vy2=ct*sp*vx2r+cp*vy2r+st*sp*vz2r +vy2;
    vz2=ct*vz2r-st*vx2r               +vz2;
    
    
    //     ***  velocity correction for inelastic collisions ***
    
    vx1=(vx1-vx_cm)*R + vx_cm;
    vy1=(vy1-vy_cm)*R + vy_cm;
    vz1=(vz1-vz_cm)*R + vz_cm;
    vx2=(vx2-vx_cm)*R + vx_cm;
    vy2=(vy2-vy_cm)*R + vy_cm;
    vz2=(vz2-vz_cm)*R + vz_cm;
    // if (s1.intersect(s2)){
    s1.vel = Vect3(vx1, vy1, vz1);
    s2.vel = Vect3(vx2, vy2, vz2);
    // }
    // return;
    
    
    //**************************************************************
    
    /*Vect3 vel1 = s1.vel;
     Vect3 vel2 = s2.vel;
     
     s1.vel = (((s2.m)*(vel2 - vel1)) + (s2.m * vel2) + (s1.m * vel1))*(1/(s1.m+s2.m));
     s2.vel = (((s1.m)*(vel1 - vel2)) + (s2.m * vel2) + (s1.m * vel1))*(1/(s1.m+s2.m));
     
     //s1.vel = ((s1.m-s2.m)*vel1 + 2 * s2.m * vel2)*(1/(s1.m+s2.m));
     //s2.vel = (2*s1.m*vel1 - (s1.m-s2.m)*vel2)*(1/(s1.m+s2.m));
     
     // s1.vel = (2 * s2.vel)*(1/(s1.m+s2.m));
     // s2.vel = (2 * s1.vel)*(1/(s1.m+s2.m));*/
     
    Vect3 pos1 = s1.pos;
     Vect3 pos2 = s2.pos;
     float diff = (s1.r+s2.r) - (s2.pos-s1.pos).getNorm();
    float delta = 0.0001;
    Vect3 deltaVector = Vect3(delta,delta,delta);
     s1.pos = pos1+(normalize(pos1-pos2+deltaVector))*((diff)/2);
     s2.pos = pos2+(normalize(pos2-pos1-deltaVector))*((diff)/2);
    
}
void collide(sphere& s1, plane& p1){
    float mag = s1.vel.getNorm();
    float d = normalize(-1*s1.vel) * (p1.n);
    Vect3 normal = p1.n;
    if (d < 0){
	normal = normal * -1;
	d = normalize(-1*s1.vel)*(p1.n*-1);
    }
    s1.vel = mag*normalize(normalize(s1.vel) + 2*d*(normal));
    //move sphere OUT of plane, if necessary.
}

void myDisplay() {
    //{ Buffers and Matrices:
    glClear(GL_COLOR_BUFFER_BIT);		    // clear the color buffer
    glClear(GL_DEPTH_BUFFER_BIT);           // clear the depth buffer

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();				        // make sure transformation is "zero'd"
    gluPerspective(60.0f,(GLfloat)viewport.w/(GLfloat)viewport.h,0.1f,100.0f);
    
    
    glMatrixMode(GL_MODELVIEW);			    // indicate we are specifying camera transformations
    glLoadIdentity();				        // make sure transformation is "zero'd"
    
    glTranslatef(transX, transY, transZ);
    glRotatef(rotY,0,1,0);  //rotate about y axis
    glRotatef(rotX,1,0,0);  //} rotate about x
    
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
        if (gravityOn)
            for (int j = 0; j < listOfSpheres.size(); j++){
                if (j == k){ continue; }
                sphere& s2 = listOfSpheres[j];
                if ((s2.pos-s1.pos).getNorm() < 0.001) continue;
                else s1.vel = s1.vel + 0.00000005*(s2.pos-s1.pos)*(s1.m+s2.m)*(1/(s1.m*(s2.pos-s1.pos).getNorm()));
            }
        
        //intersection loop
        for (int j = 0; j < listOfSpheres.size(); j++) {
            if (j == k) { continue; }
            sphere& s2 = listOfSpheres[j];
            if (s1.intersect(s2)) { collide(s1,s2); }
            // collide(s1,s2);
        }
        
        for (int j = 0; j < listOfPlanes.size(); j++) {
            plane p = listOfPlanes[j];
            if (s1.intersect(p)) { collide(s1,p); }
        }

	//applyVectorField(s1);	 
        s1.render();
    }
    
    if(loadFromFile) {
        if (fDataCounter < fData.size()) {
            vector<Vect3> thisFrame = fData[fDataCounter];
            for (int i = 0 ; i < thisFrame.size(); i++) {
                Vect3 pos = thisFrame[i];
                glTranslatef(pos.x,pos.y,pos.z);
                glutSolidSphere(0.2,sphAcc,sphAcc);
                glTranslatef(-pos.x,-pos.y,-pos.z);
            }
        } else {fDataCounter = 0;}
    }
    
    glFlush();
    glutSwapBuffers();					// swap buffers (we earlier set double buffer)
  
}
//}

int main(int argc, char *argv[]) {
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
    if (loadFromFile) {
        myParse(fname);  //}
    }
    
//{ Initialization of glut and window:  
    viewport.w = 600;
    viewport.h = 600;
    glutInit(&argc, argv);                        // This initializes glut  
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );  // Use a double-buffered window with red, green, and blue channels 
    glutInitWindowSize(viewport.w, viewport.h);   //The size and position of the window
    glutInitWindowPosition(0,0);      // x-,y- coords of the topLeft of new window.
    glutCreateWindow(argv[0]);        //} name of window.


    
//{ initScene() and callBack function:    
    initScene();
    glutDisplayFunc(myDisplay);	    // function to run when its time to draw something
    glutIdleFunc(myDisplay);	    // function to run when its time to draw something
    glutReshapeFunc(myReshape);	    // function to run when the window gets resized
    glutMouseFunc(myMsHndlr);
    glutKeyboardFunc(myKybdHndlr);
    glutSpecialFunc(myKybdHndlr);
    glutMainLoop();				    // infinite loop that will keep drawing and resizing
    return 0;                       //} never reaches here? 

}
