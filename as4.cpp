//{ Includes, Declarations of functions/classes:
#include "ColorAndVector.h"

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
#define sphAcc 20 //higher number => prettier spheres.
#define thresholdForBounce 0 //higher number => bounces happen sooner.

inline float sqr(float x) { return x*x; }
using namespace std;
void myDisplay();

class Viewport;
class pt;

bool debug = false;
int firstTime = 1;
int prevCounter = 0;
int counter = 0;

void appendToFile(string fnameParam, string toAppend){
    ofstream outfile;
    const char * fname = fnameParam.c_str();
    outfile.open(fname, ios_base::app);
    outfile << toAppend; 
    outfile.close();
}

//}

//{ Classes (Viewport, pt):
class Viewport {
  public:
    int w, h; // width and height
};

class sphere {
    public:
     Vect3 pos;
     Vect3 vel;
     float r;
     float m;
     sphere();
     sphere(Vect3, Vect3, float);
     sphere(Vect3, Vect3, float, float);
     void render();
     bool intersect(sphere);
     void move();
    };
sphere::sphere(){
pos = Vect3();
vel = Vect3();
r = 0;
m = 1;
}    
sphere::sphere(Vect3 posp, Vect3 velp, float rp){
pos = posp;
vel = velp;
r=rp;
m = 1;
}
sphere::sphere(Vect3 posp, Vect3 velp, float rp, float massp){
pos = posp;
vel = velp;
r=rp;
m = massp;
}
void sphere::render(){
    glTranslatef(pos.x,pos.y,pos.z);
    
    // glutWireSphere(r,sphAcc,sphAcc);
    glutSolidSphere(r,sphAcc,sphAcc);
    
    glTranslatef(-pos.x,-pos.y,-pos.z);
    
    if (counter != prevCounter){
        appendToFile("test1","EOF\n");
        prevCounter = counter;
    }
    
    string str = "";
    ostringstream ss;
    ss << pos.x << " " << pos.y << " " << pos.z << "\n";
    appendToFile("test1",ss.str());
    
}
bool sphere::intersect(sphere s2){
    float sumOfRadii = r + s2.r;
    float distBetweenRadii = (pos - s2.pos).getNorm();
    float surfaceDistance = distBetweenRadii - sumOfRadii;
    if (surfaceDistance < thresholdForBounce) { return true; }
    else { return false; }
}
void sphere::move(){   
pos = pos + 1*vel;
}
//}

//{ Global Variables
Viewport viewport;
GLfloat light_diffuse[] = {0.2, 0.1, 1, 1.0};  /* Red diffuse light. */
GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};  /* Infinite light location. */


vector<sphere> listOfSpheres;
int prevX, prevY;

time_t initTime;
vector<vector<Vect3> > fData;
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
      // Vect3 vel = 0.02 * normalize(Vect3(-xx,-yy,0)); //Make balls always start with vel TOWARDS CENTER.
      Vect3 vel = lenOfDrag * 0.02 * normalize(Vect3(x-prevX,-y+prevY,r)); //MOUSE DRAG decides direction of vel.
      // listOfSpheres.push_back(sphere(Vect3(xx,yy,0),Vect3(0.02*r,0.02*r,0),0.2)); //RANDOM vel dir.
      listOfSpheres.push_back(sphere(Vect3(xx,yy,0),vel,0.2));
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

} //}

void initScene(){
    glLineWidth(0.5);
    glColor3f(1,1,1);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    
    int numCubed = 0;
    for (int i = 0; i < numCubed; i++){
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
    // float t1 = s1.vel.getNorm();
    // float t2 = s2.vel.getNorm();
    // Vect3 dir = 0.04 * normalize(s1.pos-s2.pos);
    // s1.vel = t2 *  normalize(s1.vel + dir);
    // s2.vel = t1 *  normalize(s2.vel - dir);
    
    
    Vect3 vel1 = s1.vel;
    Vect3 vel2 = s2.vel;
    
    s1.vel = (((s2.m)*(vel2 - vel1)) + (s2.m * vel2) + (s1.m * vel1))*(1/(s1.m+s2.m));
    s2.vel = (((s1.m)*(vel1 - vel2)) + (s2.m * vel2) + (s1.m * vel1))*(1/(s1.m+s2.m));
    
    //s1.vel = ((s1.m-s2.m)*vel1 + 2 * s2.m * vel2)*(1/(s1.m+s2.m));
    //s2.vel = (2*s1.m*vel1 - (s1.m-s2.m)*vel2)*(1/(s1.m+s2.m));
    
    // s1.vel = (2 * s2.vel)*(1/(s1.m+s2.m));
    // s2.vel = (2 * s1.vel)*(1/(s1.m+s2.m));
    
    Vect3 pos1 = s1.pos;
    Vect3 pos2 = s2.pos;
    float diff = (s2.pos-s1.pos).getNorm() * (s1.r+s2.r);
    s1.pos = pos1+(pos1-pos2)*((diff)/2);
    s2.pos = pos2+(pos2-pos1)*((diff)/2);

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
    cout << (float) counter/(finalTime - initTime) << endl;
    }
    counter++;
    for (int k = 0; k < listOfSpheres.size(); k++){
        sphere& s1 = listOfSpheres[k];

        s1.move();

        
        // gravity loop
        for (int j = 0; j < listOfSpheres.size(); j++){
            if (j == k){ continue; }
            sphere& s2 = listOfSpheres[j];
            s1.vel = s1.vel + 0.00005*(s2.pos-s1.pos)*(s1.m+s2.m)*(1/((s2.pos-s1.pos).getNorm()));
        }
        
        //intersection loop
        for (int j = 0; j < listOfSpheres.size(); j++){
            if (j == k){ continue; }
            sphere& s2 = listOfSpheres[j];
            if (s1.intersect(s2)){ collide(s1,s2); }
        }
        
        s1.render();
    }
    
    /*
    if (counter<fData.size()){
        vector<Vect3> thisFrame = fData[counter];
        for (int i = 0 ; i < thisFrame.size(); i++){
            Vect3 pos = thisFrame[i];
            glTranslatef(pos.x,pos.y,pos.z);
            glutSolidSphere(0.2,sphAcc,sphAcc);
            glTranslatef(-pos.x,-pos.y,-pos.z);
        }
    } else {counter = 0;}
    */

    glFlush();
    glutSwapBuffers();					// swap buffers (we earlier set double buffer)
  
}
//}

int main(int argc, char *argv[]) {
//{ PARSING:
    // if (argc < 3) { cout << "Please provide the filename and subdivision_Parameter" << endl; exit(0); }
    //string fname = argv[1];
	string fname = "test1";
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
    myParse(fname);  //}
 
    //GL polygon mode (filled / line)
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
