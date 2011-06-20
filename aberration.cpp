/**       
 *	aberration.c - Richard Royal  
 *	version as of October 9, 2008
 *	PHY498A - Fall, 2008
 */

#include <GL/gl.h>   		/* OpenGL 							*/
#include <GL/glut.h> 		/* GLUT support library: windowing system, user input support.	*/
#include <GL/glui.h>
#include <math.h>	

#include <stdlib.h>		/* Input, key_assignments					*/
#include <stdio.h>		/* Output							*/

#include <iostream>

// macros - mapped to ouput sequence by preprossesor before compilation
#define SQR(x) ((x)*(x))
#define DOTPROD(x,y) ( (x[0]*y[0])+(x[1]*y[1])+(x[2]*y[2]) )
#define MAG(a) (  sqrt( (SQR(a[0]))+(SQR(a[1]))+(SQR(a[2])) )  )
#define CROSSPROD(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1]; 

// printf shortcuts
#define pV(a) ( printf("(%Lf, %Lf, %Lf)\n", a[0],a[1],a[2]) )		// print vector 
#define pS(a) ( (printf("%Lf\n", a )) ) 				// print scalar
#define pQ(a) ( printf("(%Lf, %Lf, %Lf, %Lf)\n", a[0],a[1],a[2],a[3]) )	// print quaternion


/*---------------physics variables------------------------*/
static long double pos[3] = {5, 0.0, 50.0};	// position of camera
static long double forward[3] = {0.0,0.0,-1.0};	// line of sight vector, mag = 1
long double v[3] = {0,0,0};
static long double angle=0.0;			// current angle from x axis in xz plane

double mass = 75.0;				// mass of user, kg
long double ke = 0.0;				// kenetic energy of user, J
double ke_step = 500.0;/*500*/				// add this much ke (J) at a time //500
long double c = 10.0;				// speed of light, m/s
long double dt = 10.0/1000.0;			// ms to s, numerator = update_timer
long double beta;				// relativistic beta, gamma
long double alpha;
long double cur_gamma;

float spect[1300][3];				// array of visible rbg values from spectrum on 'roygbiv' type scale

long double qpqi[4];				// quaternion triple product placeholder
long double dsCrossV[3];			// |ds||v|sin(theta)*[n] placeholder for nomalized [n] vector
long double collision = 1.0f;			// if user is within this distance of a vertex, then we call that a collision
/*--------------------------------------------------------*/

/*--------------- setup variables-------------------------*/
int mainWindow, overheadWindow, displayWindow,settingsWindow;	// Three main windows
int width = 900;                                // initial width of main viewing window
int height = 500;                               // initial height of main viewing window
static long double ratio;			// used in window resize
static long double angle_step = 0.03;		// amount of turn per input left/right

float vertices_array[5000][4];			// array of vertices in world, must be large enough to contain all 3D points
						// contains coordinates and color: index,x,y,z,freq
int num_vertices;				// stores the counted number of inputed vertices
int num_colors;					// the number of possible RGB values from EM spectrum
long double color_slope = 638.0/175.0;		// slope of line converting RGB to EM
int vert_per_polygon = 4;			// vertices per polygon, 3 for triagles, 4 for quads



int i,j; 					// generic iteration loop pointers
static unsigned int update_timer = 10.0;	// call updatePosition every x ms, 0.01s

int lookDownHeight = 25; 			// height from which the overhead display is positioned
int paused = 0;					// boolean: 1-pause geometry and position, 0-still moving
float grey[3] = {200.0/255.0,200.0/255.0,200.0/255.0};
float brown[3] = {127.0/255.0,0,0};
/*--------------------------------------------------------*/


/*-------------------Live variables-----------------------*/

GLUI *settings;

//int bkgrdColor[3] = {0,0,0};
int bkgrdColor[3] = {1,1,1};
int switchColor = 0;



int dop_on = 1;
int doppler = 1;
int usr_ke_step = ke_step;
int usr_c = c;
using namespace std;
/*--------------------------------------------------------*/





/* Opens GLUT window and loads world vertices into array */
void initScene(){

	// Load vertices and frequency into array
	float point_x, point_y, point_z,freq;
	char c;
	FILE *file;					// file to be opened, containing vertices of polygons
	file = fopen("stonehenge.txt", "r");		// opens file
	num_vertices = 0;
	// requires file to have "'x y z freq'\n" format
	while ( fscanf(file,"%f%c%f%c%f%c%f\n", &point_x,&c,&point_y,&c,&point_z,&c,&freq) != EOF){
		vertices_array[num_vertices][0] = point_x;
		vertices_array[num_vertices][1] = point_y;
		vertices_array[num_vertices][2] = point_z;
		vertices_array[num_vertices][3] = freq;
		num_vertices++;
	}
	fclose(file);

	// Load Color Map from text file
	float r, g, b;
	FILE *color;					// file to be opened, containing spectrum in RGB
	color = fopen("spectrum.txt", "r");		// opens file
	num_colors = 0;
	// requires file to have "'r g b'\n" format
	while ( fscanf(color,"%f%c%f%c%f\n", &r,&c,&g,&c,&b) != EOF){

		spect[num_colors][0]= r;
		spect[num_colors][1]= g;
		spect[num_colors][2]= b;
		num_colors++;
	}
	fclose(color);
}

/* Prints bitmap text onto window*/
void renderBitmapString(float px,float py, void *font, char *string) {  
	char *c;
	glRasterPos2f(px,py);
	glColor3f(1.0,1.0,1.0);						// Text color
	for (c=string; *c != '\0'; c++) {
		glutBitmapCharacter(font, *c);
	}
}

void headsUpDisplay(){
	glutSetWindow(displayWindow);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		// Clear The Screen And The Depth Buffer
	glClearColor(bkgrdColor[0], bkgrdColor[1], bkgrdColor[2], 0);
	glMatrixMode(GL_MODELVIEW);

	double gamma_d = cur_gamma;
	char gammaBuffer[20];							// array containing string to be printed on screen
	sprintf (gammaBuffer, "gamma: %f", gamma_d);				// combines string
	renderBitmapString(-150.0,0.0,GLUT_BITMAP_HELVETICA_18,gammaBuffer);	// calls function that draws text

	double beta_d = beta;
	char betaBuffer[20];							// array containing string to be printed on screen
	sprintf (betaBuffer, "beta: %f", beta_d);				// combines string
	renderBitmapString(-150.0,-20.0,GLUT_BITMAP_HELVETICA_18,betaBuffer);	// calls function that draws text

	int c_int = c;
	char cBuffer[20];						// array containing string to be printed on screen
	sprintf (cBuffer, "c: %d.0 m/s", c_int);			// combines string
	renderBitmapString(0.0,0.0,GLUT_BITMAP_HELVETICA_18,cBuffer);	// calls function that draws text

	int ke_int = fabs(ke);
	char keBuffer[20];							// array containing string to be printed on screen
	sprintf (keBuffer, "ke: %d J", ke_int);					// combines string
	renderBitmapString(0.0,-20.0,GLUT_BITMAP_HELVETICA_18,keBuffer);	// calls function that draws text

	char pausedBuffer[20];							// array containing string to be printed on screen
	sprintf (pausedBuffer, "PAUSED");					// combines string
	if(paused == 1){
		renderBitmapString(85.0,-20.0,GLUT_BITMAP_TIMES_ROMAN_10,pausedBuffer);	// calls function that draws text
	}

	glutSwapBuffers();
	glutPostRedisplay();
}

void pausE(){
	if (paused == 0){paused = 1;}
	else {paused = 0;}
}

void overheadDisplay(){

	glutSetWindow(overheadWindow);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		// Clear The Screen And The Depth Buffer
	glClearColor(bkgrdColor[0], bkgrdColor[1], bkgrdColor[2], 0);
	glMatrixMode(GL_MODELVIEW);
    	glLoadIdentity(); 						// reset drawing brush, required here!

   	// looking down from height above cur_pos
	gluLookAt(pos[0], pos[1]+lookDownHeight, pos[2], pos[0],pos[1]-1, pos[2], forward[0],0,forward[2]);

   	/*Draws overhead view of all polygons in vertices_array, as if infinite c*/
	int k=0;
	int l = 0;
	while( k<= (num_vertices)){
		l=0;
        	long double poly_color = vertices_array[k][3];
        	if (poly_color !=  0){
            		int cur_Color = (poly_color - 400.0)*(color_slope);
            		int rgb_Value[3] = {spect[cur_Color][0],spect[cur_Color][1],spect[cur_Color][2]};
            		glColor3f(rgb_Value[0], rgb_Value[1], rgb_Value[2]);				
        	}
		
		glBegin(GL_POLYGON);
			for(l=0;l<vert_per_polygon;l++){
				glVertex3f(vertices_array[k+l][0],  vertices_array[k+l][1],  vertices_array[k+l][2]);
			}
		glEnd();
            
		k = k+vert_per_polygon;
	}
    

    	/* Draws a colored triangle on current position, pos[x,0,z],
    	always points in direction of motion  */
	double kappa = 3.0; // scalar related to triangle size
    	double b[3]={-forward[2],0,forward[0]};   
    
    	int pointer_color_f = 850;
	glColor3f(spect[pointer_color_f][0],spect[pointer_color_f][1],spect[pointer_color_f][2]);

    	glBegin(GL_POLYGON);
        	glVertex3f(pos[0],0,pos[2]);
        	glVertex3f((pos[0]-kappa*forward[0])+b[0],0,(pos[2]-kappa*forward[2])+b[2]);
        	glVertex3f((pos[0]-kappa*forward[0])-b[0],0,(pos[2]-kappa*forward[2])-b[2]);
	glEnd();
	glFlush();	// ensures nothing gets "stuck" in frame buffer 

	glutSwapBuffers();

	glutPostRedisplay();
}

/* Deals with window resize from user */
void changeSize(int w, int h){

	if(h == 0)
		h = 1;
	ratio = 1.0f * w / h;
	// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the clipping volume
	gluPerspective(45,ratio,1,1000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	// gluLookAt(cur_position, forward vector (from cur_position), up vector)
	gluLookAt(pos[0], pos[1], pos[2], pos[0]+forward[0],pos[1]+forward[1], pos[2]+forward[2], 0.0f,1.0f,0.0f);

}


/* function returns speed of user given ke, mass, c */
long double speed(){ return c*sqrt(1-(SQR((mass*SQR(c))/(fabs(ke)+mass*SQR(c))))); }

long double getRotateAngle( long double v[3], long double ds[3], long double beta ){
	alpha = acos(  DOTPROD(v,ds) /  (  MAG(v) * MAG(ds) ) );			// angle between velocity and photon vector
	long double alpha_prime = acos( (cos(alpha)+beta)/(1+cos(alpha)*beta) );	// relativistically corrected alpha
	return fabs(alpha-alpha_prime);							// amount of rotation needed for sr correction
}

/*
 * 	Quaternion Triple Product for point p rotation ( p' = q * p * q_inverse )
 * 	p'[0,x',y',z'] = q * p[0,x,y,z] * qi
 * 	Stores product into qpqi[0,x',y',z']
 */
void quaternionTripProduct(long double ds[3], long double ax[3], long double angle){
	
	long double axisMag = MAG(ax);		
	long double n[3] = {ax[0]/axisMag, ax[1]/axisMag, ax[2]/axisMag}; 	// normalizes axis of rotation
	long double q[4] = {cos((angle/2.0)), sin((angle/2.0))*n[0],sin((angle/2.0))*n[1],sin((angle/2.0))*n[2]}; 	// quaternion q, used for rotations
	long double p[4] = {0,ds[0],ds[1],ds[2]};				// quaternion p [0,x,y,z]
	
	// quaternion math for q*p*qi = p' = [0,x',y',z']
	qpqi[0] = 0;
	qpqi[1] = q[2]*(p[3]*q[0] + p[2]*q[1] - p[1]*q[2]) - q[3]*(p[2]*q[0] - p[3]*q[1] + p[1]*q[3]) + q[0]*(p[1]*q[0] + p[3]*q[2] - p[2]*q[3]) - q[1]*(-p[1]*q[1] - p[2]*q[2] - p[3]*q[3]);

	qpqi[2] = -q[1] *(p[3]*q[0] + p[2]*q[1] - p[1]*q[2]) + q[0]*(p[2]*q[0] - p[3]*q[1] + p[1]*q[3]) + q[3]*(p[1]*q[0] + p[3]*q[2] - p[2]*q[3]) - q[2]*(-p[1]*q[1] - p[2]*q[2] - p[3]*q[3]);

	qpqi[3] =  q[0]*(p[3]*q[0] + p[2]*q[1] - p[1]*q[2]) + q[1]*(p[2]*q[0] - p[3]*q[1] + p[1]*q[3]) - q[2]*(p[1]*q[0] + p[3]*q[2] - p[2]*q[3]) - q[3]*(-p[1]*q[1] - p[2]*q[2] - p[3]*q[3]);

}

/* 
 * Main opengl drawing function
 * renderScene() is called as fast as opengl allows
 * contains all elements (quads, triangle, etc) that get displayed 
 */
void renderScene(void) {
	glutSetWindow(mainWindow);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and Depth Buffer
	glClearColor(bkgrdColor[0], bkgrdColor[1], bkgrdColor[2], 0);

	long double spd = speed();
	if (paused == 0){
		v[0] = forward[0]*spd;
 		v[2] = forward[2]*spd;
	}
	// program allows ke < 0, meaning forward and vel vector are anti-parallel
	if (ke <0){ 
		v[0] = -v[0];
		v[2] = -v[2];
	}

	beta = spd/c;
	cur_gamma = 1.0/sqrt(1-SQR(beta));
	
	i = 0;
	while( i<= (num_vertices)){
		long double polygon_frequency = vertices_array[i][3];

		if (beta!=0.0){

			glBegin(GL_POLYGON);
			j = 0;
			while(j<vert_per_polygon){
				// incoming photon vector, "infinite c"
				long double ds[3] = {vertices_array[i+j][0] - pos[0], vertices_array[i+j][1] - pos[1],vertices_array[i+j][2] - pos[2]};
				CROSSPROD(dsCrossV, ds,v); 							// store ( ds X v ) into dsCrossV[]
				long double rotate_angle = getRotateAngle(v,ds,beta);				// returns S.R. corrected angle
	
				quaternionTripProduct(ds,dsCrossV,rotate_angle); // stores quaternion triple product into qpqi, this is rotated ds

				// set the current color for the ith vertex
				if (dop_on==1){
					long double f_prime = polygon_frequency*(cur_gamma * (1+beta*cos(alpha)));
					if ( (f_prime > 400) && ( f_prime < 750 ) ){		// the polygon is visible
						int cur_color = ((f_prime - 400.0)*(color_slope));
						int rgb_value[3] = {spect[cur_color][0],spect[cur_color][1],spect[cur_color][2]};
						glColor3f(rgb_value[0], rgb_value[1], rgb_value[2]);				
					}
					else {
						if (f_prime < 400){		// too red, brown
							glColor3f(brown[0],brown[1],brown[2]);
						}
						if ( f_prime > 750){		// too blue, gray
							glColor3f(grey[0],grey[1],grey[2]);
						}
					}		
				}
				else{
					if (polygon_frequency !=  0){
						int cur_color = (polygon_frequency - 400.0)*(color_slope);
						int rgb_value[3] = {spect[cur_color][0],spect[cur_color][1],spect[cur_color][2]};
						glColor3f(rgb_value[0], rgb_value[1], rgb_value[2]);				
					}
				}		


				// vert_prime(qpqi) + cur_pos = vertex prime pos		
				glVertex3f(qpqi[1] + pos[0],qpqi[2] + pos[1],qpqi[3] + pos[2]);				
				j++;
			}	
			glEnd();					
		}
		else{		/* beta = 0, draw as if infinite c */
				// if beta = 0, then use inputed frequencies for colors directly
			if (polygon_frequency !=  0){
				int cur_color = (polygon_frequency - 400.0)*(color_slope);
				int rgb_value[3] = {spect[cur_color][0],spect[cur_color][1],spect[cur_color][2]};
				glColor3f(rgb_value[0], rgb_value[1], rgb_value[2]);				
			}

			glBegin(GL_POLYGON);
				for(j=0;j<vert_per_polygon;j++){
					glVertex3f(vertices_array[i+j][0],  vertices_array[i+j][1],  vertices_array[i+j][2]);
				}
			glEnd();	
		}
		glFlush();	// ensures nothing gets "stuck" in frame buffer 
		i = i+vert_per_polygon;
	}
	glutSwapBuffers();
	glutPostRedisplay();
}

/* Changes Forward Vector */
void changeView(long double ang) {
	// MAG(forward) = 1
	forward[0] = sin(ang);
	forward[2] = -cos(ang);
	// gluLookAt(cur_position, forward vector from cur_position - center of scene vector, up vector)
	gluLookAt(pos[0], pos[1], pos[2], pos[0]+forward[0],pos[1]+forward[1], pos[2]+forward[2], 0.0f,1.0f,0.0f);
}

/* Changes ke based on user input */
void changeKE(int sign){
	// input(sign): forward = 1, backward = -1
	if(sign==-1){
		if(ke > 0){ ke = 0.0; }
		else { ke = ke - ke_step; }
	}
	else if (ke<0){ ke = 0.0; }
	else{ ke = ke + ke_step; }	
}

/* Moves camera based on current ke */
/* function runs everytime 'update_timer' goes off */
void moveCamera(int value){

	// change in pos[x,z] =  gamma*v*dt 
	long double spd = speed();
	long double gamma = 1.0/sqrt(1-SQR(spd/c));

	if (paused == 0){	
	    // program allows ke < 0 
	    // if ke < 0, user is moving antiparallel to forward vector 
	    if (ke>0){
		    pos[0] = pos[0] + (forward[0]*spd*gamma*dt);
		    pos[2] = pos[2] + (forward[2]*spd*gamma*dt);
	    }
	    else{
		    pos[0] = pos[0] - (forward[0]*spd*gamma*dt);
		    pos[2] = pos[2] - (forward[2]*spd*gamma*dt);
	    }
	}		
	glLoadIdentity();
	
	// gluLookAt(cur_position, forward vector (from cur_position), up vector)
	gluLookAt(pos[0], pos[1], pos[2], pos[0]+forward[0],pos[1]+forward[1], pos[2]+forward[2], 0.0f,1.0f,0.0f);

	glutPostRedisplay();
	glutTimerFunc(update_timer,moveCamera,0);
}

void processNormalKeys(unsigned char key, int x, int y) {

	if (key == 27) 
		exit(0);				// esc key exits
	if (key == 32)
		pausE();                // spacebar pauses geometry
}

void idleScene(){
	headsUpDisplay();
	overheadDisplay();
	renderScene();

	//sleep(10);
}

void control_cb( int control ){

	if (control == 5){
		GLUI_Master.close_all();

		if (switchColor==0){

			bkgrdColor[0]=0;
			bkgrdColor[1]=0;
			bkgrdColor[2]=0;
		}
		else{
			bkgrdColor[0]=1;
			bkgrdColor[1]=1;
			bkgrdColor[2]=1;
		}

		if (doppler == 0){ dop_on = 0;}
		else{dop_on = 1;}
		
		c = usr_c;
		ke_step = usr_ke_step;

		//glutSetWindow(mainWindow);
		//glutShowWindow();
	}

//	printf("%d \n",control);
}

void createSettingsMenu(){

	//glutSetWindow(mainWindow);
	//glutHideWindow();
	GLUI *settings = GLUI_Master.create_glui("Settings",0);

	GLUI_Panel *color_panel = settings->add_panel("Background Color");
	GLUI_RadioGroup *group1 = settings->add_radiogroup_to_panel(color_panel,&switchColor,1,control_cb);
	settings->add_radiobutton_to_group(group1,"Black");
	settings->add_radiobutton_to_group(group1,"White");

	GLUI_Panel *doppler_panel = settings->add_panel("Doppler Effect");
	GLUI_RadioGroup *group2 = settings->add_radiogroup_to_panel(doppler_panel,&doppler,2,control_cb);
	settings->add_radiobutton_to_group(group2,"OFF");
	settings->add_radiobutton_to_group(group2,"ON");

	GLUI_Spinner *ke_spinner = settings->add_spinner("KE step (J):",GLUI_SPINNER_INT, &usr_ke_step);
	ke_spinner ->set_int_limits(500,50000);
	ke_spinner ->set_speed(.1);
	
	GLUI_Spinner *c_spinner = settings->add_spinner("C (m/s):",GLUI_SPINNER_INT, &usr_c);
	c_spinner ->set_int_limits(5,300000000);
	c_spinner ->set_speed(.1);


	new GLUI_Button( settings, "Update", 5,control_cb);
	//settings->set_main_gfx_window(mainWindow);

	/* We register the idle callback with GLUI, *not* with GLUT */
	//GLUI_Master.set_glutIdleFunc( myGlutIdle );
	GLUI_Master.set_glutIdleFunc( idleScene );	

}
/* User input function for movement */
void inputKey(int key, int x, int y) {
	
	switch (key) {
		case GLUT_KEY_LEFT  : angle -= angle_step;changeView(angle);break;
		case GLUT_KEY_RIGHT : angle += angle_step;changeView(angle);break;
		case GLUT_KEY_UP    : changeKE(1);break;
		case GLUT_KEY_DOWN  : changeKE(-1);break;
    case GLUT_KEY_F1    : createSettingsMenu();break;
	}
}

/* Main, no command line arguments accepted */
int main(int argc, char **argv){
	
	
	glutInit(&argc, argv);
	/* Select type of Display mode:   
	Double buffer, Color mode, Depth buffer, Alpha blending */  
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ALPHA);


	// Initialize Main Simulation Window
	glutInitWindowPosition(100,50);			// init position from upper lefthand corner on screen
	glutInitWindowSize(width,height);		// init size of window
	
//	glLoadIdentity(); 				// reset opengl
	
//			cout << "baller" << endl;
	mainWindow = glutCreateWindow("SR Simulation 3D");
	
	
	glutKeyboardFunc(processNormalKeys);		// Define fuctions called when glut recieves input
	glutSpecialFunc(inputKey);
	glutReshapeFunc(changeSize);

	
	glutDisplayFunc(renderScene);			// define main drawing funtion, runs recursively

	
	glutTimerFunc(update_timer,moveCamera,0);	// moves camera every (update_timer), ms	
	// Set Background Color 
	//glClearColor(bkgrdColor[0], bkgrdColor[1], bkgrdColor[2], 0);		
	glClearDepth(1.0);				// Enables Clearing Of The Depth Buffer

	glDepthFunc(GL_LESS);				// type of depth test to do.
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);			// Enables Smooth Color Shading, Default
	glMatrixMode(GL_MODELVIEW);			// 3-space

	initScene();					// Set up glut window and input vertices from file	

	displayWindow = glutCreateSubWindow(mainWindow,0,0,width/2,height/9);
	glClearColor(bkgrdColor[0], bkgrdColor[1], bkgrdColor[2], 0);
	glutDisplayFunc(headsUpDisplay);
	glutReshapeFunc(changeSize);


	overheadWindow = glutCreateSubWindow(mainWindow,width/2,0,width/2,height/9);
	glClearColor(bkgrdColor[0], bkgrdColor[1], bkgrdColor[2], 0);
	glutDisplayFunc(overheadDisplay);
	glutReshapeFunc(changeSize);


	glutMainLoop();						// Start OpenGL
	return(0);
}
