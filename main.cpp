#include "CA.h"
#include "itoa.h"

#include <sstream>
#include <pthread.h>
#include <GL/freeglut.h>

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <sys/types.h>
#ifdef HAVE_SYS_FILE_H
#  include <sys/file.h>
#endif
#include <sys/stat.h>

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#include <fcntl.h>
#include <stdio.h>
#include <errno.h>
 
#if defined (HAVE_STRING_H)
#  include <string.h>
#else 
#  include <strings.h>
#endif

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

#ifdef READLINE_LIBRARY
#  include "readline.h"
#  include "history.h"
#else
#  include <readline/readline.h>
#  include <readline/history.h>
#endif

const int DEFAULT_DMX = 50;
const int DEFAULT_DMY = 50;


#define FRAME_TITLE "Simulation Platform"

CA* automata = new CA(DEFAULT_DMX,DEFAULT_DMY);

int payoff_matrix[2][2] = {{1,5},{0,3}};//dd,dc,cd,cc
bool stop = true;
pthread_mutex_t command_mutex = PTHREAD_MUTEX_INITIALIZER;

//ofstream data_file("data.dat");

struct SimulationInfo {
  int dmx,dmy;
  int pvac, pdef;
  int pmut;
  int noise;
  bool cross;
  int epochs;//0 for infinite
};

SimulationInfo sim_info = { DEFAULT_DMY,DEFAULT_DMY,13,50,1,0,true,0 };

struct GLWindowInfo {
  int  wdmx,  wdmy,  wposx, wposy, xmin,  xmax, ymin,  ymax,  dmx, dmy, off;  
};

GLWindowInfo wndw =  {800,400,400,400,-100,100,-100,100,200,200,10};

string arr2str(int* arr, int N) {
  string arrstr = "";
  char* c = new char;
  for (int i=0; i<N; i++) {
    itoa(arr[i],c,10);
    arrstr+=c;
  }
  return arrstr;
}

string int2str(const int& x) {
  //  isstringstream iss
}

void drawLine (float x1, float y1, float x2, float y2) {
  glBegin (GL_LINES);
    glVertex2f (x1, y1);
    glVertex2f (x2, y2);
  glEnd ();
}

void drawBorder (void) {
  glLineWidth (2.0);
  drawLine ((float)(wndw.xmin+wndw.off), (float)(wndw.ymin+wndw.off),
            (float)(wndw.xmax-wndw.off), (float)(wndw.ymin+wndw.off));
  drawLine (wndw.xmax-wndw.off, wndw.ymin+wndw.off,wndw.xmax-wndw.off, wndw.ymax-wndw.off);
  drawLine (wndw.xmax-wndw.off, wndw.ymax-wndw.off,wndw.xmin+wndw.off, wndw.ymax-wndw.off);
  drawLine (wndw.xmin+wndw.off, wndw.ymax-wndw.off,wndw.xmin+wndw.off, wndw.ymin+wndw.off);
  glLineWidth (1.0);
}

void drawText (double x, double y, const char s[]) {
  int i, length;
  char c;

  glRasterPos2d (x, y);

  length = strlen (s);
  for (i=0; i<length; i++) { 
    c = s[i];
    glutBitmapCharacter (GLUT_BITMAP_8_BY_13, c);
  }
}

void drawGrid(int id=0) {
  const int dmx = sim_info.dmx;
  const int dmy = sim_info.dmy;
  const double gridoff = double(2.0*wndw.off);
  const double wminx = double(wndw.xmin)+gridoff;
  const double wminy = double(wndw.ymin)+gridoff;
  const double wmaxx = double(wndw.xmax)-gridoff;
  const double wmaxy = double(wndw.ymax)-gridoff;

  const double dx = (wmaxx-wminx)/double(dmx);
  const double dy = (wmaxy-wminy)/double(dmy);

  int M,A;
  bool damaged;
  
  for (int i=0; i<dmx; i++) {
    const double x = wminx+double(i)*dx;
    for (int j=0; j<dmy; j++) {
      const double y = wminy+double(j)*dy;
      const int nbreaks = (*automata)(i,j).numBadInput();
      double c[3] = {0.0,0.0,0.0};
      if (nbreaks > 0) {
	c[2]=1;
      }
      //else {
	M = (*automata)(i,j).gtM();
	A = (*automata)(i,j).gtS();
	
	if (A == 0) 
	  if (M == 0) 
	    c[0]=.2;
	  else
	    c[0]=1;
	else
	  if (M == 0)
	    c[1]=.2;
	  else
	    c[1]=1;
	//}
      glColor3f(c[0],c[1],c[2]);
      glRectf(x,y,x+dx,y+dy);
    }
  }
}

int timer = 1;
double avg_fitness[2] = {0.0,0.0};
int cell_bin[2][2] = {{0,0},{0,0}};
double coherence[2] = {0.0,0.0};

void drawRect(const int& x0, const int& y0, const int& w, const int& h) {
  glBegin(GL_LINES);
  glVertex2d(x0,y0);
  glVertex2d(x0+w,y0);
  glVertex2d(x0,y0);
  glVertex2d(x0,y0+h);
  glVertex2d(x0+w,y0);
  glVertex2d(x0+w,y0+h);
  glVertex2d(x0,y0+h);
  glVertex2d(x0+w,y0+h);
  glEnd();
}

void drawMonitor() {
  const double vpsize = 200.0;
  const double h = vpsize-(4*wndw.off);
  double hd = coherence[0]*h*.99;
  double hc = coherence[1]*h*.99;
  char* c = new char;

  double opr = double(cell_bin[0][1] + cell_bin[1][1]);
  double hd2 = (double(cell_bin[0][1])/opr)*h*.99;
  double hc2 = (double(cell_bin[1][1])/opr)*h*.99;

  glColor3f(0,0,1);
  drawRect(-88,-88,40,h+9);
  glColor3f(1,1,1);
  drawText(-86,80,"Coherence");

  glColor3f(1,0,0);
  glRectf(-86,-80,-70,-80+hd);
  //glRectf(-86,-80,-78,-80+hd);
  //glRectf(-78,-80,-70,-80+hd2);

  glColor3f(0,1,0);
  glRectf(-66,-80,-50,-80+hc);
  //glRectf(-66,-80,-58,-80+hc);
  //glRectf(-58,-80,-50,-80+hc2);


  glColor3f(1,1,1);
  itoa(int(coherence[0]*100.0),c,10);
  drawText(-82,-87,c);
  itoa(int(coherence[1]*100.0),c,10);
  drawText(-62,-87,c);

  //automata->accCells(cell_bin);

  
  const int textxoff = -44;
  const int textyoff = 70;
  const int texth = 8;
  const int textw = 105;

  const int N = sim_info.dmx*sim_info.dmy;

  glColor3f(0,0,1);
  drawRect(-46,-50,134,131);
  glColor3f(1,1,1);
  drawText(-44,80,"Monitor");

  drawText(textxoff,textyoff,"#cells:");
  itoa(sim_info.dmx*sim_info.dmy,c,10);
  drawText(textxoff+textw,textyoff,c);

  drawText(textxoff,textyoff-texth,"mutation-rate:");
  itoa(sim_info.pmut,c,10);
  drawText(textxoff+textw,textyoff-texth,c);

  drawText(textxoff,textyoff-2*texth,"damage:");
  itoa(sim_info.noise,c,10);
  drawText(textxoff+textw,textyoff-2*texth,c);

  drawText(textxoff,textyoff-3*texth,"matrix(dd,dc,cd,cc):");
  //itoa(1020304,c,10);
  drawText(textxoff+textw,textyoff-3*texth,"");

  drawText(textxoff,textyoff-4*texth,"epoch:");
  itoa(timer,c,10);
  drawText(textxoff+textw,textyoff-4*texth,c);

  drawText(textxoff,textyoff-5*texth,"defect-vacant:");
  itoa(int(100*(double(cell_bin[0][0])/double(N))),c,10);
  drawText(textxoff+textw,textyoff-5*texth,c);

  drawText(textxoff,textyoff-6*texth,"defect-operational:");
  itoa(int(100*(double(cell_bin[0][1])/double(N))),c,10);
  drawText(textxoff+textw,textyoff-6*texth,c);

  drawText(textxoff,textyoff-7*texth,"cooperate-vacant:");
  itoa(int(100*(double(cell_bin[1][0])/double(N))),c,10);
  drawText(textxoff+textw,textyoff-7*texth,c);

  drawText(textxoff,textyoff-8*texth,"cooperate-operational:");
  itoa(int(100*(double(cell_bin[1][1])/double(N))),c,10);
  drawText(textxoff+textw,textyoff-8*texth,c);

  drawText(textxoff,textyoff-9*texth,"defect-fitness:");
  itoa(int(avg_fitness[0]),c,10);
  drawText(textxoff+textw,textyoff-9*texth,c);

  drawText(textxoff,textyoff-10*texth,"cooperate-fitness:");
  itoa(int(avg_fitness[1]),c,10);
  drawText(textxoff+textw,textyoff-10*texth,c);

  glColor3f(0,0,1);
  drawRect(-46,-88,134,30);
  glColor3f(1,1,1);
  drawText(-44,-60,"Legend");

  drawText(0,-67,"Vacant");
  drawText(40,-67,"Operational");
  drawText(-42,-67-texth,"Defect");
  drawText(-42,-67-2*texth,"Cooperate");
  glColor3f(.5,0,0);
  glRectf(10,-70,15,-75);
  glColor3f(0,.5,0);
  glRectf(10,-78,15,-83);
  glColor3f(1,0,0);
  glRectf(60,-70,65,-75);
  glColor3f(0,1,0);
  glRectf(60,-78,65,-83);

  delete c;
  c=NULL;
}

void Display() {
  glClear (GL_COLOR_BUFFER_BIT);
  glColor3f (1.0, 0.0, 0.0);  

  pthread_mutex_lock(&command_mutex);

  cell_bin[0][0] = cell_bin[0][1] = cell_bin[1][0] = cell_bin[1][1] = 0;
  avg_fitness[0] = avg_fitness[1] = 0.0;
  coherence[0] = coherence[1] = 0.0;
  
  if (stop == false) {
    automata->coherenceModel(payoff_matrix,sim_info.pmut,cell_bin,avg_fitness,coherence);
    timer++;
  }

  //data_file << timer << " " << automata->gtFitness() << endl;
  
  
  glViewport(0,0,400,400);
  glColor3f(0,0,1);
  drawBorder();
  drawGrid();

  glViewport(400,0,400,400);
  glColor3f(0,0,1);
  drawBorder();
  drawMonitor();


  pthread_mutex_unlock(&command_mutex);

  glutSwapBuffers ();
  glutPostRedisplay();
}

int interval = 50;
void timerFunc(int val) {
  glutPostRedisplay();
  glutTimerFunc(interval, timerFunc, 1);
}

void Reshape (int w, int h) {
  glViewport (0, 0, (GLint) w, (GLint) h);

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  gluOrtho2D (wndw.xmin,
              wndw.xmin + (wndw.xmax-wndw.xmin) * ((GLdouble) w)/((GLdouble) wndw.wdmx),
              wndw.ymin,
              wndw.ymin + (wndw.ymax-wndw.ymin) * ((GLdouble) h)/((GLdouble) wndw.wdmy));

  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity ();
}

void Keyboard (unsigned char key, int x, int y) {
  switch (key) { 
  case 'p': stop = stop ? false : true;
    break;
    //case 't': timerFunc(1);
    //  cout 
    //break;
  }
}

void Mouse (int button, int state, int x, int y) {
  const double vpsize = 400;
  const double scale = vpsize/(wndw.xmax-wndw.xmin);
  const double grid_xmin = scale*2*wndw.off;
  const double grid_xmax = vpsize - grid_xmin;
  const double grid_ymin = scale*2*wndw.off;
  const double grid_ymax = vpsize - grid_ymin;

  if (!(x>=grid_xmin && x<=grid_xmax)) return;
  if (!(y>=grid_ymin && y<=grid_ymax)) return;

  y = vpsize-y;
  x-=grid_xmin;
  y-=grid_ymin;

  const double dx = (grid_xmax-grid_xmin)/double(sim_info.dmx);
  const double dy = (grid_ymax-grid_ymin)/double(sim_info.dmy);

  int cx = (int)(double(x)/dx);
  int cy = (int)(double(y)/dy);
  
 
  switch (button) { 
  case GLUT_LEFT_BUTTON:   
    if (state == GLUT_DOWN) { 
      (*automata)(cx,cy).print();
    }
      break;
  case GLUT_RIGHT_BUTTON:  
    if (state == GLUT_DOWN)
      (*automata)(cx,cy).damageWire(sim_info.noise);
      break;
  }
}

void glutInit(int argc, char** argv) {
  glutInit (&argc, argv);
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize (wndw.wdmx, wndw.wdmy); 
  glutInitWindowPosition (wndw.wposx,wndw.wposy);
  glutCreateWindow (FRAME_TITLE);
  glClearColor (0.0, 0.0, 0.0, 0.0);
  glShadeModel (GL_FLAT);
  glutDisplayFunc (Display); 
  glutReshapeFunc (Reshape);
  glutKeyboardFunc (Keyboard);
  glutMouseFunc (Mouse);
}



////////////////////////////////////////////////////////////////////////////////////////////

static char *progname = NULL;
static bool done = false;

int comHelp PARAMS((char *));
int comQuit PARAMS((char *));
int comStop PARAMS((char *));
int comStart PARAMS((char *));
int comRandomnizeStates PARAMS((char *));
int comRandomnizeModes PARAMS((char *));
int comRandomnizeGenome PARAMS((char *));
int comRandomnizeStrategies PARAMS((char *));
int comPayoffMatrix PARAMS((char *));
int comMutate PARAMS((char *));
int comBreakWires PARAMS((char *));
int comNoise PARAMS((char *));
int comDim PARAMS((char *));
int comMatrix PARAMS((char *));
int comTimerSpeed PARAMS((char *));

struct COMMAND {
  char *name;
  rl_icpfunc_t *func;
  char *doc;
};

COMMAND commands[] = {
  { (char*)"help", comHelp,
    (char*)"Help. Type 'name' to get help on a specific command."},
  { (char*)"quit", comQuit,
    (char*)"Exit the program"},
  { (char*)"stop", comStop,
    (char*)"Stop the simulation. Use 'start' to continue the simulation."},
  { (char*)"start", comStart,
    (char*)"Start the simulation. Use 'stop' to pause the simulation."},
  { (char*)"rnd-states", comRandomnizeStates,
    (char*)"Randomnize all cell-states.\n\tThe command takes the percentage of defective cells as argument, remaining cells will be set to cooperate."},
  { (char*)"rnd-mode", comRandomnizeModes,
    (char*)"Randomnize all cell-modes.\n\tThe command takes the percentage of vacant cells as argument, remaining cells will be set to operational."},
  { (char*)"rnd-genome", comRandomnizeGenome,
    (char*)"Randomnize the genome. Takes no arguments. Sets each bit in the genome to 1 or 0 with 5 percent chance."},
  { (char*)"rnd-strategies", comRandomnizeStrategies,
    (char*)"Randomnize all cell strategies. Randomnizes the genome of a cell according to certain strategies.\n\tThe command takes the number of defectors as argument, remaining cells will be set to be a cooperator.\n\tFor example, 'rnd-strategies 25' sets 25 percent of all cells to be a defector and 75 percent to be cooperator."},
  { (char*)"payoff-matrix", comPayoffMatrix,
    (char*)"Set the payoff matrix. Takes the rewards for combinations of two actions as arguments.\n\tIn total there are four possible scenarios, DD,DC,CD,DD, arguments are applied in respective order.\n\tFor example, 'payoff-matrix 1 5 0 3' defines the payoff-matrix for the prisoners dilemma. "},
  { (char*)"mutate", comMutate,
    (char*)"Set the mutation rate. Takes one value from 0 to 1000 as argument.\n\tDuring normal simulations a rate of 1 is gangbaar."},
  { (char*)"damage", comBreakWires,
    (char*)"Apply damage to the wires which cells use to communicate actions.\n\tTakes one arguments which determines the percentage percentage of broken wires.\n\tFor example, 'damage 50' breaks 50% of all wires."},
  { (char*)"noise", comNoise,
    (char*)"When a wire is damaged, the noise setting determines the percentage of information that is transferred incorrectly. For example 'noise 1' transfers 99% of all information correctly."},
  { (char*)"grid-size", comDim,
    (char*)"Specify x and y dimensions of the grid.\n\tFor example: 'grid-size 50 50' defines a grid of 50x50 cells." },
 { (char*)"timer-speed", comTimerSpeed,
    (char*)"Under construction."},
  { (char*)NULL, (rl_icpfunc_t *)NULL, (char*)NULL }
};

int comPayoffMatrix (char* arg) {
  istringstream iss(arg);
  int dd,dc,cd,cc;
  iss >> dd >> dc >> cd >> cc;
  pthread_mutex_lock(&command_mutex);
  payoff_matrix[0][0]=dd;
  payoff_matrix[0][1]=dc;
  payoff_matrix[1][0]=cd;
  payoff_matrix[1][1]=cc;
  pthread_mutex_unlock(&command_mutex);

  cout << "> Set the payoff matrix to dd: " << payoff_matrix[0][0] 
       << ", dc: " << payoff_matrix[0][1] 
       << ", cd: " << payoff_matrix[1][0]
       << ", cc: " << payoff_matrix[1][1] << "." << endl;
  return 0;
}

int comStop (char* arg) {
  pthread_mutex_lock(&command_mutex);
  stop = true;
  //cout << "> "
  pthread_mutex_unlock(&command_mutex);  
  return 0;
}

int comStart (char* arg) {
  stop = false;
  return 0;
}

 int comMutate (char* arg) {
  int rate = atoi(arg);

  if (rate >= 0 && rate <= 1000) {
    pthread_mutex_lock(&command_mutex);
    sim_info.pmut = rate;
    pthread_mutex_unlock(&command_mutex);  
    cout << "> Set mutation rate to " << sim_info.pmut << "." << endl;
  }
  else cout << "> Invalid value. Enter rate between 0 and 1000." << endl;
  return 0;
}

int comBreakWires(char* arg) {
  int perc = atoi(arg);

  pthread_mutex_lock(&command_mutex);
  automata->breakWires(sim_info.noise,perc);
  pthread_mutex_unlock(&command_mutex);
  return 0;
}

int comNoise(char* arg) {
  int noise = atoi(arg);

  pthread_mutex_lock(&command_mutex);
  sim_info.noise = noise;
  automata->stNoise(noise);
  pthread_mutex_unlock(&command_mutex);
  return 0;
}

int comDim(char* arg) {
  istringstream iss(arg);
  int dmx,dmy;
  iss >> dmx >> dmy;

  if (dmx <= 0 || dmy <= 0) {
    cout << "> Invalid dimensions. Values must be larger than 0." << endl;
    return 1;
  }
  pthread_mutex_lock(&command_mutex);
  delete automata;
  automata = new CA(dmx,dmy);
  automata->connectCells();
  automata->mapWires();
  sim_info.dmx = dmx;
  sim_info.dmy = dmy;
  pthread_mutex_unlock(&command_mutex);
  cout << "> Set grid dimensions to " << dmx << "x" << dmy << "." << endl;
  return 0;
}

int comRandomnizeStates(char* arg) {
  int pdef = atoi(arg);

  if (pdef >= 0 && pdef <= 100) {
    pthread_mutex_lock(&command_mutex);
    sim_info.pdef = pdef;
    automata->rndState(sim_info.pdef);
    pthread_mutex_unlock(&command_mutex);  
    cout << "> Randomnized " << sim_info.pdef << "% of cells to defect." << endl;
  }
  else cout << "> Invalid value. Enter rate between 0 and 100." << endl;
  return 0;
}

int comRandomnizeModes(char* arg) {
  int pvac = atoi(arg);

  if (pvac >= 0 && pvac <= 100) {
    pthread_mutex_lock(&command_mutex);
    sim_info.pvac = pvac;
    automata->rndMode(sim_info.pvac);
    pthread_mutex_unlock(&command_mutex);  
    cout << "> Randomnized " << sim_info.pvac << "% of cells to vacant." << endl;
  }
  else cout << "> Invalid value. Enter rate between 0 and 100." << endl;
  return 0;
}

int comRandomnizeGenome(char* arg) {
  pthread_mutex_lock(&command_mutex);
  automata->rndGenome();
  pthread_mutex_unlock(&command_mutex);
  cout << "> Randomnized cell-genome." << endl;
  return 0;
}

int comRandomnizeStrategies(char* arg) {
  int pdef = atoi(arg);
  pthread_mutex_lock(&command_mutex);
  automata->rndGenome(pdef);
  pthread_mutex_unlock(&command_mutex);
  cout << "> Randomnized cell-strategies." << endl;
  return 0;
}

int comHelp (char* arg) {
  register int i;
  int printed = 0;
  cout << "Command overview." << endl << endl;
  for (i = 0; commands[i].name; i++) {
      if (!*arg || (strcmp (arg, commands[i].name) == 0)) {
          printf ("%s\n\t%s.\n\n", commands[i].name, commands[i].doc);
          printed++;
      }
  }
  if (!printed) {
      printf ("No commands match `%s'.  Possibilties are:\n", arg);

      for (i = 0; commands[i].name; i++) {
	if (printed == 6) {
	  printed = 0;
	  printf ("\n");
	}
	
	printf ("%s\t", commands[i].name);
	printed++;
      }
      
      if (printed)
        printf ("\n");
  }
  return (0);
}

int comTimerSpeed(char* arg) {
  int speed = atoi(arg);

  pthread_mutex_lock(&command_mutex);
  stop = true;
  interval = speed;
  stop = false;
  pthread_mutex_unlock(&command_mutex);
  timerFunc(speed);
}

int comQuit (char* arg) {
  done = 1;
  return (0);
}


char* dupstr (char* s) {
  char *r;
  r = new char[strlen(s)+1];
  strcpy (r, s);
  return (r);
}

char* command_generator (const char* text, int state) {
  static int list_index, len;
  char *name;

  if (!state) {
    list_index = 0;
    len = strlen (text);
  }

  while (name = commands[list_index].name) {
    list_index++;
    
    if (strncmp (name, text, len) == 0)
      return (dupstr(name));
  }
  return ((char *)NULL);
}

char** fileman_completion (const char* text, int start, int end) {
  char **matches;

  matches = (char **)NULL;
  if (start == 0)
    matches = rl_completion_matches (text, command_generator);

  return (matches);
}


void initialize_readline () {
  rl_readline_name = "ca";
  rl_attempted_completion_function = fileman_completion;
}

char* stripwhite (char* string) {
  register char *s, *t;
  for (s = string; whitespace (*s); s++);
    
  if (*s == 0)
    return (s);

  t = s + strlen (s) - 1;
  while (t > s && whitespace (*t))
    t--;
  *++t = '\0';

  return s;
}

COMMAND* find_command (char* name) {
  register int i;

  for (i = 0; commands[i].name; i++)
    if (strcmp (name, commands[i].name) == 0)
      return (&commands[i]);

  return ((COMMAND *)NULL);
}


int execute_line (char* line) {
  register int i;
  COMMAND *command;
  char *word;

  i = 0;
  while (line[i] && whitespace (line[i])) i++;
  word = line + i;
  
  while (line[i] && !whitespace (line[i]))
    i++;
  
  if (line[i]) line[i++] = '\0';
  
  command = find_command (word);

  if (!command) {
    fprintf (stderr, "%s: No such command.\n", word);
    return (-1);
  }
  while (whitespace (line[i])) i++;
  word = line + i;
  return ((*(command->func)) (word));
}


void* commandThread(void* arg) {
  char *line, *s;
  initialize_readline();
  
  for ( ; done == 0; ) {
    line = readline ("# ");
    if (!line)
      break;
    
    s = stripwhite (line);
    if (*s) {
      add_history (s);
      execute_line (s);
    }
    delete (line);
  }
  return NULL;
}
/////////////////////////////////////////////////////////////////////////////////////////////////

void printBool(bool* arr, int N) {
  for (int i=0; i<N; i++)
    cout << arr[i];
  cout << endl;
}

inline int twoDec(bool* arr, const int& wrap) {
  unsigned int res = 0;
  for(int i=0; i<wrap; i++)
    res = res*2 + arr[i];
  return res;
}

void simulate() {
  const int nsim    = 20;
  const int simlen  = 900;
  const int breakpoint = 300;
  const int bailpoint =  1800;
  const double threshold = .95;

  automata->connectCells();
  automata->mapWires();

  cout.precision(2);

  const bool bad = false;
  const bool good = true;
  const bool dfct = false;
  const bool cprt = true;

  double data[2][2][simlen];

  for (int i=0; i<simlen; i++) {
    data[dfct][bad][i]=0;
    data[dfct][good][i]=0;
    data[cprt][bad][i]=0;
    data[cprt][good][i]=0;
  }

  int num_bad = 0;

  for (int s=0; (s-num_bad)<nsim; s++) {
    automata->rndState(sim_info.pdef);
    automata->rndMode(sim_info.pvac);
    automata->stGenome("rules.dat");
    
    cout << "Starting simulation " << (s-num_bad) << endl;

    cell_bin[0][0] = cell_bin[0][1] = cell_bin[1][0] = cell_bin[1][1] = 0;
    avg_fitness[0] = avg_fitness[1] = 0.0;
    coherence[0] = coherence[1] = 0.0;

    int t=0;
    while ((t++<bailpoint) && ((coherence[0]<threshold)&&(coherence[1]<threshold))) {
      automata->coherenceModel(payoff_matrix,0,cell_bin,avg_fitness,coherence);
      if (!(t%50)) cout << "epoch: " << (s-num_bad) << ":" << t
			  << "\tc[0]: " << fixed << coherence[0] 
			  << "\tc[1]: " << fixed <<coherence[1] << endl;
    }
    /*
    if (coherence[0]>threshold) {
      cout << "Designed a defective system with a coherence of " << coherence[0] 
	   << " under defective cells, and a coherence of " << coherence[1] 
	   << " under cooperative cells in " << t << " epochs." << endl;      
    }
    */
    if (coherence[0]>threshold) {
      cout << "Designed a cooperative system with a coherence of " << coherence[0] 
	   << " under defective cells, and a coherence of " << coherence[1] 
	   << " under cooperative cells in " << t << " epochs." << endl;
      automata->saveModel("base.mdl");

      cell_bin[0][0] = cell_bin[0][1] = cell_bin[1][0] = cell_bin[1][1] = 0;
      avg_fitness[0] = avg_fitness[1] = 0.0;
      coherence[0] = coherence[1] = 0.0;
      srand(0);

      int e;

      cout << endl << "Simulating system under normal conditions." << endl;
      for (e=0; e<simlen; e++) {
	automata->coherenceModel(payoff_matrix,sim_info.pmut,cell_bin,avg_fitness,coherence);
	data[dfct][good][e] += coherence[dfct];
	data[cprt][good][e] += coherence[cprt];
	if (!(e%50)) cout << "epoch: " << (s-num_bad) << ":" << e
			  << "\tc[0]: " << fixed << coherence[0] 
			  << "\tc[1]: " << fixed <<coherence[1] << endl;
      }

      CA replica;
      replica.loadModel("base.mdl");

      cell_bin[0][0] = cell_bin[0][1] = cell_bin[1][0] = cell_bin[1][1] = 0;
      avg_fitness[0] = avg_fitness[1] = 0.0;
      coherence[0] = coherence[1] = 0.0;
      srand(0);

      cout << endl << "Simulating system under faulty conditions." << endl;
      for (e=0; e<breakpoint; e++) {
	replica.coherenceModel(payoff_matrix,sim_info.pmut,cell_bin,avg_fitness,coherence);
	data[dfct][bad][e] += coherence[dfct];
	data[cprt][bad][e] += coherence[cprt];
	if (!(e%50)) cout << "epoch: " << (s-num_bad) << ":" << e
			  << "\tc[0]: " << fixed << coherence[0] 
			  << "\tc[1]: " << fixed <<coherence[1] << endl;
      }

      replica.breakWires(1,100);

      for (e; e<simlen; e++){
	replica.coherenceModel(payoff_matrix,sim_info.pmut,cell_bin,avg_fitness,coherence);
	data[dfct][bad][e] += coherence[dfct];
	data[cprt][bad][e] += coherence[cprt]; 
	if (!(e%50)) cout << "epoch: " << (s-num_bad) << ":" << e
			  << "\tc[0]: " << fixed << coherence[0] 
			  << "\tc[1]: " << fixed <<coherence[1] << endl;
      }
      
    }
    else {
      cout << "Bad system. Next round, new chances." << endl;
      num_bad++;
    }
  }

  ofstream osc("osc2.dat");
  ofstream osd("osd2.dat");

  for (int e=0; e<simlen; e++) {
    data[dfct][bad][e]/=double(nsim);
    data[cprt][bad][e]/=double(nsim);
    data[dfct][good][e]/=double(nsim);
    data[cprt][good][e]/=double(nsim);

    osc << e << " " << data[cprt][bad][e] << " " << data[cprt][good][e] << endl;
    osd << e << " " << data[dfct][bad][e] << " " << data[dfct][good][e] << endl;
  }
  osc.close();
  osd.close();
}


int main1(int argc, char** argv) {
  cout << endl << "main1" << endl;
  return 0;
}

int main2(int argc, char** argv) {
  int seed;
  if (argc == 3)
    seed = atoi(argv[2]);
  else seed = time(NULL);
  srand(seed);
  cout << endl << endl << "Seed: " << seed << endl;
  cout << endl << "> type help for more information" << endl << endl;
  pthread_t command_thread;
  

  automata->connectCells();
  automata->mapWires();
  //automata->rndState(sim_info.pdef);
  //automata->stUniform(853);
  automata->rndGenome();
  //automata->stGenome("rules.dat");
  //automata->rndMode(sim_info.pvac); 
  //automata->breakWires(50,500);

  glutInit(argc,argv);
  pthread_mutex_init(&command_mutex,NULL);
  pthread_create(&command_thread,NULL,commandThread,NULL);

  for (; done==0; )
    glutMainLoopEvent();

  cout << endl << endl;
  pthread_join(command_thread,NULL);
  pthread_exit(NULL);
  //data_file.close();
  return 0;
}

int main(int argc, char** argv) {
  return main2(argc,argv);
  /*
  if (argc==1)
    return main1(argc,argv);
  if (atoi(argv[1])==1)
    return main1(argc,argv);
  if (atoi(argv[1])==2)
    return main2(argc,argv);
  */
  return 0;
}


/*

class Wire {
protected:
  int wid;
  bool signal;
  Action* data;
public:
  Wire();
  ~Wire();
  
  void writeData(Action* a);
  void readData(Action* a);
  const bool& gtSignal() const;
  void stSignal(const bool& s);
  void print() const;
};


Wire::Wire() {
  wid=wid_counter++;
  signal = 0;
  data = new Action;
}

Wire::~Wire() {
  delete data;
}

void Wire::writeData(Action* a) {
  *data = *a;
  stSignal(true);
}

void Wire::readData(Action* a) {
  *a = *data;
  stSignal(false);
}

const bool& Wire::gtSignal() const {
  return signal;
}

void Wire::stSignal(const bool& s) {
  signal = s;
}

void Wire::print() const {
  cout << "wid: " << wid << endl;
}

int cid_counter = 0;

class Cell {
 private:
  int cid;
  Action A;
  Wire** iput;
  Wire** oput;
 public:
  Cell();
  ~Cell();

  void connectCell(Cell& C);
  Wire* freeWire(Wire** wires);

  void writeWires();
  void readWires();
  void print() const;

  const static int R = 3;
  const static int K = 2;
};

Cell::Cell() {
  cid=cid_counter++;
  iput = new Wire*[R];
  oput = new Wire*[R];

  for (int r=0; r<R; r++) {
    oput[r] = NULL;
    iput[r] = NULL;
  }
  A = Action(rand()%2);
}

Cell::~Cell() {}

Wire* Cell::freeWire(Wire** wires) {
  for (int r=0; r<R; r++)
    if (!wires[r]) {
      wires[r] = new Wire;
      return wires[r];
    }
  return NULL;
}

void Cell::connectCell(Cell& C) {
 for (int r=0; r<R; r++)
   if (!iput[r]) {
     iput[r] = freeWire(C.oput);
     return;
   }
}

void Cell::writeWires() {
  for (int r=0; r<R; r++)
    oput[r]->writeData(&A);
}

void Cell::readWires() {
  Action* tmp = new Action;
  for (int r=0; r<R; r++) {
    iput[r]->readData(tmp);
  }
}

void Cell::print() const {
  cout << "cid: " << cid << endl;

  cout << "wires out: " << endl;
  for (int r=0; r<R; r++) {
    if(oput[r]) oput[r]->print();
    else cout << "NULL" << endl;
  }

  cout << "wires in: " << endl;
  for (int r=0; r<R; r++) {
    if(iput[r]) iput[r]->print();
    else cout << "NULL" << endl;
  }    
}

*/


 /*
  Cell cl,cc,cr;

  if (!cl.pushEnvIn(cl)) cout << "error" << endl;
  if (!cl.pushEnvIn(cr)) cout << "error" << endl;
  if (!cl.pushEnvIn(cc)) cout << "error" << endl;

  if (!cc.pushEnvIn(cc)) cout << "error" << endl;
  if (!cc.pushEnvIn(cl)) cout << "error" << endl;
  if (!cc.pushEnvIn(cr)) cout << "error" << endl;

  if (!cr.pushEnvIn(cr)) cout << "error" << endl;
  if (!cr.pushEnvIn(cc)) cout << "error" << endl;
  if (!cr.pushEnvIn(cl)) cout << "error" << endl;
  */
