#include "CA.h"
#include "Vector.h"
#include "fstream"

CA::CA() : R(5), K(2), dmx(0), dmy(0), model(NULL), gene_bin(new int[1024]) {
  for (int i=0; i<1024; i++)
    gene_bin[i] = 0;
  cum_fitness = fitness = 0;

  const int N = dmx*dmy;
  S = new bool[N];
  M = new bool[N];
  Z = new bool*[N];

  for (int n=0; n<N; n++) 
    Z[n] = new bool[2*R];
}

CA::CA(const int& x, const int& y) : R(5), K(2), dmx(x), dmy(y), model(NULL), gene_bin(new int[1024]) {
  const int N = dmx*dmy;
  model = new Cell[N];
  for (int i=0; i<1024; i++)
    gene_bin[i] = 0;
  cum_fitness = fitness = 0;
  S = new bool[N];
  M = new bool[N];
  Z = new bool*[N];
  for (int n=0; n<N; n++) 
    Z[n] = new bool[2*R];
}

CA::CA(const CA& ca) : R(ca.R), K(ca.K), dmx(ca.dmx), dmy(ca.dmy), model(NULL), gene_bin(new int[1024]) {
  const int N=dmx*dmy;
  model = new Cell[N];

  for (int i=0; i<N; i++) 
    model[i] = ca.model[i];

  for (int i=0; i<1024; i++)
    gene_bin[i] = ca.gene_bin[i];
  cum_fitness = ca.cum_fitness;
  fitness = ca.fitness;
  S = new bool[N];
  M = new bool[N];
  Z = new bool*[N];
  for (int n=0; n<N; n++) { 
    Z[n] = new bool[2*R];
    for (int r=0; r<2*R; r++) {
      Z[n][r] = ca.Z[n][r];
    }
    S[n]=ca.S[n];
    M[n]=ca.S[n];
  }   
}

CA::~CA() {
  delete[] model;
  delete[] gene_bin;
  model = NULL;
  gene_bin = NULL;

  delete[] S;
  delete[] M;
  for (int i=0; i<dmx*dmy; i++)
    delete[] Z[i];
  delete[] Z;
  S=NULL;
  M=NULL;
  Z=NULL;
}

void CA::print() const {
  const int N=dmx*dmy;
  for (int i=0; i<N; i++)
    model[i].print();
  //cout << "model: ";printBool(model,HH);
}

void CA::printGeneBin() const {
  for (int i=0; i<1024; i++)
    cout << i << "\t" << gene_bin[i] << endl;

  ofstream os("genebin.dat");

  bool* g = new bool[2*R];

  for (int i=0; i<1024; i++) {
    os << i << "\t" << gene_bin[i] << "\t";  
    fromDec(i,2,g,2*R);
    for (int r=0; r<2*R;r++)
      os << g[r];
    os << endl;
  }
  

  os.close();
}

CA& CA::operator = (const CA& ca) {
  R=ca.R;
  K=ca.K;

  if (!((dmx==ca.dmx) && (dmy==ca.dmy))) {
    dmx=ca.dmx;
    dmy=ca.dmy;
    delete[] model;
    model = new Cell[dmx*dmy];
  }

  const int N=dmx*dmy;
  for (int i=0; i<N; i++)
    model[i] = ca.model[i];

  cum_fitness = ca.cum_fitness;
  fitness = ca.fitness;
  for (int n=0; n<N; n++) {
    for (int r=0; r<2*R; r++)
      Z[n][r] = ca.Z[n][r];
    S[n]=ca.S[n];
    M[n]=ca.S[n];
  }  
  return *this;
}

const int& CA::Dmx() const {
  return dmx;
}

const int& CA::Dmy() const {
  return dmy;
}

bool CA::connectCells() {
  static int map[5][2] = {{0,0},{0,-1},{1,0},{0,1},{-1,0}};

  int i=0;
  for (int x=0; x<dmx; x++) 
    for (int y=0; y<dmy; y++) 
      for (int r=0; r<R; r++) {
	int rx = x+map[r][0];
	int ry = y+map[r][1];
	if (rx<0) 
	  rx+=dmx;
	else if (rx>=dmx) 
	  rx-=dmx;
	if (ry<0) 
	  ry+=dmy;
	else if (ry>=dmy) 
	  ry-=dmy;
	if (!model[y*dmx+x].queueNeighbor(&model[ry*dmx+rx])) {
	  cout << "Error pushing cell [" << rx << " " << ry 
	       << "] to cell [" << x << " " << y << "]" << endl;
	  return false;
	}
      }
  return true;
}

void CA::mapWires() {
  const int N=dmx*dmy;
  for (int i=0; i<N; i++) 
    model[i].mapWires();
}

void CA::readWires() {
  const int N=dmx*dmy;
  for (int i=0; i<N; i++) 
    model[i].readCable();
}

void CA::breakCells(const int& bias, const int& num1) {
  const int N = dmx*dmy;
  Vector<int> open_cells(N);

  const int num = int(double(N)*(double(num1)/100.0));

  cout << "num: " << num << endl;

  for (int i=0; i<N; i++) 
    open_cells(i)=i;

  int rnd_cell,rnd_wire;
  int destruction;

  for (int i=0; i<num; i++) {
    int rnd_c_idx = rand()%open_cells.Dim();
    rnd_cell = open_cells(rnd_c_idx);
    destruction = (rand()%(R-1))+1;

    Vector<int> open_wires(R-1);
    for (int j=1; j<R; j++)
      open_wires(j-1)=j;
    
    for (int k=0; k<destruction; k++) {
      int rnd_w_idx = rand()%open_wires.Dim();
      rnd_wire = open_wires(rnd_w_idx);
      model[rnd_cell].damageWire(rnd_wire,bias);
      open_wires.removeElement(rnd_w_idx);
    }
    open_cells.removeElement(rnd_c_idx);
  }
}

void CA::breakWires(const int& bias, const int& pbreak) {
  const int N = dmx*dmy;
  const int num_wires = N*(R-1);
  const int nbreaks = int((double(pbreak)/100.0)*double(num_wires));

  Vector<int> open(num_wires);
  Vector<int> closed;

  int m=0;
  for (int i=0; i<num_wires+N; i++)
    if (i%R)
      open(m++)=i;

  for (int i=0; i<nbreaks; i++) {
    int rnd = rand()%open.Dim();
    int val = open(rnd);
    int cell = val/R;
    int wire = val%R;
    model[cell].damageWire(wire,bias);
    open.removeElement(rnd);
    closed.insertElement(i,val);
  }
  cout << "broken " << nbreaks << " out of " << num_wires << " wires with a bias-value of " << bias << endl;
}

CA CA::makeReplica(const CA& origional) {
  const int dmx = origional.Dmx();
  const int dmy = origional.Dmy();
  const int N=origional.Dmx()*origional.Dmy();

  CA replica(dmx,dmy);
  replica.connectCells();
  replica.mapWires();
  
  for (int j=0; j<dmy; j++) {
    for (int i=0; i<dmx; i++) {
      replica(i,j).stS(origional(i,j).gtS());
      replica(i,j).stM(origional(i,j).gtM());
      replica(i,j).stG(origional(i,j).gtG());
      replica(i,j).updS();
      replica(i,j).updM();
      replica(i,j).updG();
    }
  }
  return replica;
}

void CA::rndState(const int& def) {
  const int N=dmx*dmy;
  
  for (int i=0; i<N; i++) {
    const int rnd = rand()%100;
    if (rnd<def)
      model[i].stS(false);
    else
      model[i].stS(true);
    model[i].updS();
  }
}

void CA::rndMode(const int& vac) {
  const int N=dmx*dmy;
  
  for (int i=0; i<N; i++) {
    const int rnd = rand()%100;
    if (rnd<vac)
      model[i].stM(false);
    else
      model[i].stM(true);
    model[i].updM();
  }
}

void CA::rndGenome() {
  const int N=dmx*dmy;
  Genome g;  

  for (int i=0; i<N; i++) {
    g.randomnize();
    model[i].stG(g);
    model[i].updG();
  }
}

void CA::rndGenome(const int& pdef) {
  const int N=dmx*dmy;
  Genome g;  

  for (int i=0; i<N; i++) {
    if ((rand()%100) < pdef) 
      g.set(false);
    else 
      g.set(true);
    model[i].stG(g);
    model[i].updG();
  }
}

void CA::stGenome(const string& fname) {
  ifstream is(fname.c_str());

  if (!is) {
    cout << "Cannot find rule file." << endl;
    return;
  }

  int dec,n=0;

  while (is>>dec) n++;

  is.clear();
  is.seekg(0);

  int rules[n];
  n=0;
  while (is>>rules[n]) n++;

  is.close();

  const int N=dmx*dmy;

  for (int i=0; i<N; i++) 
    model[i].stGenome(rules[rand()%n]);
}

void CA::stNoise(const int& noise){
  const int N=dmx*dmy;
  
  for (int i=0; i<N; i++) 
    model[i].stNoise(noise);
}

const int& CA::gtCumFitness() const {
  return cum_fitness;
}

const int& CA::gtFitness() const {
  return fitness;
}

void CA::propModel(int matrix[2][2], const int& pmut) {
  const int N=dmx*dmy;
  
  int gi=0,gdec=0;

  for (int i=0; i<N; i++) {
    model[i].execute(gi,gdec);
    gene_bin[gdec]++;
  }

  for (int i=0; i<N; i++) {
    model[i].pickS();
    model[i].pickG();
    model[i].updM();
  }

  for (int i=0; i<N; i++) 
    model[i].evaluate(matrix);
  for (int i=0; i<N; i++)
    model[i].die();

  for (int i=0; i<N; i++) {
    model[i].updM();
    model[i].updF();
  }
  
  int max0 = matrix[0][0] > matrix[0][1] ? matrix[0][0] : matrix[0][1];
  int max1 = matrix[1][0] > matrix[1][1] ? matrix[1][0] : matrix[1][1];
  int max = max0 > max1 ? max0 : max1;
  max = (R-1)*max*2;
  for (int i=0; i<N; i++) 
    model[i].breed(pmut,max);
    
  for (int i=0; i<N; i++) {
    if (model[i].gtM())
      model[i].updG();
  }
  
  for (int i=0; i<N; i++) {
    model[i].clearSx_q();
    model[i].clearCx_q();
  }
    
}

void CA::readCable() {
  const int N=dmx*dmy;
  for (int i=0; i<N; i++) {
    model[i].readCable();
  }
}

void CA::execute() {
  const int N=dmx*dmy;
  int gi=0,gdec=0;
  for (int i=0; i<N; i++) {
    model[i].executeDamage(gi,gdec);
    gene_bin[gdec]++;
  }
}

int CA::evaluate(int matrix[2][2]) {
  const int N=dmx*dmy;
  int poff = 0;
  for (int i=0; i<N; i++) 
    poff += model[i].evaluateDamage(matrix);
  return poff;
}


void CA::die() {
  const int N=dmx*dmy;
  for (int i=0; i<N; i++)
    model[i].die();
}

void CA::breed(const int& pmut, int matrix[2][2]) {
  const int N=dmx*dmy;
  int max0 = matrix[0][0] > matrix[0][1] ? matrix[0][0] : matrix[0][1];
  int max1 = matrix[1][0] > matrix[1][1] ? matrix[1][0] : matrix[1][1];
  int max = max0 > max1 ? max0 : max1;
  max = (R-1)*max*2;

  //breed
  for (int i=0; i<N; i++) 
    model[i].breed(pmut,max);
}

int tot_poff = 0;

void CA::damageModel(int matrix[2][2], const int& pmut) {
  const int N=dmx*dmy;
  
  readCable();
  execute();

  for (int i=0; i<N; i++) {
    model[i].pickS();
    model[i].pickG();
    model[i].updM();
  }

  readCable();
  fitness = evaluate(matrix);
  cum_fitness += fitness;
  die();

  for (int i=0; i<N; i++) {
    model[i].updM();
    model[i].updF();
  }
  
  breed(pmut,matrix);
    
  for (int i=0; i<N; i++) {
    if (model[i].gtM())
      model[i].updG();
  }
  
  for (int i=0; i<N; i++) {
    model[i].clearSx_q();
    model[i].clearCx_q();
  }

  for (int i=0; i<N; i++) {
    model[i].gtCumPayoff();
  }
}

void CA::coherenceModel(int matrix[2][2], const int& pmut, int bin[2][2], double avg_fit[2], double chere[2]) {
  const int N=dmx*dmy;

  for (int i=0; i<N; i++) {
    S[i]=model[i].gtS();
    M[i]=model[i].gtM();
  }

  for (int i=0; i<N; i++) 
    model[i].readCable();

  int z_I, z_dec;

  for (int i=0; i<N; i++) 
    model[i].executeDamage(z_I,z_dec,Z[i]);

  for (int i=0; i<N; i++) {
    model[i].pickS();
    model[i].pickG();
    model[i].updM();
  }

  bool* nbhs = new bool[R];
  double opr = 0;
  for (int i=0; i<N; i++) {
    if (M[i]) {
      opr++;
      model[i].readEnv(nbhs);
      double ctmp=0;
      for (int r=0; r<R; r++) {
	if (Z[i][2*r] == nbhs[r]) 
	  ctmp++;
      }
      ctmp/=double(R);
      chere[S[i]]+=ctmp;
    }
  }

  if (opr>0) {
    chere[0]/=opr;
    chere[1]/=opr;
  }

  delete[] nbhs;
  nbhs = NULL;

  for (int i=0; i<N; i++) 
    M[i]=false;

  for (int i=0; i<N; i++) 
    model[i].readCable();

  int fit = 0;
  fitness = 0;
  for (int i=0; i<N; i++) {
    fit = model[i].evaluateDamage(matrix);
    fitness+=fit;
    avg_fit[model[i].gtS()]+=fit;
  }
  cum_fitness += fitness;
  avg_fit[0]/=double(N);
  avg_fit[1]/=double(N);

  for (int i=0; i<N; i++)
    model[i].die();

  for (int i=0; i<N; i++) {
    model[i].updM();
    model[i].updF();
  }

  int max0 = matrix[0][0] > matrix[0][1] ? matrix[0][0] : matrix[0][1];
  int max1 = matrix[1][0] > matrix[1][1] ? matrix[1][0] : matrix[1][1];
  int max = max0 > max1 ? max0 : max1;
  max = (R-1)*max*2;

  //breed
  for (int i=0; i<N; i++) 
    model[i].breed(pmut,max);
    
  for (int i=0; i<N; i++) {
    if (model[i].gtM())
      model[i].updG();
  }
  
  for (int i=0; i<N; i++) {
    model[i].clearSx_q();
    model[i].clearCx_q();
    bin[model[i].gtM()][model[i].gtS()]++;
  }
}

void CA::fitnessIOerror(double fitness[5][5]) {
  const int N = dmx*dmy;

  double bin[R][R];

  for (int x=0; x<R; x++)
    for (int y=0; y<R; y++)
      bin[x][y]=fitness[x][y]=0;

  int errin,errout;
  for (int n=0; n<N; n++) {
    errin = model[n].numBadInput();
    errout = model[n].numBadOutput();
    fitness[errin][errout] += model[n].gtCumPayoff();
    bin[errin][errout]++;
  }

  cout << "bin: " << endl;
  for (int y=0; y<R; y++) {
    for (int x=0; x<R; x++) 
      cout << bin[x][y] << " ";
    cout << endl;
  }

  cout << "fitness: " << endl;
  for (int y=0; y<R; y++) {
    for (int x=0; x<R; x++) 
      cout << fitness[x][y] << " ";
    cout << endl;
  }  

  for (int y=0; y<R; y++) 
    for (int x=0; x<R; x++) 
      fitness[x][y]/=bin[x][y];

  cout << "fitness: " << endl;
  for (int y=0; y<R; y++) {
    for (int x=0; x<R; x++) 
      cout << fitness[x][y] << " ";
    cout << endl;
  } 
}

const Cell& CA::operator() (const int& x, const int& y) const {
  return model[y*dmx+x];
}

Cell& CA::operator() (const int& x, const int& y) {
  return model[y*dmx+x];
}

void CA::stUniform(const int& rule) {
  const int N=dmx*dmy;
  for (int i=0; i<N; i++)  
    model[i].stGenome(rule);
}

void CA::accCells(int bin[2][2]) {
  const int N=dmx*dmy;
  for (int i=0; i<N; i++)  
    bin[model[i].gtM()][model[i].gtS()]++;
}

void CA::saveModel(const string& fname) {
  const int N=dmx*dmy;
  ofstream os(fname.c_str());

  os << dmx << " " << dmy << " ";

  const int gsize = 320;

  for (int n=0; n<N; n++) {
    os << model[n].gtS() << " " << model[n].gtM() << " ";
    for (int i=0; i<gsize; i++)
      os << model[n].gtG()(i) << " ";
  }
  os.close();
}

void CA::loadModel(const string& fname) {
  cum_fitness = fitness = 0;

  if (S) {
    delete[] S;
    S=NULL;
  }

  if (M) {
    delete[] M;
    M=NULL;
  }

  if (Z) {
    for (int i=0; i<dmx*dmy; i++) 
      delete[] Z[i];
    delete[] Z;
    Z=NULL;
  }

  ifstream is(fname.c_str());

  if (!is) {
    cout << "> Cannot find file." << endl;
    return;
  }

  if (model) {
    delete[] model;
    model = NULL;
  }

  bool val;

  is >> dmx >> dmy;

  const int N = dmx*dmy;
  const int gsize = 320;

  model = new Cell[N];

  connectCells();
  mapWires();

  Genome g;

  for (int n=0; n<N; n++) {
    is >> val;
    model[n].stS(val);
    is >> val;
    model[n].stM(val);

    int k=0;
    for (int j=0;j<g.gtDmy(); j++) {
      for (int i=0;i<g.gtDmx(); i++) {
	is >> val;
	g(k++) = val;
      }
    }
    model[n].stG(g);
    model[n].updS();
    model[n].updM();
    model[n].updG();
  }
  is.close();

  S = new bool[N];
  M = new bool[N];
  Z = new bool*[N];

  for (int i=0; i<N; i++) 
    Z[i]=new bool[2*R];

  for (int i=0; i<1024; i++)
    gene_bin[i] = 0;
}
