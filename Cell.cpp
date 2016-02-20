#include "Cell.h"

static int cid_counter = 0;

Cell::Cell() : R(RADIUS), K(NSTATES), cum_payoff(0) {
  cid = cid_counter++;
  Z = new bool[2*R];
  
  for (int r=0; r<2*R; r++)
    Z[r] = 0;

  Sx_q = new bool*[R];
  Cx_q = new Genome*[R];
  
  for (int r=0; r<R; r++) {
    Sx_q[r] = NULL;
    Cx_q[r] = NULL;
  }
  n_Sx_q = n_Cx_q = 0;

  nbh = NULL;
  cstructEnv(nbh);

  S[CRR]=S[NXT]=0;
  M[CRR]=M[NXT]=1;
  F[CRR]=F[NXT]=0;

  cable = new Wire[R];
}

Cell::Cell(const Cell& C) {
  nbh=NULL;
  cid = C.cid;
  R = C.R;
  K = C.K;
  n_Sx_q = C.n_Sx_q;
  n_Cx_q = C.n_Cx_q;
  S[CRR] = C.S[CRR];
  S[NXT] = C.S[NXT];
  M[CRR] = C.M[CRR];
  M[NXT] = C.M[NXT];
  G[CRR] = C.G[CRR];
  G[NXT] = C.G[NXT];
  F[CRR] = C.F[CRR];
  F[NXT] = C.F[NXT];
  cum_payoff = C.cum_payoff;

  Z = new bool[2*R];
  for (int i=0; i<2*R; i++)
    Z[i]=C.Z[i];

  cstructEnv(nbh);
  Sx_q = new bool*[R];
  Cx_q = new Genome*[R];

  for (int r=0; r<R; r++) {
    nbh[r]=C.nbh[r];
    Sx_q[r]=C.Sx_q[r];
    Cx_q[r]=C.Cx_q[r];
  }

  cable = new Wire[R];
  for (int r=0; r<R; r++)
    cable[r] = C.cable[r];
}

Cell::~Cell() {
  delete[] Z;
  delete[] Sx_q;
  delete[] Cx_q;
  delete[] cable;
  Z = NULL;
  Sx_q = NULL;
  Cx_q = NULL;
  n_Sx_q = n_Cx_q = 0;
}

Cell& Cell::operator = (const Cell& C) {
  cid = C.cid;
  R = C.R;
  K = C.K;
  n_Sx_q = C.n_Sx_q;
  n_Cx_q = C.n_Cx_q;
  S[CRR] = C.S[CRR];
  S[NXT] = C.S[NXT];
  M[CRR] = C.M[CRR];
  M[NXT] = C.M[NXT];
  G[CRR] = C.G[CRR];
  G[NXT] = C.G[NXT];
  F[CRR] = C.F[CRR];
  F[NXT] = C.F[NXT];
  cum_payoff = C.cum_payoff;

  for (int i=0; i<2*R; i++)
    Z[i]=C.Z[i];

  for (int r=0; r<R; r++) {
    nbh[r]=C.nbh[r];
    Sx_q[r]=C.Sx_q[r];
    Cx_q[r]=C.Cx_q[r];
  }
  for (int r=0; r<R; r++)
    cable[r] = C.cable[r];
  return *this;
}

void Cell::cstructEnv(Env& tmp) {
  dstructEnv(tmp);
  tmp = new Cell*[R];
  for (int r=0; r<R; r++)
    tmp[r]=NULL;
}

void Cell::dstructEnv(Env& tmp) {
  if (tmp)
    delete[] tmp;
  tmp=NULL;
}

bool Cell::cpEnv(const Env& src, Env& tgt) {
  if (!src || !tgt) return false;
  for (int r=0; r<R; r++)
    tgt[r]=src[r];
  return true;
}

bool Cell::queueNeighbor(Cell* C) {
  for (int r=0; r<R; r++)
    if (!nbh[r]) {
      nbh[r] = C;
      return true;
    }
  return false;
}

Wire* Cell::findWire(Cell* nbr) {
  for (int r=0; r<R; r++)
    if (nbh[r] == nbr) 
      return &cable[r];
  return NULL;
}

void Cell::mapWires() {
  for (int r=0; r<R; r++) {
    Wire* tmp_wire = nbh[r]->findWire(this);
    if (!tmp_wire) {
      cout << "serious error" << endl;
      return;
    }
    else {
      cable[r].connect(*this);
      cable[r].connect(*tmp_wire);
    }
  }
}

void Cell::readCable() {
  for (int r=0; r<R; r++) {
    cable[r].setV_in();
  }
}

void Cell::damageWire(const int& r, const int& damage) {
  if (r<0 || r>=R) {
    cout << "Cant damage wire, bad wire index." << endl;
    return;
  }
  cable[r].stDamage(damage);
}

void Cell::damageWire(const int& noise) {
  int ngood = (R-1)-numBadInput();

  if (ngood > 0) {
    int good[ngood];
    cout << "ngood: " << ngood << endl;
    int i=0;
    for (int r=1; r<R; r++) 
      if (!cable[r].damage) 
	good[i++]=r;
    int rnd = good[rand()%ngood];
    cable[rnd].stDamage(noise);
    cout << "set: " << rnd << " to " << noise << endl;
  }
}

void Cell::stNoise(const int& noise) {
  for (int r=0; r<R; r++) {
    if (cable[r].damage)
      cable[r].damage = noise;
  }
}

int Cell::numBadOutput() const {
  int num = 0;
  for (int r=0; r<R; r++) 
    if (cable[r].oput->gtDamage())
      num++;
  return num;
}

int Cell::numBadInput() const {
  int num = 0;
  for (int r=0; r<R; r++) 
    if (cable[r].gtDamage())
      num++;
  return num;
}

int Cell::rndOprNbr() {
  if (gtM() == false)
    cout << "error" << endl;
  int nopr=0;
  int rid[R];

  for (int r=0; r<R; r++) 
    if (nbh[r]->gtM())
      rid[nopr++]=r;

  return rid[rand()%nopr];
}

bool Cell::queueSx(bool* x) {
  if (!(n_Sx_q<R)) return false;
  Sx_q[n_Sx_q++]=x;
  return true;
}

bool Cell::queueCx(Genome* x) {
  if (!(n_Cx_q<R)) return false;
  Cx_q[n_Cx_q++]=x;
  return true;
}

void Cell::clearSx_q() {
  //for (int i=0; i<n_Sx_q; i++)
  //   Sx_q[i]=NULL;
  n_Sx_q = 0;
}

void Cell::clearCx_q() {
  //for (int i=0; i<n_Cx_q; i++)
  //  Cx_q[i]=NULL;
  n_Cx_q = 0;
}

const int& Cell::gtCid() const {
  return cid;
}

const int& Cell::gtR() const {
  return R;
}

const int& Cell::gtK() const {
  return K;
}

const bool& Cell::gtS() const {
  return S[CRR];
}

const bool& Cell::gtM() const {
  return M[CRR];
}

const int& Cell::gtF() const {
  return F[CRR];
}

const int& Cell::gtCumPayoff() const {
  int cum_wires = 0;

  for(int r=1; r<R; r++) 
    cum_wires += cable[r].total_payoff;
  
  if (cum_wires != cum_payoff) {
    cout << "Inconsistent total payoffs: " << endl;
  }

  return cum_payoff;
}

const Genome& Cell::gtG() const {
  return G[CRR];
}

void Cell::stS(const bool& s) {
  S[NXT] = s;
}

void Cell::stM(const bool& m) {
  M[NXT] = m;
}

void Cell::stF(const int& f) {
  F[NXT] = f;
}

void Cell::stG(const Genome& g) {
  G[NXT] = g;
}

void Cell::stGenome(const int& rule) {
  Genome g;
  g.set(rule);
  G[CRR] = G[NXT] = g;
}

void Cell::stGenome(const bool& s) {
  Genome g;
  g.set(s);
  G[CRR] = G[NXT] = g;
}

void Cell::updS() {
  S[CRR]=S[NXT];
}

void Cell::updM() {
  M[CRR]=M[NXT];
}

void Cell::updF() {
  F[CRR]=F[NXT];
}

void Cell::updG() {
  G[CRR]=G[NXT];
}

void Cell::pickS() {
  if (n_Sx_q == 0) return;

  stS(*Sx_q[rand()%n_Sx_q]);
  updS();
  clearSx_q();
}

void Cell::pickG() {
  if (n_Cx_q == 0) return;
  if (gtM()) return;

  stG(*Cx_q[rand()%n_Cx_q]);
  stM(true);
  updG();
  updM();
  clearCx_q();
}

void Cell::print(const bool& print_g) const {
  cout << "------------------------------" << endl;
  cout << "cid: " << cid << endl
       << "S: " << S[CRR] << " " << S[NXT] << endl
       << "M: " << M[CRR] << " " << M[NXT] << endl
       << "F: " << F[CRR] << " " << F[NXT] << endl
       << "n_Sx_q: " << n_Sx_q << endl
       << "n_Cx_q: " << n_Cx_q << endl;

  cout << "nbh: ";
  for (int r=0; r<R;r++)
    if (nbh[r])
      cout << nbh[r]->gtCid() << " ";
    else 
      cout << "N ";

  cout << endl << "Z: ";
  for (int x=0; x<2*R; x++)
    cout << Z[x];

  bool* tmp = new bool[2*R];
  int dec=0;
  cout << endl;

  if (print_g) {
    for (int y=0; y<G[CRR].gtDmy(); y++) {
      G[CRR].lookupGene(tmp,y,dec);
      cout << dec << "\t";
      for (int x=0; x<2*R; x++)
	cout << tmp[x];
      cout << endl;
    }
  }
  for (int x=0; x<R; x++)
    cable[x].print();
  cout << "------------------------------" << endl;
}


bool Cell::call(const int& r, const int& noise) {
  return *cable[r].V_in;

  if ((rand()%100) < noise)
    return nbh[r]->gtS() ? false : true;
  return nbh[r]->gtS();
}

void Cell::readEnv(bool*& s) {
  for (int r=0; r<R; r++)
    s[r]=nbh[r]->gtS();
}

void Cell::execute(int& z_I, int& z_dec) {
  if (!gtM()) return;

  int pow2 = 1;
  z_I = 0;
  for (int r=R-1; r>=0; r--) {
    //const bool s = *cable[r].V_in;
    const bool s = nbh[r]->gtS();
    z_I = z_I | (pow2*s);
    pow2 = pow2 << 1;
  }

  gtG().lookupGene(Z,z_I,z_dec);

  for (int r=0; r<R; r++) {
    nbh[r]->queueSx(&Z[2*r]);
    if (Z[2*r+1])
      nbh[r]->queueCx(&G[CRR]);
  }
}

void Cell::executeDamage(int& z_I, int& z_dec) {
  if (!gtM()) return;

  int pow2 = 1;
  z_I = 0;
  for (int r=R-1; r>=0; r--) {
    const bool s = *cable[r].V_in;
    z_I = z_I | (pow2*s);
    pow2 = pow2 << 1;
  }

  gtG().lookupGene(Z,z_I,z_dec);

  for (int r=0; r<R; r++) {
    nbh[r]->queueSx(&Z[2*r]);
    if (Z[2*r+1])
      nbh[r]->queueCx(&G[CRR]);
  }
}

void Cell::executeDamage(int& z_I, int& z_dec, bool*& z) {
  if (!gtM()) return;

  int pow2 = 1;
  z_I = 0;
  for (int r=R-1; r>=0; r--) {
    const bool s = *cable[r].V_in;
    z_I = z_I | (pow2*s);
    pow2 = pow2 << 1;
  }

  gtG().lookupGene(Z,z_I,z_dec);

  for (int r=0; r<2*R; r++)
    z[r]=Z[r];

  for (int r=0; r<R; r++) {
    nbh[r]->queueSx(&Z[2*r]);
    if (Z[2*r+1])
      nbh[r]->queueCx(&G[CRR]);
  }
}

int Cell::evaluate(int matrix[2][2]) {
  if (!gtM()) return 0;
  int poff = 0;
  for (int r=1; r<R; r++) 
    if (nbh[r]->gtM()) {
      poff+=matrix[gtS()][nbh[r]->gtS()];
    }
  //  cout << "poff: " << poff << endl;
  stF(poff);
  updF();
  cum_payoff+=poff;
  return poff;
}

int Cell::evaluateDamage(int matrix[2][2]) {
  if (!gtM()) return 0;
  int poff = 0;
  for (int r=1; r<R; r++) 
    if (nbh[r]->gtM())
      poff+=cable[r].payState(matrix);
  stF(poff);
  updF();
  cum_payoff+=poff;
  return poff;
}

void Cell::die() {
  if (!gtM()) return;
  for (int r=1; r<R; r++)
    if (gtF() < nbh[r]->gtF()) {
      stM(false);
      stF(0);
      return;
    }
}

void Cell::breed(const int& pmut, const int& max) {
  if (!gtM()) return;
  int rnd_nbr = rndOprNbr();
  const double pcross = 1.0-double(gtF()+nbh[rnd_nbr]->gtF())/double(max);
  G[NXT].crossover(gtG(),nbh[rnd_nbr]->gtG(),pcross);
  G[NXT].mutate(pmut);
  stF(0);
}
