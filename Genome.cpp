#include "Genome.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <cassert>

using namespace std;

inline unsigned int toDec(bool* arr, const int& wrap) {
  unsigned int res = 0;
  for(int i=0; i<wrap; i++)
    res = res*2 + arr[i];
  return res;
}

inline void fromDec(int dec, int base, bool*& arr, int wrap) {
  for (int i=0; i<wrap; i++)
    arr[i]=false;
  while (dec > 0) {
    wrap-=1;
    arr[wrap] = dec%base;       
    dec = dec/base;
  }     
}

Genome::Genome() : dmx(10), dmy(32), Z(NULL) {
 const int N=dmx*dmy;
  Z = new bool[dmx*dmy];
  for (int i=0; i<N; i++)
    Z[i]=false;
}

Genome::Genome(const Genome& G) : dmx(10), dmy(32), Z(NULL) {
  const int N=dmx*dmy;
  Z = new bool[N];
  for (int i=0; i<N; i++)
    Z[i]=G.Z[i];
}

Genome::Genome(const Genome& G1, const Genome& G2, double pcross) : dmx(10), dmy(32), Z(NULL) {
  const int N=dmx*dmy;
  Z = new bool[N];

  if (pcross<0.0) {
    cout << "Warning set pcross to 0" << endl;
    pcross=0;
  }
  if (pcross>1.0) {
    cout << "Warning set pcross to 1" << endl;
    pcross=1;
  }

  int p = (int)(pcross*(double)N);

  int i;
  for (i=0; i<p; i++)
    Z[i]=G1.Z[i];

  for (i; i<N; i++)
    Z[i]=G2.Z[i];
}

Genome::~Genome() {
  delete[] Z;
  Z=NULL;
}

void Genome::lookupGene(bool*& g, const int& y, int& g_dec) const {
  int k = (y+1)*dmx;

  int pow2 = 1;
  g_dec = 0;
  for (int x=dmx-1; x>=0; x--) {
    g[x]=Z[--k];
    g_dec = g_dec | (pow2*g[x]);
    pow2 = pow2 << 1;
  }

  //for (int i=0; i<10; i++)
  //  cout << g[i];
  //cout << endl;

}

void Genome::mutate(const int& rate) {
  const int N=dmx*dmy;
  int rnd;

  double n = double(N)*(double(rate)/1000.0);
  double fn = floor(n);
  double rest = n-fn;

  int in = (int)n;
  int irest = (int)(rest*10);

  for (int i=0; i<in; i++) {
    rnd = rand()%N;
    Z[rnd] = Z[rnd] ? false : true;
  }

  rnd = rand()%10;

  if (rnd < irest) {
    rnd = rand()%N;
    Z[rnd] = Z[rnd] ? false : true;
  }
    
  return;
  
  
  for (int i=0; i<N; i++) {
    rnd = rand()%1000;
    if (rnd<rate)
      Z[i] = Z[i] ? false : true;
  }
}

void Genome::randomnize() {
  const int N=dmx*dmy;
  for (int i=0; i<N; i++)
    Z[i]=rand()%2;
}

void Genome::set(const int& dec) {
  bool* bin = new bool[dmx];

  for (int r=0; r<dmx; r++) 
    bin[r]=false;

  fromDec(dec,2,bin,dmx);

  int i=0;
  for (int y=0; y<dmy; y++)
    for (int x=0; x<dmx; x++)
      Z[i++]=bin[x];
  delete[] bin;
  bin=NULL;
}

void Genome::set(const bool& S) {
  const int N = dmx*dmy;
  
  for (int n=0; n<N/2; n++) {
    Z[n]=S;
    Z[n+1]=rand()%2;
  }
}

const int& Genome::gtDmx() const {
  return Genome::dmx;
}

const int& Genome::gtDmy() const {
  return Genome::dmy;
}


Genome& Genome::operator = (const Genome& G) {
  const int N=dmx*dmy;
  for (int i=0; i<N; i++)
    Z[i]=G.Z[i];
  return *this;
}

void Genome::crossover(const Genome& G1, const Genome& G2, double pcross) {
  const int N=dmx*dmy;

  if (pcross<0.0) {
    cout << "Warning set pcross to 0" << endl;
    pcross=0;
  }
  if (pcross>1.0) {
    cout << "Warning set pcross to 1" << endl;
    pcross=1;
  }

  int p = (int)(pcross*(double)N);

  int i;
  for (i=0; i<p; i++)
    Z[i]=G1.Z[i];

  for (i; i<N; i++)
    Z[i]=G2.Z[i];
}

const bool& Genome::operator () (const int& i) const {
  assert(i<dmx*dmy);
  return Z[i];
}

bool& Genome::operator () (const int& i) {
  assert(i<dmx*dmy);
  return Z[i];
}


ostream& operator << (ostream& os, const Genome& G) {
  const int N=G.dmx*G.dmy;
  int i=0;
  for (int y=0; y<G.dmy; y++){
    for (int x=0; x<G.dmx; x++)
      os << G(i++) << " ";
    os << endl;
  }
  return os;
}
