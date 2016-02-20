#ifndef __CELL_H__
#define __CELL_H__

#include "Genome.h"
#include "Wire.h"

const static int RADIUS = 5;
const static int NSTATES = 2;

enum Temp {CRR,NXT};

class Cell;
class CA;
typedef Cell** Env;

class Wire;

class Cell {
private:
  int cid;
  int R,K;
  bool* Z;
  bool** Sx_q;
  Genome** Cx_q;
  int n_Sx_q, n_Cx_q;
  Env nbh;
  Wire* cable;
  bool S[2];
  bool M[2];
  int F[2];
  int cum_payoff;
  Genome G[2];

  void cstructEnv(Env& tmp);
  void dstructEnv(Env& tmp);
  bool cpEnv(const Env& src, Env& tgt);
  int rndOprNbr();

  friend class Wire;
  friend class CA;
public:
  Cell();
  Cell(const Cell& C);
  ~Cell();

  bool queueNeighbor(Cell* C);
  Wire* findWire(Cell* nbr);
  void mapWires();
  void damageWire(const int& r, const int& damage);
  void damageWire(const int& noise);
  void stNoise(const int& noise);
  void damageCable(const int& n, const int& damage);
  int numBadOutput() const;
  int numBadInput() const;
  void readCable();
  bool queueSx(bool* x);
  bool queueCx(Genome* x);
  void clearSx_q();
  void clearCx_q();
  bool call(const int& r, const int& noise=0);
  void readEnv(bool*& s);
  void execute(int& z_I, int& z_dec);
  void executeDamage(int& z_I, int& z_dec);
  void executeDamage(int& z_I, int& z_dec, bool*& z);
  int evaluate(int matrix[2][2]);
  int evaluateDamage(int matrix[2][2]);
  void die();
  void breed(const int& pmut, const int& max);
  void sumErrorBin(int matrix[2][2][2]) const;
  void stGenome(const int& rule);
  void stGenome(const bool& s);
  const int& gtCid() const;
  const int& gtR() const;
  const int& gtK() const;
  const bool& gtS() const;
  const bool& gtM() const;
  const int& gtF() const;
  const int& gtCumPayoff() const;
  const Genome& gtG() const;
  void stS(const bool& s);
  void stM(const bool& m);
  void stF(const int& f);
  void stG(const Genome& g);
  void updS();
  void updM();
  void updF();
  void updG();
  void pickS();
  void pickG();
  void print(const bool& print_g=false) const;

  Cell& operator = (const Cell& C);
};

#endif
