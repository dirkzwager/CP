#ifndef __CA_H__
#define __CA_H__

#include "Cell.h"

class CA {
private:
  int R,K,dmx,dmy;
  Cell* model;
  int* gene_bin;
  int cum_fitness, fitness;
  bool* S;
  bool* M;
  bool** Z;

  void breed(const int& pmut, int matrix[2][2]);
  void die();
  int evaluate(int matrix[2][2]);
  void execute();
  void readCable();
public:
  CA();
  CA(const int& x, const int& y);
  CA(const CA& ca);
  ~CA();

  const int& Dmx() const;
  const int& Dmy() const;
  bool connectCells();
  void mapWires();
  void readWires();
  void breakCells(const int& bias, const int& num);
  void breakWires(const int& bias, const int& pbreak);
  static CA makeReplica(const CA& origional);
  void rndState(const int& def);
  void rndMode(const int& vac);
  void rndGenome();
  void rndGenome(const int& pdef);
  void stGenome(const string& fname);
  void stNoise(const int& noise);
  const int& gtCumFitness() const;
  const int& gtFitness() const;
  void updModel();
  void propModel(int matrix[2][2], const int& pmut);
  void damageModel(int matrix[2][2], const int& pmut);
  void coherenceModel(int matrix[2][2], const int& pmut, int bin[2][2], double avg_fit[2], double chere[2]);
  void fitnessIOerror(double matrix[5][5]);
  void stUniform(const int& rule);
  void accCells(int bin[2][2]);
  void print() const;
  void printGeneBin() const;
  void saveModel(const string& fname);
  void loadModel(const string& fname);
  CA& operator = (const CA& ca);
  const Cell& operator() (const int& x, const int& y) const;
  Cell& operator() (const int& x, const int& y);
};


#endif
