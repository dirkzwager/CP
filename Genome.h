#ifndef __GENOME_H__
#define __GENOME_H__


#include "Vector.h"

class Genome;

const static int MUTATION_LIMIT = 1000;

void fromDec(int dec, int base, bool*& arr, int wrap);
unsigned int toDec(bool* arr, const int& wrap);

class Genome;
ostream& operator << (ostream& os, const Genome& G);

class Genome {
private:
  bool* Z;
  int dmx, dmy;
public:
  Genome();
  Genome(const Genome& G);
  Genome(const Genome& G1, const Genome& G2, double pcross);
  ~Genome();

  void lookupGene(bool*& g, const int& y, int& g_dec) const;
  void crossover(const Genome& G1, const Genome& G2, double pcross);
  void mutate(const int& rate);
  void randomnize();
  void set(const int& dec);
  void set(const bool& S);
  const int& gtDmx() const;
  const int& gtDmy() const;

  Genome& operator = (const Genome& G);
  const bool& operator () (const int& i) const;
  bool& operator () (const int& i);
  friend ostream& operator << (ostream& os, const Genome& G);
};


#endif
