#include "Cell.h"
#include "Wire.h"

Wire::Wire() : V_in(NULL), V_out(new bool), iput(NULL), oput(NULL), damage(0), total_payoff(0) {
}

Wire::Wire(const int& pe) : V_in(NULL), V_out(new bool), iput(NULL), oput(NULL), damage(pe), total_payoff(0) {
}

Wire::Wire(const Wire& W) : V_in(NULL), V_out(new bool), iput(NULL), oput(NULL) {
  damage = W.damage;
  total_payoff = W.total_payoff;
  V_in = W.V_in;
  *V_out = *W.V_out;
  iput = W.iput;
  oput = W.oput;
}

Wire::~Wire() {
  delete[] V_out;
  V_out = NULL;
}

void Wire::setV_in() {
  bool make_error = ((rand()%MAX_DAMAGE) < damage) ? true : false;
  oput->setV_out();
  
  if (make_error) {
    *V_in = *V_in ? false : true;
    //error_bin[*V_out][*V_in][true]++;
  }
  //else error_bin[*V_out][*V_in][false]++;
}

void Wire::setV_out() {
  *V_out = iput->gtS();
}

void Wire::stDamage(const int& pe) {
  if (pe < 0) {
    cout << "Invalid argument, please provide values from 0 to 100." << endl;
    return;
  }
  if (pe > 100) {
    cout << "Invalid argument, please provide values from 0 to 100." << endl;
    return;
  }
  damage = pe;
}

const int& Wire::gtDamage() const {
  return damage;
}

int Wire::payState(int matrix[2][2])  {
  total_payoff += matrix[*V_out][*V_in];
  return matrix[*V_out][*V_in];
}

void Wire::connect(Cell& C_in) {
  iput = &C_in;
}

void Wire::connect(Wire& W_out) {
  oput = &W_out;
  V_in = W_out.V_out;
}

Wire& Wire::operator = (const Wire& W) {
  damage = W.damage;
  total_payoff = W.total_payoff;
  V_in = W.V_in;
  if (V_out) {
    delete V_out;
    V_out=NULL;
  }
  V_out = W.V_out;
  iput = W.iput;
  oput = W.oput;
}

void Wire::print() const {
  cout << "damage: " << damage << endl
       << "total_payoff: " << total_payoff << endl;
  cout << iput->gtCid() << "->" << oput->iput->gtCid() << ":" << *V_out << " | "
       << oput->iput->gtCid() << "->" << iput->gtCid() << ":" << *V_in << "." << endl;
}

/*

  cout << "error -> ddf: " << error_bin[0][0][0] 
       << " dcf: " << error_bin[0][1][0] 
       << " cdf: " << error_bin[1][0][0] 
       << " ccf: " << error_bin[1][1][0] 
       << endl;
  cout << "error -> ddt: " << error_bin[0][0][1] 
       << " dct: " << error_bin[0][1][1] 
       << " cdt: " << error_bin[1][0][1] 
       << " cct: " << error_bin[1][1][1] 
       << endl;

*/
