#ifndef __WIRE__H__
#define __WIRE__H__

enum Mode {Vac,Opr};

const static int MAX_DAMAGE = 100;

class Cell;

class Wire {
 protected:
  bool* V_in;
  bool* V_out;
  Cell* iput;
  Wire* oput;
  int damage;
  int total_payoff;

  friend class Cell;
 public:
  Wire();
  Wire(const int& pe);
  Wire(const Wire& W);
  ~Wire();

  void setV_in();
  void setV_out();
  void stDamage(const int& pe);
  const int& gtDamage() const;
  int payState(int matrix[2][2]);
  void connect(Cell& C_in);
  void connect(Wire& W_out);
  Wire& operator = (const Wire& W);
  void print() const;
};

#endif
