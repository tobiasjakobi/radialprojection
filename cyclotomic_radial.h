#ifndef _CYCLOTOMIC_RADIAL_H_
#define _CYCLOTOMIC_RADIAL_H_

namespace SingleMachine {

  /* Default main routine for single machine execution */
  int main(int argc, char* argv[]);

};

namespace MultiMachine {

  /* multimachine routine that is executed on master machine */
  int master(int argc, char* argv[]);

  /* multimachine routine that is executed on the slave machines */
  int slave(int argc, char* argv[]);

};

#endif

