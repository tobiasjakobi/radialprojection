#include "common.h"

#include <sstream>

int main(int argc, char* argv[]) {

  uint mode = 0; // 0 = bulk, 1 = tail
  double stepsize;
  double param = 3.0; // histogram range (bulk mode), histogram start (tail mode)

  cerr << "Histogram envelope generation mode: stdin = histogram data / stdout = envelope data\n";

  // Parse mode
  if (argc >= 2) {
    stringstream ss(argv[1]);
    ss >> mode;

    if (mode != 0 && mode != 1) {
      cerr << "Unsupported mode selected (0 = bulk, 1 = tail).\n";
      return 0;
    }
  }

  // Parse step size
  if (argc >= 3) {
    stringstream ss(argv[2]);
    ss >> stepsize;

    if (stepsize < 0.0) {
      cerr << "Step size has to be a positive value.\n";
      return 0;
    }
  } else {
    stepsize = (mode == 0) ? 0.01 : 0.1;
  }

  // Parse (range/start) parameter
  if (argc >= 4) {
    stringstream ss(argv[3]);
    ss >> param;

    if (param <= 0.0) {
      cerr << "Range/start parameter has to be strictly positive.\n";
      return 0;
    }
  }

  if (mode == 0) {
    cerr << "Creating envelope for bulk section (output data type is 64-bit IEEE float, 8 byte alignment).\n";
    Common::histogramEnvelope(0.0, param, stepsize);
  } else {
    cerr << "Creating envelope for tail section (output data type is 80-bit IEEE float, 16 byte alignment).\n";
    Common::histoTailEnvelope(param, stepsize);
  }

  return 0;
}

