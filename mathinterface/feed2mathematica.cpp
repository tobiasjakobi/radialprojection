#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

using namespace std;

typedef unsigned int uint;

// Disable C++ name mangling so that exporting to mathematica works
extern "C" {
  #include "mathlink.h"

  struct vec2d {
    double a[2];
  };

  void ReadDoubleVec2Data(const char* inputfile) {
    vec2d buffer[512];
    unsigned temp;
    int dims[2];
    vector<vec2d> data;

    ifstream input(inputfile, ios::in | ios::binary);

    if (!input.is_open()) {
      cerr << "error: opening file failed\n";
      MLPutMessage(stdlink, MLAbortMessage);
      return;
    }

    cerr << "status: reading from file " << inputfile << endl;
    input.seekg(0, ios::beg);

    input.read(reinterpret_cast<char*>(&temp), sizeof(unsigned));
    if (input.eof() && input.fail()) goto readfail;
    if (temp != sizeof(double)) goto signfail;

    input.read(reinterpret_cast<char*>(&temp), sizeof(unsigned));
    if (input.eof() && input.fail()) goto readfail;
    if (temp != 2) goto signfail;

    input.read(reinterpret_cast<char*>(&temp), sizeof(unsigned));
    if (input.eof() && input.fail()) goto readfail;

    data.reserve(temp);

    while (true) {
      input.read(reinterpret_cast<char*>(buffer), sizeof(vec2d) * 512);

      if (!input) {
        const int bytesread = input.gcount();

        if (bytesread % sizeof(vec2d) != 0) {
          cerr << "error: failed to read full double vec2ds\n.";
          break;
        }

        for (uint k = 0; k < (bytesread/sizeof(vec2d)); ++k)
          data.push_back(buffer[k]);

        break;

      } else {
        for (uint k = 0; k < 512; ++k)
          data.push_back(buffer[k]);
      }
    }

    input.close();
    dims[0] = data.size();
    dims[1] = 2;

    if (data.size() != temp)
      cerr << "warning: length mismatch after reading data.\n";

    // Send data to mathematica
    MLPutReal64Array(stdlink, reinterpret_cast<double*>(&(*data.begin())), dims, 0, 2);
    return;

  signfail:
    cerr << "error: verifying signature failed\n";
    MLPutMessage(stdlink, MLAbortMessage);
    return;

  readfail:
    cerr << "error: reading signature failed\n";
    MLPutMessage(stdlink, MLAbortMessage);
    return;
  }

  void ReadDoubleData(const char* inputfile) {
    double buffer[1024];
    vector<double> data;

    ifstream input(inputfile, ios::in | ios::binary);

    if (!input.is_open()) {
      cerr << "error: opening file failed\n";
      MLPutMessage(stdlink, MLAbortMessage);
      return;
    }

    cerr << "status: reading from file " << inputfile << endl;
    input.seekg(0, ios::beg);

    while (true) {
      input.read(reinterpret_cast<char*>(buffer), sizeof(double) * 1024);

      if (!input) {
        const int bytesread = input.gcount();

        if (bytesread % sizeof(double) != 0) {
          cerr << "error: failed to read full doubles\n.";
          break;
        }

        for (uint k = 0; k < (bytesread/sizeof(double)); ++k) {
          data.push_back(buffer[k]);
        }
        break;

      } else {
        for (uint k = 0; k < 1024; ++k) {
          data.push_back(buffer[k]);
        }
      }
    }

    input.close();
    int dims[1] = {data.size()};

    // Send data to mathematica
    MLPutReal64Array(stdlink, &(*data.begin()), dims, 0, 1);
  }

  /* Read extended precision floating point data from a file:                                     *
   * Despite Mathematica calling the function MLPutReal128Array, the data passed to the kernel    *
   * is of type "long double", which corresponds to 80-bit IEEE floating point data packed        *
   * in either 96 bits (12 bytes) or 128 bits (16 bytes) (on x86-64).                             */
  void ReadExtFloatData(const char* inputfile) {
    long double buffer[1024];
    vector<long double> data;

    ifstream input(inputfile, ios::in | ios::binary);

    if (!input.is_open()) {
      cerr << "error: opening file failed\n";
      MLPutMessage(stdlink, MLAbortMessage);
      return;
    }

    cerr << "status: reading from file " << inputfile << endl;
    input.seekg(0, ios::beg);

    while (true) {
      input.read(reinterpret_cast<char*>(buffer), sizeof(long double) * 1024);

      if (!input) {
        const int bytesread = input.gcount();

        if (bytesread % sizeof(long double) != 0) {
          cerr << "error: failed to read full extended precision floats\n.";
          break;
        }

        for (uint k = 0; k < (bytesread/sizeof(long double)); ++k) {
          data.push_back(buffer[k]);
        }
        break;

      } else {
        for (uint k = 0; k < 1024; ++k) {
          data.push_back(buffer[k]);
        }
      }
    }

    input.close();
    int dims[1] = {data.size()};

    MLPutReal128Array(stdlink, &(*data.begin()), dims, 0, 1);
  }

  // Read long 64-bit integer data from file and pass to Mathematica
  void ReadLongIntData(const char* inputfile) {
    long int buffer[1024];
    vector<long int> data;

    cerr << "debug: sizeof(long int) = " << sizeof(long int) << endl;

    ifstream input(inputfile, ios::in | ios::binary);

    if (!input.is_open()) {
      cerr << "error: opening file failed\n";
      MLPutMessage(stdlink, MLAbortMessage);
      return;
    }

    cerr << "status: reading from file " << inputfile << endl;
    input.seekg(0, ios::beg);

    while (true) {
      input.read(reinterpret_cast<char*>(buffer), sizeof(long int) * 1024);

      if (!input) {
        const int bytesread = input.gcount();

        if (bytesread % sizeof(long int) != 0) {
          cerr << "error: failed to read full extended precision floats\n.";
          break;
        }

        for (uint k = 0; k < (bytesread/sizeof(long int)); ++k) {
          data.push_back(buffer[k]);
        }
        break;

      } else {
        for (uint k = 0; k < 1024; ++k) {
          data.push_back(buffer[k]);
        }
      }
    }

    input.close();
    int dims[1] = {data.size()};

    MLPutInteger64Array(stdlink, &(*data.begin()), dims, 0, 1);
  }
}

int main(int argc, char* argv[]) {

  return MLMain(argc, argv);
}

