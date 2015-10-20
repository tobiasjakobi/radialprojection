/*  radialprojection - tools to numerically compute the radial projection of point sets
 *  Copyright (C) 2012-2014 - Tobias Jakobi <tjakobi at math dot uni dash bielefeld dot de>
 *
 *  radialprojection is free software: you can redistribute it and/or modify it under the terms
 *  of the GNU General Public License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  radialprojection is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with radialprojection.
 *  If not, see <http://www.gnu.org/licenses/>.
 */

#include "common.h"

#include <sstream>

void print_usage() {
  cerr << "histogram: usage:" << endl;

  cerr << "histogram --normal: selects normal main mode" << endl;
  cerr << "histogram --largedata: selects large data main mode" << endl;
  cerr << "(both modes share the same set of parameters" << endl;

  cerr << "\tparameter 1: mode (0 = bulk; 1 = tail)" << endl;
  cerr << "\tparameter 2: step size (binning resolution)" << endl;
  cerr << "\tparameter 3: range (for 'bulk') / start (for 'tail')" << endl;

  cerr << endl;

  cerr << "histogram --2d: selects two-dimensional main mode" << endl;
  cerr << "(this mode only supports bulk mode)" << endl;

  cerr << "\tparameter 1: step size in x-direction" << endl;
  cerr << "\tparameter 2: step size in y-direction" << endl;
  cerr << "\tparameter 3: range in x-direction" << endl;
  cerr << "\tparameter 4: range in y-direction" << endl;
}

void print_info(bool bulk) {
  cerr << (bulk ?
    "info: creating envelope for bulk section (output data type is 64-bit IEEE float, 8 byte alignment).\n" :
    "info: creating envelope for tail section (output data type is 80-bit IEEE float, 16 byte alignment).\n");
}

int main_1d(int argc, char* argv[], bool largedata) {
  stringstream parser;

  uint mode = 0;
  double stepsize = -1.0;
  double param = 3.0;

  if (argc >= 2) {
    // mode
    parser.str(argv[1]);
    parser.clear();
    parser >> mode;

    if (argc >= 3) {
      // stepsize
      parser.str(argv[2]);
      parser.clear();
      parser >> stepsize;

      if (argc >= 4) {
        // param (either 'range' or 'start' depending on mode)
        parser.str(argv[3]);
        parser.clear();
        parser >> param;
      }
    }
  }

  if (mode > 1) {
    cerr << "error: unsupported mode selected.\n";
    return 1;
  }

  if (stepsize <= 0.0) {
    if (stepsize == -1.0) {
      stepsize = (mode == 0) ? 0.01 : 0.1;
    } else {
      cerr << "error: step size has to be strictly positive.\n";
      return 2;
    }
  }

  if (param <= 0.0) {
    cerr << "error: range / start has to be strictly positive.\n";
    return 3;
  }

  if (largedata) {
    print_info(mode == 0);

    if (mode == 0)
      Common::histogramEnvelope(0.0, param, stepsize, true);
    else
      Common::histoTailEnvelope(param, stepsize, false);
  } else {
    cerr << "info: processing large data (streaming).\n";
    print_info(mode == 0);

    if (mode == 0)
      Common::histogramEnvelopeLD(0.0, param, stepsize);
    else
      Common::histoTailEnvelopeLD(param, stepsize);
  }

  return 0;
}

int main_2d(int argc, char* argv[]) {
  stringstream parser;

  // The defaults create a 30x30 binning grid
  vec2d stepsize(0.1, 0.1);
  vec2d range(3.0, 3.0);

  if (argc >= 2) {
    // x-step
    parser.str(argv[1]);
    parser.clear();
    parser >> stepsize[0];

    if (argc >= 3) {
      // y-step
      parser.str(argv[2]);
      parser.clear();
      parser >> stepsize[1];

      if (argc >= 4) {
        // x-range
        parser.str(argv[3]);
        parser.clear();
        parser >> range[0];

        // y-range
        if (argc >= 5) {
          parser.str(argv[4]);
          parser.clear();
          parser >> range[1];
        }
      }
    }
  }

  if (stepsize[0] <= 0.0 || stepsize[1] <= 0.0) {
    cerr << "error: {x,y}-stepsize has to be strictly positive.\n";
    return 1;
  }

  if (range[0] <= 0.0 || range[1] <= 0.0) {
    cerr << "error: {x,y}-range has to be strictly positive.\n";
    return 2;
  }

  print_info(true);

  Common::histogramEnvelope2D(vec2d(0.0, 0.0), range, stepsize, false);

  return 0;
}

int main(int argc, char* argv[]) {
  stringstream parser;
  string tempstr;

  int main_mode = -1;

  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> tempstr;
  }

  if (tempstr == "--normal")
    main_mode = 0;
  else if (tempstr == "--largedata")
    main_mode = 1;
  else if (tempstr == "--2d")
    main_mode = 2;

  if (main_mode == -1) {
    print_usage();
    return 0;
  }

  cerr << "info: histogram envelope generation mode.\n"
       << "(stdin = histogram data / stdout = envelope data)\n";

  switch (main_mode) {
  case 0:
  case 1:
    return main_1d(argc - 1, argv + 1, main_mode == 1);

  case 2:
    return main_2d(argc - 1, argv + 1);

  default:
    assert(false);
    return 0;
  }
}
