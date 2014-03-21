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

