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
#include <fstream>
#include <algorithm>

typedef unsigned long int ulong;

void highprecision_poisson(uint steps) {
  const uint n = steps * steps;

  vector<__float128> output;

  {
    vector<__float128> angles;
    angles.resize(n);

    {
      ulong* mem = new ulong[n];
      __float128* ptr = &(angles[0]);
      const __float128 scaler = __float128(n) / __float128(numeric_limits<ulong>::max());

      ifstream urandom("/dev/urandom", ios::in|ios::binary);
      urandom.read(reinterpret_cast<char*>(mem), n * sizeof(ulong));
      urandom.close();  

      for (uint i = 0; i < n; ++i) {
        ptr[i] = __float128(mem[i]) * scaler;
      }
      delete [] mem;
    }

    sort(angles.begin(), angles.end());

    output.reserve(n - 1);

    vector<__float128>::const_iterator j = angles.begin();
    __float128 mean = *j;

    vector<__float128>::const_iterator k = j;
    ++k;
    while (k != angles.end()) {
      output.push_back(*k - *j);
      ++j;
      ++k;
    }

    --j;
    mean = (*j - mean) / __float128(angles.size() - 1);

    cerr << "mean = " << double(mean) << endl;
  }

  uint num_bin = uint(4.0 / 0.01);
  uint in_bin = 0;

  uint* bins = new uint[num_bin];
  for (uint i = 0; i < num_bin; ++i)
    bins[i] = 0;

  for (vector<__float128>::const_iterator j = output.begin(); j != output.end(); ++j) {
    const __float128 cur = *j;

    if (cur < 0.0W || cur >= 4.0W) continue;

    ++bins[uint(cur / 0.01W)];
    ++in_bin;
  }

  Common::dlist envelopeData;

  cerr << "Computing histogram with " << num_bin << " bins (interval = ["
       << 0.0 << ',' << 4.0 << "); step width = " << 0.01 << ")\n";
  cerr << "statistics: " << in_bin << " data points (from " << output.size()
         << ") fall into the binning area\n";

  const double scaler = 1.0 / (double(output.size()) * 0.01);
  for (uint i = 0; i < num_bin; ++i) {
    envelopeData.push_back(double(bins[i]) * scaler);
  }
  delete [] bins;

  Common::writeRawConsole(envelopeData);
}

int main(int argc, char* argv[]) {
  uint steps = 100;
  uint mode = 0;

  if (argc >= 2) {
    stringstream ss(argv[1]);
    ss >> steps;
  }

  if (argc >= 3) {
    stringstream ss(argv[2]);
    ss >> mode;
  }

  // Integer lattice
  if (mode == 0) {
    using namespace Common;

    vector<vec2i> vertices;
    dlist angles, output;
    double meandist;

    cerr << "info: reserving space for " << uint(double(steps*steps) * 0.3)
         << " vertices.\n";
    vertices.reserve(double(steps*steps) * 0.3);

    for (uint y = 0; y <= steps; ++y) {
      for (uint x = y; x <= steps; ++x) {
        if (x == 0 && y == 0) continue;
        if (x*x + y*y > steps*steps) continue;
        if (Coprime::gcdZ(x, y) > 1) continue;

        vertices.push_back(vec2i(x, y));
      }
    }

    cerr << "Constructed " << vertices.size()
         << " vertices of the visible lattice points.\n";

    angles.reserve(vertices.size());
    output.reserve(vertices.size() - 1);

    for (vector<vec2i>::const_iterator i = vertices.begin();
         i != vertices.end(); ++i) {
      angles.push_back(i->angle());
    }

    sort(angles.begin(), angles.end());
    neighbourDiff(angles, output, meandist);
    normalizeAngDists(output, meandist);

    writeRawConsole(output);
  }

  // Poisson case
  if (mode == 1) {
    using namespace Common;

    dlist angles, output;
    double meandist;

    angles.resize(steps * steps);
    output.reserve(steps * steps - 1);

    // construct n randomly distributed points in the interval [0,n]
    // this way the mean distance between two consecutive points
    // automatically becomes one (n = steps*steps here).
    random(steps * steps, double(steps * steps), &angles[0]);

    cerr << "Constructed " << angles.size()
         << " Poisson distributed vertices.\n";

    sort(angles.begin(), angles.end());
    neighbourDiff(angles, output, meandist);

    cerr << "mean distance = " << meandist << endl;

    writeRawConsole(output);
  }

  if (mode == 2) {
    highprecision_poisson(steps);
  }

  return 0;
}

