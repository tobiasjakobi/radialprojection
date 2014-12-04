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
  stringstream parser;

  const double rndNorm = 1.0 / double(RAND_MAX);

  uint steps = 100;
  double prob = 0.5;

  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> steps;

    if (argc >= 3) {
      parser.str(argv[2]);
      parser.clear();
      parser >> prob;
    }
  }

  Common::vec2ilist vertices, randomized;

  for (uint y = 0; y <= steps; ++y) {
    for (uint x = y; x <= steps; ++x) {
      if (x*x + y*y > steps*steps) continue;
      if (Coprime::gcdZFast(x, y) != 1) continue;

      vertices.push_back(vec2i(x, y));
    }
  }

  Common::srandExt();

  for (Common::vec2ilist::const_iterator i = vertices.begin(); i != vertices.end(); ++i) {
    const double p = double(rand()) * rndNorm;
    if (p >= prob) randomized.push_back(*i);
  }

  cout << randomized << endl;

  return 0;
}
