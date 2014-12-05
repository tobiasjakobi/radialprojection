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

#include "random.h"

#include <sstream>

void RandomVis::vVisibleRandom(uint steps, double prob, Common::vec2ilist& out) {
  using namespace Common;

  assert((prob >= 0.0) && (prob <= 1.0));

  const double rndNorm = 1.0 / double(RAND_MAX);

  vec2ilist vertices;
  vertices.reserve(((steps + 1) * (steps + 1)) / 2);

  for (uint y = 0; y <= steps; ++y) {
    for (uint x = y; x <= steps; ++x) {
      if (x*x + y*y > steps*steps) continue;
      if (Coprime::gcdZFast(x, y) != 1) continue;

      vertices.push_back(vec2i(x, y));
    }
  }

  srandExt();

  out.clear();
  out.reserve(lround(double(vertices.size()) * (1.0 - prob)));

  for (vec2ilist::const_iterator i = vertices.begin(); i != vertices.end(); ++i) {
    const double p = double(rand()) * rndNorm;
    if (p >= prob) out.push_back(*i);
  }
}

void RandomVis::vRandomVisible(uint steps, double prob, Common::vec2ilist& out) {
  using namespace Common;

  assert((prob >= 0.0) && (prob <= 1.0));

  const double rndNorm = 1.0 / double(RAND_MAX);
  const uint est_vertices = lround(double(((steps + 1) * (steps + 1)) / 2) * (1.0 - prob));

  vec2ilist vertices;
  vertices.reserve(est_vertices);

  srandExt();

  for (uint y = 0; y <= steps; ++y) {
    for (uint x = y; x <= steps; ++x) {
      if (x == 0 && y == 0) continue;
      if (x*x + y*y > steps*steps) continue;

      const double p = double(rand()) * rndNorm;
      if (p >= prob) vertices.push_back(vec2i(x, y));
    }
  }

  vec2ielist ext;

  for (vec2ilist::const_iterator i = vertices.begin(); i != vertices.end(); ++i)
    ext.push_back(*i);

  sort(ext.begin(), ext.end());
  normalize(ext);
  ext.erase(unique(ext.begin(), ext.end()), ext.end());

  out.clear();
  out.reserve(ext.size());

  for (vec2ielist::const_iterator i = ext.begin(); i != ext.end(); ++i)
    out.push_back(*i);
}

void RandomVis::radialProjVisRnd(uint steps, double prob, Common::dlist& out) {
  using namespace Common;

  assert((prob >= 0.0) && (prob <= 1.0));

  const double rndNorm = 1.0 / double(RAND_MAX);

  vec2ilist vertices;
  vertices.reserve(((steps + 1) * (steps + 1)) / 2);

  for (uint y = 0; y <= steps; ++y) {
    for (uint x = y; x <= steps; ++x) {
      if (x*x + y*y > steps*steps) continue;
      if (Coprime::gcdZFast(x, y) != 1) continue;

      vertices.push_back(vec2i(x, y));
    }
  }

  srandExt();

  dlist angles;
  double meandist;

  angles.reserve(lround(double(vertices.size()) * (1.0 - prob)));
  for (vec2ilist::const_iterator i = vertices.begin(); i != vertices.end(); ++i) {
    const double p = double(rand()) * rndNorm;
    if (p >= prob) angles.push_back(i->angle());
  }

  sort(angles.begin(), angles.end());

  out.clear();
  out.reserve(angles.size() - 1);
  neighbourDiff(angles, out, meandist);
  normalizeAngDists(out, meandist);
}

void RandomVis::radialProjRndVis(uint steps, double prob, Common::dlist& out) {
  using namespace Common;

  assert((prob >= 0.0) && (prob <= 1.0));

  const double rndNorm = 1.0 / double(RAND_MAX);
  const uint est_vertices = lround(double(((steps + 1) * (steps + 1)) / 2) * (1.0 - prob));

  vec2ilist vertices;
  vertices.reserve(est_vertices);

  srandExt();

  for (uint y = 0; y <= steps; ++y) {
    for (uint x = y; x <= steps; ++x) {
      if (x == 0 && y == 0) continue;
      if (x*x + y*y > steps*steps) continue;

      const double p = double(rand()) * rndNorm;
      if (p >= prob) vertices.push_back(vec2i(x, y));
    }
  }

  vec2ilist visible;
  visible.reserve(vertices.size());

  for (vec2ilist::const_iterator i = vertices.begin(); i != vertices.end(); ++i)
    visible.push_back(i->primitive());

  vertices.resize(0);

  sort(visible.begin(), visible.end());
  visible.erase(unique(visible.begin(), visible.end()), visible.end());

  dlist angles;
  double meandist;

  angles.reserve(visible.size());

  for (vec2ilist::const_iterator i = visible.begin(); i != visible.end(); ++i)
    angles.push_back(i->angle());

  visible.resize(0);

  sort(angles.begin(), angles.end());

  out.clear();
  out.reserve(angles.size() - 1);
  neighbourDiff(angles, out, meandist);
  normalizeAngDists(out, meandist);
}

int main(int argc, char* argv[]) {
  stringstream parser;

  uint mode = 0;
  uint steps = 100;
  double prob = 0.5;

  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> mode;

    if (argc >= 3) {
      parser.str(argv[2]);
      parser.clear();
      parser >> steps;

      if (argc >= 4) {
        parser.str(argv[3]);
        parser.clear();
        parser >> prob;
      }
    }
  }

  Common::vec2ilist randomized;
  Common::dlist spacings;

  if ((prob < 0.0) || (prob > 1.0)) {
    cerr << "error: probability value " << prob << " not in [0,1].\n";
    return 1;
  }

  if (mode > 3) {
    cerr << "error: unknown mode (" << mode <<  ") selected.\n";
    return 1;
  }

  switch (mode) {
  case 0:
    RandomVis::vVisibleRandom(steps, prob, randomized);
  break;

  case 1:
    RandomVis::radialProjVisRnd(steps, prob, spacings);
  break;

  case 2:
    RandomVis::vRandomVisible(steps, prob, randomized);
  break;

  case 3:
    RandomVis::radialProjRndVis(steps, prob, spacings);
  break;

  default:
    assert(false);
  break;
  }

  if (mode % 2 == 0) {
    cout << randomized << endl;
  } else {
    Common::writeRawConsole(spacings);
  }

  return 0;
}
