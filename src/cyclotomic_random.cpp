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

#include "cyclotomic_random.h"

#include "cyclotomic_octagonal.h"
#include "cyclotomic_decagonal.h"
#include "cyclotomic_dodecagonal.h"
#include "cyclotomic_rhombic_penrose.h"

#include <sstream>

// See SingleMachine::apply_shift (cyclotomic_radial) for comments.
void CyclotomicRandom::apply_shift(uint mode) {
  using namespace Common;

  if (mode == octagonal_visrnd || mode == octagonal_rndvis)
    return;

  if (mode == rhmbpenrose_visrnd || mode == rhmbpenrose_rndvis) {
    vec4i::shift.set(1.0e-4, 1.0e-4);
    return;
  }

  if (circularWindow) {
    cerr << "Using (simplified) circular window with same area.\n";
    return;
  }

  if (windowBookOrientation) {
    cerr << "Using decagon/dodecagon window orientation from the book, applying slight shift to it.\n";
    vec4i::shift.set(1.0e-4, 1.0e-4);
  }
}

template <typename T>
void CyclotomicRandom::randomize(const vector<T>& input,
                          vector<T>& output, double prob) {
  Common::srandExt();
  const double rndNorm = 1.0 / double(RAND_MAX);

  for (typename vector<T>::const_iterator i = input.begin();
       i != input.end(); ++i) {
    const double p = double(rand()) * rndNorm;
    if (p >= prob) output.push_back(*i);
  }
}

static inline void randomization_stats_msg(const Common::vec4ilist& l) {
  cerr << "info: after randomization " << l.size()
       << " vertices remain\n";
}

int main(int argc, char* argv[]) {
  using namespace CyclotomicRandom;

  const vec4i init(0, 0, 0, 0);

  stringstream parser;

  uint mode = 0;
  uint steps = 100;
  double prob = 0.5;

  Common::vec4ilist tiling, visible;
  Common::dlist spacings;
  double mean;

  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> mode;

    if (argc >= 3) {
      parser.str(argv[2]);
      parser.clear();
      parser >> steps;
    }

    if (argc >= 4) {
        parser.str(argv[3]);
        parser.clear();
        parser >> prob;
      }
  }

  const double probeps = 0.001;

  if ((prob < probeps) || (prob > 1.0 - probeps)) {
    cerr << "error: probability value " << prob << " not in [0,1].\n";
    return 1;
  }

  if (check_mode(mode)) {
    cerr << "error: unknown mode (" << mode <<  ") selected.\n";
    return 1;
  }

  apply_shift(mode);

  switch (mode) {
    case octagonal_visrnd:
      Octagonal::projTilingVisLocal(init, steps, false, tiling, visible);
      tiling.clear(); /* original tiling vertices are not used */
      randomize(visible, tiling, prob);
      randomization_stats_msg(tiling);
      Octagonal::radialProj(visible, spacings, mean, false);
    break;

    case octagonal_rndvis:
      Octagonal::projTiling(init, steps, tiling);
      randomize(tiling, visible, prob);
      tiling.swap(visible);
      randomization_stats_msg(tiling);
      Octagonal::extractVisible(init, true, tiling, visible);
      /* TODO: put message here */
      Octagonal::radialProj(visible, spacings, mean, false);
    break;

    case decagonal_visrnd:
      Decagonal::projTilingVisLocal(init, steps, tiling, visible);
      tiling.clear();
      randomize(visible, tiling, prob);
      randomization_stats_msg(tiling);
      Decagonal::radialProj(visible, spacings, mean);
    break;

    case decagonal_rndvis:
      Decagonal::projTiling(init, steps, tiling);
      randomize(tiling, visible, prob);
      tiling.swap(visible);
      /* TODO: implement */
      assert(false);
    break;

    case dodecagonal_visrnd:
      Dodecagonal::projTilingVisLocal(init, steps, tiling, visible);
      tiling.clear();
      randomize(visible, tiling, prob);
      randomization_stats_msg(tiling);
      Dodecagonal::radialProj(visible, spacings, mean);
    break;

    case dodecagonal_rndvis:
      Dodecagonal::projTiling(init, steps, tiling);
      randomize(tiling, visible, prob);
      tiling.swap(visible);
      /* TODO: implement */
      assert(false);
    break;

    case rhmbpenrose_visrnd:
      RhombicPenrose::projTilingVis(init, init, steps, false, tiling, visible);
      tiling.clear();
      randomize(visible, tiling, prob);
      randomization_stats_msg(tiling);
      RhombicPenrose::radialProj(visible, spacings, mean);
    break;

    case rhmbpenrose_rndvis:
      RhombicPenrose::projTiling(init, steps, tiling);
      randomize(tiling, visible, prob);
      tiling.swap(visible);
      /* TODO: implement */
      assert(false);
    break;

    default:
      assert(false);
    break;
  }

  cerr << "mean distance " << mean
       << " during radial projection of " << (spacings.size() + 1)
       << " vertices.\n";

  Common::writeRawConsole(spacings);

  return 0;
}
