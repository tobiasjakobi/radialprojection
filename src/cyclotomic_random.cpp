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

namespace CyclotomicRandom {
  const static RadialFunc octagonalRF(
    Octagonal::projTiling, Octagonal::projTilingVisLocal2,
    Octagonal::extractVisible, Octagonal::radialProj2,
    Octagonal::estimateGrowth);

  const static RadialFunc decagonalRF(
    Decagonal::projTiling, Decagonal::projTilingVisLocal,
    Decagonal::extractVisible, Decagonal::radialProj,
    Decagonal::estimateGrowth);

  const static RadialFunc dodecagonalRF(
    Dodecagonal::projTiling, Dodecagonal::projTilingVisLocal,
    Dodecagonal::extractVisible, Dodecagonal::radialProj,
    Dodecagonal::estimateGrowth);

  const static RadialFunc rhombicPenroseRF(
    RhombicPenrose::projTiling, RhombicPenrose::projTilingVis2,
    RhombicPenrose::extractVisible2, RhombicPenrose::radialProj,
    RhombicPenrose::estimateGrowth);
};

static inline void randomization_stats_msg(const Common::vec4ilist& l) {
  cerr << "info: after randomization " << l.size()
       << " vertices remain\n";
}

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
void CyclotomicRandom::randomize(const T& input, T& output, double prob) {
  Common::srandExt();
  const double rndNorm = 1.0 / double(RAND_MAX);

  for (typename T::const_iterator i = input.begin(); i != input.end(); ++i) {
    const double p = double(rand()) * rndNorm;
    if (p >= prob) output.push_back(*i);
  }
}

CyclotomicRandom::RadialFunc::RadialFunc(const RadialFunc& rf) {
  assert(false); // disallow use of copy-constructor
}

void CyclotomicRandom::RadialFunc::call(random_mode mode, uint steps,
                    double prob, Common::dlist& spacings) const {
  static const vec4i init(0, 0, 0, 0);

  Common::vec4ilist tiling, visible;
  double meandist;

  const uint target_size = uint(
    double(estimateGrowth(steps, false)) / (1.0 - prob));
  const uint steps_p = estimateGrowth(target_size, true);

  if (mode == cyclotomic_visrnd) {
    projTilingVis(init, steps_p, tiling, visible);
    tiling.clear(); /* original tiling vertices are not used */
    randomize<Common::vec4ilist>(visible, tiling, prob);
    tiling.swap(visible);
    randomization_stats_msg(visible);
  } else if (mode == cyclotomic_rndvis) {
    projTiling(init, steps_p, tiling);
    randomize<Common::vec4ilist>(tiling, visible, prob);
    tiling.swap(visible);
    randomization_stats_msg(tiling);
    extractVisible(init, true, tiling, visible);
  } else {
    assert(false);
    return;
  }

  radialProj(visible, spacings, meandist);
  Common::meanDistanceMessage(spacings.size() + 1, meandist);
}

int main_normal(int argc, char* argv[]) {
  using namespace CyclotomicRandom;

  const vec4i init(0, 0, 0, 0);

  stringstream parser;

  uint mode = 0;
  uint steps = 100;
  double prob = 0.5;

  const CyclotomicRandom::RadialFunc* rfunc;

  Common::dlist spacings;

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

  switch (mode) {
    case octagonal_visrnd:
    case octagonal_rndvis:
      rfunc = &octagonalRF;
    break;

    case decagonal_visrnd:
    case decagonal_rndvis:
      rfunc = &decagonalRF;
    break;

    case dodecagonal_visrnd:
    case dodecagonal_rndvis:
      rfunc = &dodecagonalRF;
    break;

    case rhmbpenrose_visrnd:
    case rhmbpenrose_rndvis:
      rfunc = &rhombicPenroseRF;
    break;

    default:
      assert(false);
      rfunc = NULL;
    break;
  }

  apply_shift(mode);

  const CyclotomicRandom::random_mode rfunc_mode = get_random_mode(mode);
  rfunc->call(rfunc_mode, steps, prob, spacings);

  Common::writeRawConsole(spacings);

  return 0;
}

int main_statistics(int argc, char* argv[]) {
  using namespace CyclotomicRandom;

  stringstream parser;

  uint mode = 0;
  uint steps = 100;
  double probstep = 0.1;

  const CyclotomicRandom::RadialFunc* rfunc;

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
        parser >> probstep;
      }
    }
  }

  switch (mode) {
    case octagonal_visrnd:
    case octagonal_rndvis:
      rfunc = &octagonalRF;
    break;

    case decagonal_visrnd:
    case decagonal_rndvis:
      rfunc = &decagonalRF;
    break;

    case dodecagonal_visrnd:
    case dodecagonal_rndvis:
      rfunc = &dodecagonalRF;
    break;

    case rhmbpenrose_visrnd:
    case rhmbpenrose_rndvis:
      rfunc = &rhombicPenroseRF;
    break;

    default:
      cerr << "error: unknown mode (" << mode <<  ") selected.\n";
      return 1;
    break;
  }

  apply_shift(mode);

  const CyclotomicRandom::random_mode rfunc_mode = get_random_mode(mode);

  const double probeps = 0.001;

  if ((probstep < probeps) || (probstep > 1.0 - probeps)) {
    cerr << "error: probability step value " << probstep << " not in [0,1].\n";
    return 1;
  }

  Common::dlist spacings;
  Common::BinningStats stats;

  /* We currently hardcode the binning parameters here. */
  stats.range[0] = 0.0;
  stats.range[1] = 3.0;
  stats.step = 0.002;
  stats.tail = false;

  cerr << "info: output format: {" << "discard probability, "
       << "minimum input, " << "maximum input, "
       << "position of largest bin" << "}\n";

  double prob = probstep;
  cout << '{';
  while (true) {
    rfunc->call(rfunc_mode, steps, prob, spacings);

    histogramStatistics(spacings, stats);

    cout << '{' << prob << ',' << stats.min << ',' << stats.max
         << ',' << stats.maxbin_position << '}';

    prob += probstep;
    if (prob > 1.0 - probeps)
      break;
    else
      cout << ',' << endl;
  }
  cout << '}' << endl;

  return 0;
}

int main_single(int argc, char* argv[]) {
  /* TODO: implement */
  assert(false);
  return -1;
}

void print_usage() {
  cerr << "cyclotomic_random: usage:" << endl;

  cerr << "cyclotomic_random --normal: selects normal main mode" << endl;
  cerr << "\tparameter 1: mode (even = visible-random; odd = random-visible)" << endl;
    cerr << "\t\t" << "0/1 = octagonal (Ammann-Beenker)" << endl;
    cerr << "\t\t" << "2/3 = decagonal (Tübingen triangle)" << endl;
    cerr << "\t\t" << "4/5 = dodecagonal (Gähler shield)" << endl;
    cerr << "\t\t" << "6/7 = rhombic Penrose (cyclotomic multi-window)" << endl;
  cerr << "\tparameter 2: steps (number of vertices depends on mode)" << endl;
  cerr << "\tparameter 2: probability step (how fine [0,1] is sampled)" << endl;

  cerr << "cyclotomic_random --statistics: selects statistics main mode" << endl;
  cerr << "\tparameter 1: mode (see \"normal\" mode)" << endl;
  cerr << "\tparameter 2: steps (number of vertices depends on mode)" << endl;
  cerr << "\tparameter 2: probability step (how fine [0,1] is sampled)" << endl;
}

int main(int argc, char* argv[]) {
  stringstream parser;
  string main_mode;
  int ret = 0;

  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> main_mode;
  }

  if (main_mode == "--normal")
    ret = main_normal(argc - 1, argv + 1);
  else if (main_mode == "--statistics")
    ret = main_statistics(argc - 1, argv + 1);
  else if (main_mode == "--single")
    ret = main_single(argc - 1, argv + 1);
  else
    print_usage();

  return ret;
}
