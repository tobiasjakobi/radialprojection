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

  cerr << "info: creating vertices in visible/random mode.\n";

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

  cerr << "info: constructed patch of " << vertices.size()
       << " visible elements.\n";

  srandExt();

  out.clear();
  out.reserve(lround(double(vertices.size()) * (1.0 - prob)));

  for (vec2ilist::const_iterator i = vertices.begin(); i != vertices.end(); ++i) {
    const double p = double(rand()) * rndNorm;
    if (p >= prob) out.push_back(*i);
  }

  cerr << "info: after randomization " << out.size() << " elements remain.\n";
}

void RandomVis::vRandomVisible(uint steps, double prob, Common::vec2ilist& out) {
  using namespace Common;

  cerr << "info: creating vertices in random/visible mode.\n";

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

  cerr << "info: constructed random patch of " << vertices.size()
       << " elements.\n";

  vec2ielist ext;

  for (vec2ilist::const_iterator i = vertices.begin(); i != vertices.end(); ++i)
    ext.push_back(*i);

  sort(ext.begin(), ext.end());
  normalize(ext);
  ext.erase(unique(ext.begin(), ext.end()), ext.end());

  cerr << "info: after visibility computation " << ext.size()
       << " elements remain.\n";

  out.clear();
  out.reserve(ext.size());

  for (vec2ielist::const_iterator i = ext.begin(); i != ext.end(); ++i)
    out.push_back(*i);
}

void RandomVis::radialProjVisRnd(uint steps, double prob, Common::dlist& out) {
  using namespace Common;

  cerr << "info: creating vertices in visible/random mode.\n";

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

  cerr << "info: constructed patch of " << vertices.size()
       << " visible elements.\n";

  srandExt();

  dlist angles;
  double meandist;

  angles.reserve(lround(double(vertices.size()) * (1.0 - prob)));
  for (vec2ilist::const_iterator i = vertices.begin(); i != vertices.end(); ++i) {
    const double p = double(rand()) * rndNorm;
    if (p >= prob) angles.push_back(i->angle());
  }

  cerr << "info: applying radial projection to " << angles.size() << " vertices.\n";

  sort(angles.begin(), angles.end());

  out.clear();
  out.reserve(angles.size() - 1);
  neighbourDiff(angles, out, meandist);
  normalizeAngDists(out, meandist);

  cerr << "info: mean distance = " << meandist << endl;
}

void RandomVis::radialProjRndVis(uint steps, double prob, Common::dlist& out) {
  using namespace Common;

  cerr << "info: creating vertices in random/visible mode.\n";

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

  cerr << "info: constructed random patch of " << vertices.size()
       << " elements.\n";

  vec2ilist visible;
  visible.reserve(vertices.size());

  for (vec2ilist::const_iterator i = vertices.begin(); i != vertices.end(); ++i)
    visible.push_back(i->primitive());

  vertices.resize(0);

  sort(visible.begin(), visible.end());
  visible.erase(unique(visible.begin(), visible.end()), visible.end());

  cerr << "info: applying radial projection to " << visible.size() << " vertices.\n";

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

  cerr << "info: mean distance = " << meandist << endl;
}

int main_normal(int argc, char* argv[]) {
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

  const double probeps = 0.001;

  if ((prob < probeps) || (prob > 1.0 - probeps)) {
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

int main_statistics(int argc, char* argv[]) {
  stringstream parser;

  uint mode = 0;
  uint steps = 100;
  double probstep = 0.1;

  RandomVis::radialfunc rfunc;

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

  if (mode == 0) {
    rfunc = RandomVis::radialProjVisRnd;
  } else if (mode == 1) {
    rfunc = RandomVis::radialProjRndVis;
  } else {
    cerr << "error: unknown mode (" << mode <<  ") selected.\n";
    return 1;
  }

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
    rfunc(uint(double(steps) / sqrt(1 - prob)), prob, spacings);

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
  stringstream parser;

  // parameters for the radial projection part
  uint mode = 0;
  uint steps = 100;
  double prob = 0.5;

  // parameters for the histogram/binning part
  uint hmode = 0;
  double hstep = 0.002;
  double hparam = 1.4;

  RandomVis::radialfunc rfunc;

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

  if (argc >= 5) {
    parser.str(argv[4]);
    parser.clear();
    parser >> hmode;

    if (argc >= 6) {
      parser.str(argv[5]);
      parser.clear();
      parser >> hstep;

      if (argc >= 7) {
        parser.str(argv[6]);
        parser.clear();
        parser >> hparam;
      }
    }
  }

  if (mode == 0) {
    rfunc = RandomVis::radialProjVisRnd;
  } else if (mode == 1) {
    rfunc = RandomVis::radialProjRndVis;
  } else {
    cerr << "error: unknown mode (" << mode <<  ") selected.\n";
    return 1;
  }

  const double probeps = 0.001;

  if ((prob < probeps) || (prob > 1.0 - probeps)) {
    cerr << "error: probability value " << prob << " not in [0,1].\n";
    return 1;
  }

  if ((hmode != 0) && (hmode != 1)) {
    cerr << "error: unknown histogram mode (" << hmode <<  ") selected.\n";
    return 1;
  }

  if ((hstep <= 0.0) || (hparam <= 0.0)) {
    cerr << "error: histogram step/parameter not strictly greater than zero.\n";
    return 1;
  }

  Common::dlist spacings;

  rfunc(uint(double(steps) / sqrt(1 - prob)), prob, spacings);

  Common::BinningData binData;
  binData.range[0] = (hmode == 0 ? 0.0 : hparam);
  binData.range[1] = (hmode == 0 ? hparam : 0.0);
  binData.step = hstep;
  binData.tail = (hmode == 1);
  histogramBinning(spacings, binData);

  Common::dlist envelopeData;
  histogramScale(binData, envelopeData, 1.0 / (double(spacings.size()) * hstep));

  Common::writeRawConsole(envelopeData);

  return 0;
}

void print_usage() {
  cerr << "random: usage:" << endl;

  cerr << "random --normal: selects normal main mode" << endl;
  cerr << "\tparameter 1: mode (even = point set; odd = radial projection)" << endl;
    cerr << "\t\t" << "0/1 = visible-random; 2/3 = random-visible;" << endl;
  cerr << "\tparameter 2: steps" << endl;
  cerr << "\tparameter 3: discard probability" << endl;

  cerr << "random --statistics: selects statistics main mode" << endl;
  cerr << "\tparameter 1: mode (0 = visible-random; 1 = random-visible)" << endl;
  cerr << "\tparameter 2: range (number of vertices depends on mode)" << endl;
  cerr << "\tparameter 2: probability step (how fine [0,1] is sampled)" << endl;

  cerr << "random --single: selects single main mode" << endl;
  cerr << "(creates a histogram for a single random realisation)" << endl;
  cerr << "\tparameter 1: mode (0 = visible-random; 1 = random-visible)" << endl;
  cerr << "\tparameter 2: range (number of vertices depends on mode)" << endl;
  cerr << "\t\t(internally the range is expanded to compensate "
       << "for the loss of vertices due to the discarding process)" << endl;
  cerr << "\tparameter 3: discard probability" << endl;
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
