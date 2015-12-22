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

#include "cyclotomic_radial.h"

#include "cyclotomic_octagonal.h"
#include "cyclotomic_decagonal.h"
#include "cyclotomic_dodecagonal.h"
#include "cyclotomic_rhombic_penrose.h"

#include <sstream>



void SingleMachine::apply_shift(uint mode) {
  using namespace Common;

  // non-singular octagonal tiling doesn't need a window shift
  if (mode == octagonal_tiling || mode == octagonal_radprj) {
    return;
  }

  // The pentagon windows in the rhombic Penrose case
  // always need to be shifted into a generic position:
  if (mode == rhmbpenrose_tiling || mode == rhmbpenrose_radprj) {
    vec4i::shift.set(1.0e-4, 1.0e-4);
    return;
  }

  // In the pure cyclotomic case the window can be replaced by a circle.
  // This of course affects the tiling, but has only negligible impact on the radial projection.
  if (circularWindow) {
    cerr << "Using (simplified) circular window with same area.\n";
    return;
  }

  // When using the decagon/dodecagon window from the book, shift it slightly.
  // A singular patch can be constructed by omitting the shift.
  if (windowBookOrientation) {
    cerr << "Using decagon/dodecagon window orientation from the book, applying slight shift to it.\n";
    vec4i::shift.set(1.0e-4, 1.0e-4);
  }
}

int SingleMachine::main(int argc, char* argv[]) {
  const vec4i init(0, 0, 0, 0);

  stringstream parser;

  // inputs
  uint steps = 40;
  bool sector = false;
  uint mode = 0;
  vec4i origin(0, 0, 0, 0);

  bool force_nonlocal_test = false;
  bool second_order = false;
  bool zero_origin;

  // outputs
  Common::vec4ilist tiling, visible;
  Common::dlist spacings;
  double mean;

  if (argc >= 2) {
    const string arg(argv[1]);

    if (arg == "--second-order") {
      second_order = true;

      argc--;
      argv++;
    }
  }

  /* create dodecagonal tiling with edges
  {
    vec4i::shift.set(1.0e-4, 1.0e-4);
    Common::edgelist edges;
    DodecagonalRadial::projTilingPatch(init, 15, tiling, edges);
    cout << '{' << tiling << ',' << edges << '}' << endl;
    return 0;
  }
  */

  /*
   * Special mode to force usage of the non-local test.
   * Mainly useful for debugging the visibility algos.
   */
  if (argc == 5) {
    const string arg(argv[4]);
    if (arg == "nonlocal")
      force_nonlocal_test = true;
  }

  if (argc > 8) argc = 8;
  if (argc < 8) argc = min(argc, 4);

  switch (argc) {
    default:
      for (uint k = 0; k < 4; ++k) {
        parser.str(argv[4 + k]);
        parser.clear();
        parser >> origin[k];
      }

    /* fall-through is intended in all cases */
    case 4:
      parser.str(argv[3]);
      parser.clear();
      parser >> sector;

    case 3:
      parser.str(argv[2]);
      parser.clear();
      parser >> steps;

    case 2:
      parser.str(argv[1]);
      parser.clear();
      parser >> mode;

    case 1:
    case 0:
      /* do nothing */
      break;
  }

  if (check_mode(mode)) {
    cerr << "error: unsupported processing mode selected!\n";
    return 1;
  }

  // Check for zero origin and sanitize 'sector' input.
  zero_origin = origin.isZero();
  if (!zero_origin && sector) {
    cerr << "warn: reduction to sector only available with zero origin.\n";
    sector = false;
  }

  apply_shift(mode);

  switch (mode) {
    case octagonal_tiling:
      if (zero_origin && !force_nonlocal_test)
        Octagonal::projTilingVisLocal(init, steps, sector, tiling, visible);
      else
        Octagonal::projTilingVis(init, origin, steps, false, tiling, visible);
    break;

    case octagonal_radprj:
      if (zero_origin && !force_nonlocal_test)
        Octagonal::projTilingVisLocal(init, steps, sector, tiling, visible);
      else
        Octagonal::projTilingVis(init, origin, steps, true, tiling, visible);

      Octagonal::radialProj(visible, spacings, mean, sector);
    break;

    case decagonal_tiling:
      if (zero_origin && !force_nonlocal_test) {
        Decagonal::projTilingVisLocal(init, steps, tiling, visible);

        if (sector) {
          Common::vec4ilist vistilSector;
          Decagonal::extractSector(visible, vistilSector);
          visible.swap(vistilSector);
          cerr << "Reduced visible tiling to a sector containing "
               << visible.size() << " vertices.\n";
        }
      } else
        Decagonal::projTilingVis(init, origin, steps, false, tiling, visible);
    break;

    case decagonal_radprj:
      if (zero_origin && !force_nonlocal_test) {
        Common::vec4ilist vistilSector;

        Decagonal::projTilingVisLocal(init, steps, tiling, visible);
        Decagonal::extractSector(visible, vistilSector);
        Decagonal::radialProj(vistilSector, spacings, mean);
      } else {
        Decagonal::projTilingVis(init, origin, steps, true, tiling, visible);
        Decagonal::radialProj(visible, spacings, mean);
      }
    break;

    case dodecagonal_tiling:
      if (zero_origin && !force_nonlocal_test) {
        Dodecagonal::projTilingVisLocal(init, steps, tiling, visible);

        if (sector) {
          Common::vec4ilist vistilSector;
          Dodecagonal::extractSector(visible, vistilSector);
          visible.swap(vistilSector);
          cerr << "Reduced visible tiling to a sector containing "
               << visible.size() << " vertices.\n";
        }
      } else
        Dodecagonal::projTilingVis(init, origin, steps, false, tiling, visible);
    break;

    case dodecagonal_radprj:
      if (zero_origin && !force_nonlocal_test) {
        Common::vec4ilist vistilSector;

        Dodecagonal::projTilingVisLocal(init, steps, tiling, visible);
        Dodecagonal::extractSector(visible, vistilSector);
        Dodecagonal::radialProj(vistilSector, spacings, mean);
      } else {
        Dodecagonal::projTilingVis(init, origin, steps, true, tiling, visible);
        Dodecagonal::radialProj(visible, spacings, mean);
      }
    break;

    case rhmbpenrose_tiling:
      RhombicPenrose::projTilingVis(init, origin, steps,
        sector ? Common::proj_tiling_onlysector : Common::proj_tiling_none,
        tiling, visible);
    break;

    case rhmbpenrose_radprj:
      RhombicPenrose::projTilingVis(init, origin, steps,
                                    Common::proj_tiling_radialprojection,
                                    tiling, visible);
      RhombicPenrose::radialProj(visible, spacings, mean);
    break;

    default:
      assert(false);
    break;
  }

  if (is_tiling_mode(mode)) {
    cerr << "Size of visible point data is around "
         << uint(double(visible.size() * sizeof(vec4i)) / 1024.0)
         << " kilobytes.\n";

    cout << visible;
  } else {
    // Don't output the tiling in radial projection mode:
    // We're only interested in the angular data, which is written
    // to the console in raw mode in this case.

    Common::meanDistanceMessage(spacings.size() + 1, mean);

    if (second_order) {
      cerr << "info: computing second-order spacings." << endl;

      Common::vec2dlist spacings2;
      spacings2.reserve(spacings.size() - 1);

      Common::secondOrderSpacings(spacings, spacings2);
      Common::writeRawConsole(spacings2);
    } else {
      Common::writeRawConsole(spacings);
    }
  }

  return 0;
}

int MultiMachine::master(int argc, char* argv[]) {
  using namespace Common;

  stringstream parser;
  vec4ilist fulltiling, tiling, origins;
  double inner, outer, srad;

  uint steps, samples;
  uint percentage; /* value in the range [0, 100] */

  const vec4i init(0, 0, 0, 0);

  if (argc != 4) return 1;

  parser.str(argv[1]);
  parser.clear();
  parser >> steps;

  parser.str(argv[2]);
  parser.clear();
  parser >> samples;

  parser.str(argv[3]);
  parser.clear();
  parser >> percentage;

  SingleMachine::apply_shift(SingleMachine::decagonal_tiling);
  Decagonal::projTiling(init, steps, fulltiling);
  Decagonal::innerOuterRadius(fulltiling, inner, outer);

  RadiusSelector::radiusSq = inner * inner;

  srad = inner * (double(percentage) / 100.0);

  selectVertices<vec4ilist, RadiusSelector>(fulltiling, tiling);
  selectOrigins<vec4ilist, Decagonal::LengthSelector>(tiling, origins, samples, srad, inner);

  cerr << "srad: " << reinterpret_double_to_ullong(srad) << endl;
  cerr << "origins: " << origins << endl;

  writeRawConsole(tiling);

  return 0;
}

int MultiMachine::slave(int argc, char* argv[]) {
  using namespace Common;

  stringstream parser;
  string token;
  ullong temp;
  uint idx = 0;

  vec4i origin;
  vec4ilist tiling;
  dlist spacings;
  double radius, mdist;

  if (argc != 3) return 1;

  parser.str(argv[1]);
  parser.clear();
  parser >> temp;
  radius = reinterpret_ullong_to_double(temp);

  parser.str(argv[2]);
  parser.clear();
  while (idx < 4) {
    std::getline(parser, token, ',');

    stringstream ss(token);
    ss >> origin[idx];
    ++idx;
  }

  readRawConsole(tiling);
  Decagonal::radialProj(tiling, origin, radius, spacings, mdist);

  writeRawConsole(spacings);

  return 0;
}

void print_usage() {
  cerr << "cyclotomic_radial: usage:" << endl;

  cerr << "cyclotomic_radial --single: selects single main mode" << endl;

  cerr << "Passing --second-order as second argument switches from first to second"
       << endl << "order spacings (this only affects the radial projection modes)."
       << endl;

  cerr << "\tparameter 1: mode (even = point set; odd = radial projection)" << endl;
    cerr << "\t\t" << "0/1 = octagonal (Ammann-Beenker)" << endl;
    cerr << "\t\t" << "2/3 = decagonal (Tübingen triangle)" << endl;
    cerr << "\t\t" << "4/5 = dodecagonal (Gähler shield)" << endl;
    cerr << "\t\t" << "6/7 = rhombic Penrose (cyclotomic multi-window)" << endl;
  cerr << "\tparameter 2: steps" << endl;
  cerr << "\tparameter 3: sector (if possible reduce point set to sector)" << endl;
  cerr << "\tparameter {4,5,6,7}: custom origin (for visibility computation)" << endl;

  cerr << endl;

  cerr << "cyclotomic_radial --multi: selects multi main mode" << endl;
  cerr << "\tparameter 1: multi type (0 = master; 1 = slave)" << endl;
  cerr << "\t\t" << "(parameters below are shifted)" << endl;

  cerr << "\tmaster mode (tiling mode currently hardcoded):" << endl;
  cerr << "\t\t" << "parameter 1: steps (of initial tiling)" << endl;
  cerr << "\t\t" << "parameter 2: samples (number of subpatches to create)" << endl;
  cerr << "\t\t" << "parameter 3: percentage (determines size of subpatch)" << endl;

  cerr << "\tslave mode:" << endl;
  cerr << "\t\t" << "parameter 1: radius (double reinterpret_casted as ulong):" << endl;
  cerr << "\t\t" << "parameter 2: origin (formatted as \"x,y,z,w\")" << endl;
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

  if (main_mode == "--single") {
    ret = SingleMachine::main(argc - 1, argv + 1);
  } else if (main_mode == "--multi") {
    uint multi_mode;

    if (argc >= 3) {
      parser.str(argv[2]);
      parser.clear();
      parser >> multi_mode;

      if (multi_mode == 0) {
        ret = MultiMachine::master(argc - 2, argv + 2);
      } else {
        if (multi_mode == 1)
          ret = MultiMachine::slave(argc - 2, argv + 2);
      }
    }
  } else {
    print_usage();
  }

  return ret;
}
