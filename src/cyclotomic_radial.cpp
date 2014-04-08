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

  // The pentagon windows in the rhombic Penrose case
  // always need to be shifted into a generic position:
  if (mode == 6 || mode == 7) {
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
  if (windowBookOrientation && (mode >= 2 && mode <= 5)) {
    cerr << "Using decagon/dodecagon window orientation from the book, applying slight shift to it.\n";
    vec4i::shift.set(1.0e-4, 1.0e-4);
  }
}

int SingleMachine::main(int argc, char* argv[]) {
  const vec4i init(0, 0, 0, 0);
  uint steps = 40;
  bool sector = false;
  uint mode = 0;
  Common::vec4ilist tiling, visible;

  vec4i origin;
  bool use_default_origin = true;

  /* create dodecagonal tiling with edges
  {
    vec4i::shift.set(1.0e-4, 1.0e-4);
    Common::edgelist edges;
    DodecagonalRadial::projTilingPatch(init, 15, tiling, edges);
    cout << '{' << tiling << ',' << edges << '}' << endl;
    return 0;
  }
  */

  if (argc >= 2) {
    stringstream ss(argv[1]);
    ss >> steps;
  }

  if (argc >= 3) {
    stringstream ss(argv[2]);
    ss >> mode;
  }

  if (argc >= 4) {
    stringstream ss(argv[3]);
    ss >> sector;
  }

  // Parse (non-default) origin. Currently only the octogonal code uses this!
  if (argc >=8) {
    for (uint k = 0; k < 4; ++k) {
      stringstream ss(argv[4 + k]);
      ss >> origin[k];
    }

    if (origin != vec4i(0, 0, 0, 0)) use_default_origin = false;
  }

  apply_shift(mode);

  /* mode 0 and 1 : octogonal tiling (L8 lattice)           *
   * mode 2 and 3 : decagonal tiling (L5 lattice)           *
   * mode 4 and 5 : dodecagonal tiling (L12 lattice)        *
   * mode 6 and 7 : rhombic Penrose tiling (L5 + 4 windows) *
   * Even modes create vertex data in the Z4 lattice space. *
   * Uneven modes are radial projection modes and           *
   * output double-precision floating point data.           */
  switch (mode) {
    case 0: if (use_default_origin)
              Octogonal::projTilingVisLocal(init, steps, sector, tiling, visible);
            else
              Octogonal::projTilingVis(init, origin, steps, false, tiling, visible); // onlySector is ignored
            break;

    case 1: if (use_default_origin)
              Octogonal::projTilingVisLocal(init, steps, sector, tiling, visible);
            else
              Octogonal::projTilingVis(init, origin, steps, true, tiling, visible); // onlySector is ignored
            {
              Common::dlist output;
              double mean;

              Octogonal::radialProj(visible, output, mean, sector);

              cerr << "mean distance " << mean
                   << " during radial projection of " << (output.size() + 1)
                   << " vertices.\n";

              Common::writeRawConsole(output);
            }
            break;

    case 2: Decagonal::projTilingVisLocal(init, steps, tiling, visible);
            if (sector) {
              Common::vec4ilist visSector;
              Decagonal::extractSector(visible, visSector);
              visible = visSector;
              cerr << "Reduced visible tiling to a sector containing "
                   << visible.size() << " vertices.\n";
            }
            break;

    case 3: Decagonal::projTilingVisLocal(init, steps, tiling, visible);
            {
              Common::vec4ilist visSector;
              Common::dlist output;
              double mean;

              Decagonal::extractSector(visible, visSector);
              Decagonal::radialProj(visSector, output, mean);

              cerr << "mean distance " << mean
                   << " during radial projection of " << visSector.size()
                   << " vertices.\n";

              Common::writeRawConsole(output);
            }
            break;

    case 4: {
              bool robustVisibilityTest = true;
              if (argc >= 5) {
                stringstream ss(argv[4]);
                string vismode;
                ss >> vismode;

                if (vismode == "local") robustVisibilityTest = false;
              }
              
              if (robustVisibilityTest)
                Dodecagonal::projTilingVis(init, steps, false, tiling, visible);
              else
                Dodecagonal::projTilingVisLocal(init, steps, tiling, visible);
            }
            if (sector) {
              Common::vec4ilist visSector;
              Dodecagonal::extractSector(visible, visSector);
              visible = visSector;
              cerr << "Reduced visible tiling to a sector containing "
                   << visible.size() << " vertices.\n";
            }
            break;

    case 5: Dodecagonal::projTilingVis(init, steps, true, tiling, visible);
            {
              Common::dlist output;
              double mean;

              Dodecagonal::radialProj(visible, output, mean);

              cerr << "mean distance " << mean
                   << " during radial projection of " << visible.size()
                   << " vertices.\n";

              Common::writeRawConsole(output);
            }
            break;

    case 6: RhombicPenrose::projTilingAll(init, steps, tiling);
            if (sector) {
              Common::vec4ilist tilingSector;
              RhombicPenrose::extractSector(tiling, tilingSector);
              tiling.swap(tilingSector);
              cerr << "Reduced tiling to a sector containing "
                   << tiling.size() << " vertices.\n";
            }
            RhombicPenrose::selectVisible(tiling, visible, false);
            break;

    case 7: RhombicPenrose::projTilingAll(init, steps, tiling);
            {
              Common::vec4ilist tilingSector;
              Common::dlist output;
              double mean;

              RhombicPenrose::extractSector(tiling, tilingSector);
              RhombicPenrose::selectVisible(tilingSector, visible, true);

              RhombicPenrose::radialProj(visible, output, mean);

              cerr << "mean distance " << mean
                   << " during radial projection of " << (output.size() + 1)
                   << " vertices.\n";

              Common::writeRawConsole(output);
            }
            break;

    default: cerr << "error: unsupported mode selected!\n";
             return 0;
  }

  // Don't output the tiling in radial projection mode:
  // We're only interested in the angular data, which is written
  // to the console in raw mode in this case.
  if (mode % 2 == 0) {
    cerr << "Size of visible point data is around "
         << uint(double(visible.size() * sizeof(vec4i)) / 1024.0)
         << " kilobytes.\n";

    cout << visible;
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

  Octogonal::projTiling(init, steps, fulltiling);
  Octogonal::innerOuterRadius(fulltiling, inner, outer);

  RadiusSelector::radiusSq = inner * inner;

  srad = inner * (double(percentage) / 100.0);

  selectVertices<vec4i, RadiusSelector>(fulltiling, tiling);
  selectOrigins<vec4i, Octogonal::LengthSelector>(tiling, origins, samples, srad, inner);

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
  dlist output;
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
  Octogonal::radialProj(tiling, origin, radius, output, mdist);

  writeRawConsole(output);

  return 0;
}

int main(int argc, char* argv[]) {

  return SingleMachine::main(argc, argv);

  /*if (argc >= 2) {
    stringstream parser(argv[1]);
    uint mode;
    parser >> mode;

    if (mode == 0) {
      return MultiMachine::master(argc - 1, argv + 1);
    } else {
      if (mode == 1) {
        return MultiMachine::slave(argc - 1, argv + 1);
      }
    }
  }*/
}
