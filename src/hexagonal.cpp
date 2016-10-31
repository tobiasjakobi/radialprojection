/*  radialprojection - tools to numerically compute the radial projection of point sets
 *  Copyright (C) 2012-2016 - Tobias Jakobi <tjakobi at math dot uni dash bielefeld dot de>
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

#include "hexagonal.h"

#include <sstream>

enum processing_mode {
  triangular_tiling  = 0, /* triangular tiling */
  triangular_radprj  = 1,
  hexagonal_tiling   = 2, /* hexagonal tiling */
  hexagonal_radprj   = 3,
  generic_tiling     = 4, /* based on triangular tiling, but using generic lattice */
  generic_radprj     = 5,
  processing_mode_end
};

void Triangular::tiling(const vec2i& initpoint, uint maxstep,
              Common::vec2ilist& tilingpoints) {
  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  const double radius = double(maxstep) * radiusFactor;
  const int steps = maxstep;

  for (int i = -steps; i < steps; ++i) {
    for (int j = -steps; j < steps; ++ j) {
      const vec2i vertex(i, j);

      if (vertex.transTriToR2().length() <= radius)
        tilingpoints.push_back(initpoint + vertex);
    }
  }

  cerr << "Constructed patch of triangular tiling with "
       << tilingpoints.size() << " vertices.\n";
}

void Triangular::tilingVisLocal(const vec2i& initpoint, uint maxstep,
                      Common::vec2ilist& tilingpoints,
                      Common::vec2ilist& visiblepoints) {
  tilingpoints.clear();
  visiblepoints.clear();

  tilingpoints.push_back(initpoint);

  const double radius = double(maxstep) * radiusFactor;
  const int steps = maxstep;

  for (int i = -steps; i < steps; ++i) {
    for (int j = -steps; j < steps; ++ j) {
      const vec2i vertex(i, j);

      if (vertex.transTriToR2().length() <= radius) {
        tilingpoints.push_back(initpoint + vertex);

        if (!vertex.isZero() && vertex.coprime())
          visiblepoints.push_back(initpoint + vertex);
      }
    }
  }

  cerr << "Constructed patch of triangular tiling with "
       << tilingpoints.size() << " vertices and "
       << visiblepoints.size() << " visible ones.\n";
}

void Triangular::extractSector(const Common::vec2ilist& input,
                     Common::vec2ilist& output) {
  using namespace Common;

  output.clear();
  output.reserve(input.size() / 6);

  for (vec2ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec2d phys(i->transTriToR2());

    if (phys.inFirstQuadrant() && phys.inSectorL3())
      output.push_back(*i);
  }
}

void Triangular::radialProj(const Common::vec2ilist& input,
                     Common::dlist& output, double& meandist) {
  using namespace Common;

  dlist angles;
  angles.reserve(input.size());

  for (vec2ilist::const_iterator i = input.begin(); i != input.end(); ++i)
    angles.push_back(i->transTriToR2().angle());

  output.clear();
  output.reserve(input.size() - 1);

  sort(angles.begin(), angles.end());
  neighbourDiff(angles, output, meandist);
  normalizeAngDists(output, meandist);
}

void Hexagonal::tiling(const vec2i& initpoint, uint maxstep,
                     Common::vec2ilist& tilingpoints) {
  using namespace Common;

  vec2ilist helper;

  /* We first construct a triangular helper tiling, and then build *
   * honeycombs around each vertex of the tiling.                  */

  const double radius = double(2 * maxstep) * Triangular::radiusFactor;

  helper.reserve((maxstep + 1) * (maxstep + 1));

  for (int i = -int(maxstep); i < int(maxstep); ++i)
    for (int j = -int(maxstep); j < int(maxstep); ++j)
      helper.push_back(initpoint + vec2i(2 * i, 2 * j));

  cerr << "Constructed helper patch of triangular tiling with "
       << helper.size() << " vertices.\n";

  tilingpoints.clear();

  const uint numhex = 6;
  const vec2i hexsteps[6] = {
    vec2i(1, 0),  vec2i(0, 1),  vec2i(-1, 1),
    vec2i(-1, 0), vec2i(0, -1), vec2i(1, -1)
  };

  for (vec2ilist::const_iterator i = helper.begin(); i != helper.end(); ++i) {
    for (uint j = 0; j < numhex; ++j) {
      const vec2i vertex(*i + hexsteps[j]);

      if (vertex.transTriToR2().length() <= radius)
        tilingpoints.push_back(vertex);
    }
  }

  // Remove duplicate vertices
  sort(tilingpoints.begin(), tilingpoints.end());
  tilingpoints.erase(unique(tilingpoints.begin(),
                            tilingpoints.end()), tilingpoints.end());

  cerr << "Constructed patch of hexagonal tiling with "
       << tilingpoints.size() << " vertices.\n";
}

void Hexagonal::tilingVis(const vec2i& initpoint, uint maxstep,
                     Common::vec2ilist& tilingpoints,
                     Common::vec2ilist& visiblepoints) {
  using namespace Common;

  vec2ilist helper;
  vec2ielist ext;

  const double radius = double(2 * maxstep) * Triangular::radiusFactor;

  helper.reserve((maxstep + 1) * (maxstep + 1));

  for (int i = -int(maxstep); i < int(maxstep); ++i)
    for (int j = -int(maxstep); j < int(maxstep); ++ j)
      helper.push_back(initpoint + vec2i(2 * i, 2 * j));

  cerr << "Constructed helper patch of triangular tiling with "
       << helper.size() << " vertices.\n";

  const uint numhex = 6;
  const vec2i hexsteps[6] = {
    vec2i(1, 0),  vec2i(0, 1),  vec2i(-1, 1),
    vec2i(-1, 0), vec2i(0, -1), vec2i(1, -1)
  };

  tilingpoints.clear();

  for (vec2ilist::const_iterator i = helper.begin(); i != helper.end(); ++i) {
    for (uint j = 0; j < numhex; ++j) {
      const vec2i vertex(*i + hexsteps[j]);

      /* we don't need to test for zero here (zero is part of the helper patch) */
      if (vertex.transTriToR2().length() <= radius) {
        tilingpoints.push_back(vertex);
        ext.push_back(vertex);
      }
    }
  }

  // This deallocates the memory of 'helper'.
  vec2ilist().swap(helper);

  sort(ext.begin(), ext.end());
  normalize(ext);
  ext.erase(unique(ext.begin(), ext.end()), ext.end());

  visiblepoints.clear();
  visiblepoints.reserve(ext.size());

  for (vec2ielist::const_iterator i = ext.begin(); i != ext.end(); ++i)
    visiblepoints.push_back(*i);

  cerr << "Constructed patch of hexagonal tiling with "
       << tilingpoints.size() << " vertices and "
       << visiblepoints.size() << " visible ones.\n";
}

void Hexagonal::extractSector(const Common::vec2ilist& input,
                     Common::vec2ilist& output) {
  /* hexagonal and triangular tiling have the same symmetry */
  Triangular::extractSector(input, output);
}

void Hexagonal::radialProj(const Common::vec2ilist& input,
                     Common::dlist& output, double& meandist) {
  Triangular::radialProj(input, output, meandist);
}

void GenericLattice::tiling(const vec2i& initpoint, const vec2d& lattice,
                     uint maxstep, Common::vec2ilist& tilingpoints) {
  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  const double radius = double(maxstep) * sin(lattice.angle());
  const int steps = maxstep;

  for (int i = -steps; i < steps; ++i) {
    for (int j = -steps; j < steps; ++ j) {
      const vec2i vertex(i, j);

      if (vertex.transGenericToR2(lattice).length() <= radius)
        tilingpoints.push_back(initpoint + vertex);
    }
  }

  cerr << "Constructed patch of generic tiling (lattice vector = "
       << lattice << ") with "
       << tilingpoints.size() << " vertices.\n";
}

void GenericLattice::tilingVisLocal(const vec2i& initpoint, const vec2d& lattice,
                     uint maxstep, Common::vec2ilist& tilingpoints,
                     Common::vec2ilist& visiblepoints) {
  tilingpoints.clear();
  visiblepoints.clear();

  tilingpoints.push_back(initpoint);

  const double radius = double(maxstep) * sin(lattice.angle());
  const int steps = maxstep;

  for (int i = -steps; i < steps; ++i) {
    for (int j = -steps; j < steps; ++ j) {
      const vec2i vertex(i, j);

      if (vertex.transGenericToR2(lattice).length() <= radius) {
        tilingpoints.push_back(initpoint + vertex);

        if (!vertex.isZero() && vertex.coprime())
          visiblepoints.push_back(initpoint + vertex);
      }
    }
  }

  cerr << "Constructed patch of generic tiling (lattice vector = "
       << lattice << ") with "
       << tilingpoints.size() << " vertices and "
       << visiblepoints.size() << " visible ones.\n";
}

void GenericLattice::radialProj(const Common::vec2ilist& input,
                     const vec2d& lattice, Common::dlist& output,
                     double& meandist) {
  using namespace Common;

  dlist angles;
  angles.reserve(input.size());

  for (vec2ilist::const_iterator i = input.begin(); i != input.end(); ++i)
    angles.push_back(i->transGenericToR2(lattice).angle());

  output.clear();
  output.reserve(input.size() - 1);

  sort(angles.begin(), angles.end());
  neighbourDiff(angles, output, meandist);
  normalizeAngDists(output, meandist);
}


int main(int argc, char* argv[]) {
  const vec2i init(0, 0);

  stringstream parser;

  uint mode = 0;
  uint steps = 40;
  bool sector = false;
  vec2d lattice;

  Common::vec2ilist tiling, visible;
  Common::dlist output;
  double mean;

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
        parser >> sector;
      }
    }
  }

  /* Parse the (second) lattice vector for the generic tiling/lattice mode. */
  if (argc >= 6) {
    double x, y;

    parser.str(argv[4]);
    parser.clear();
    parser >> x;

    parser.str(argv[5]);
    parser.clear();
    parser >> y;

    lattice.set(x, y);
  } else {
    lattice.set(0.0, 1.0);
  }

  if (mode >= processing_mode_end) {
    cerr << "error: unsupported processing mode selected!\n";
    return 1;
  }

  switch (mode) {
    case triangular_tiling:
      Triangular::tilingVisLocal(init, steps, tiling, visible);

      if (sector) {
        Triangular::extractSector(visible, tiling);
        visible.swap(tiling);
        cerr << "Reduced visible tiling to a sector containing "
             << visible.size() << " vertices.\n";
      }
    break;

    case triangular_radprj:
      Triangular::tilingVisLocal(init, steps, tiling, visible);
      Triangular::extractSector(visible, tiling);
      Triangular::radialProj(tiling, output, mean);
    break;

    case hexagonal_tiling:
      Hexagonal::tilingVis(init, steps, tiling, visible);

      if (sector) {
        Hexagonal::extractSector(visible, tiling);
        visible.swap(tiling);
        cerr << "Reduced visible tiling to a sector containing "
             << visible.size() << " vertices.\n";
      }
    break;

    case hexagonal_radprj:
      /* radial projection is the same as for the triangular tiling, since *
       * the missing points (the center of the hexagons) are all invisible *
       * points (they all have coordinates of the type 2 * (i,j)).         */
      Hexagonal::tilingVis(init, steps, tiling, visible);
      Hexagonal::extractSector(visible, tiling);
      Hexagonal::radialProj(tiling, output, mean);
    break;

    case generic_tiling:
      /* The 'sector' parameter is ignored here. */
      GenericLattice::tilingVisLocal(init, lattice, steps, tiling, visible);
    break;

    case generic_radprj:
      GenericLattice::tilingVisLocal(init, lattice, steps, tiling, visible);
      GenericLattice::radialProj(visible, lattice, output, mean);
    break;

    default:
      assert(false);
    break;
  }

  if (mode % 2 == 0) {
    cerr << "Size of visible point data is around "
         << uint(double(visible.size() * sizeof(vec2i)) / 1024.0)
         << " kilobytes.\n";

    cout << visible;
  } else {
    Common::meanDistanceMessage(output.size() + 1, mean);
    Common::writeRawConsole(output);
  }

  return 0;
}

