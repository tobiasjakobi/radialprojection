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

  for (int i = -steps + 1; i < steps; ++i) {
    for (int j = -steps + 1; j < steps; ++ j) {
      const vec2i vertex(i, j);

      if (vertex.transTriToR2().length() > radius)
        continue;

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

  for (int i = -steps + 1; i < steps; ++i) {
    for (int j = -steps + 1; j < steps; ++ j) {
      const vec2i vertex(i, j);

      if (vertex.transTriToR2().length() > radius)
        continue;

      tilingpoints.push_back(initpoint + vertex);

      if (!vertex.isZero() && vertex.coprime())
        visiblepoints.push_back(initpoint + vertex);
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

  output.clear();
  output.reserve(input.size());

  dlist angles;
  angles.reserve(input.size());

  for (vec2ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec2d phys(i->transTriToR2());
    angles.push_back(phys.angle());
  }

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
  const int steps = maxstep;

  helper.push_back(initpoint);

  for (int i = -steps + 1; i < steps; ++i) {
    for (int j = -steps + 1; j < steps; ++ j) {
      const vec2i vertex(2 * i, 2 * j);

      if (vertex.transTriToR2().length() > radius)
        continue;

      helper.push_back(initpoint + vertex);
    }
  }

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

      if (find(tilingpoints.rbegin(), tilingpoints.rend(), vertex) == tilingpoints.rend()) {
        tilingpoints.push_back(vertex);
      }
    }
  }

  cerr << "Constructed patch of hexagonal tiling with "
       << tilingpoints.size() << " vertices.\n";
}

void Hexagonal::tilingVis(const vec2i& initpoint, uint maxstep,
                     Common::vec2ilist& tilingpoints,
                     Common::vec2ilist& visiblepoints) {
  // TODO: implement
}

void Hexagonal::extractSector(const Common::vec2ilist& input,
                     Common::vec2ilist& output) {
  // TODO: implement
}

void Hexagonal::radialProj(const Common::vec2ilist& input,
                     Common::dlist& output, double& meandist) {
  // TODO: implement
}

void GenericLattice::tiling(const vec2i& initpoint, const vec2d& lattice,
                     uint maxstep, Common::vec2ilist& tilingpoints) {
  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  const double radius = double(maxstep) * sin(lattice.angle());
  const int steps = maxstep;

  for (int i = -steps + 1; i < steps; ++i) {
    for (int j = -steps + 1; j < steps; ++ j) {
      const vec2i vertex(i, j);

      if (vertex.transGenericToR2(lattice).length() > radius)
        continue;

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

  for (int i = -steps + 1; i < steps; ++i) {
    for (int j = -steps + 1; j < steps; ++ j) {
      const vec2i vertex(i, j);

      if (vertex.transGenericToR2(lattice).length() > radius)
        continue;

      tilingpoints.push_back(initpoint + vertex);

      if (!vertex.isZero() && vertex.coprime())
        visiblepoints.push_back(initpoint + vertex);
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

  output.clear();
  output.reserve(input.size());

  dlist angles;
  angles.reserve(input.size());

  for (vec2ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec2d phys(i->transGenericToR2(lattice));
    angles.push_back(phys.angle());
  }

  sort(angles.begin(), angles.end());
  neighbourDiff(angles, output, meandist);
  normalizeAngDists(output, meandist);
}


int main(int argc, char* argv[]) {
  const vec2i init(0, 0);

  uint steps = 40;
  uint mode = 0;
  bool sector = false;
  vec2d lattice;

  Common::vec2ilist tiling, visible;
  Common::dlist output;
  double mean;

  bool output_vertices = true;

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

  /* Parse the (second) lattice vector for the generic tiling/lattice mode. */
  if (argc >= 6) {
    double x, y;
    stringstream ss_x(argv[4]), ss_y(argv[5]);

    ss_x >> x;
    ss_y >> y;
    lattice.set(x, y);
  } else {
    lattice.set(0.0, 1.0);
  }

  switch (mode) {
    case triangular_tiling:
    {
      Triangular::tilingVisLocal(init, steps, tiling, visible);

      if (sector) {
        Common::vec2ilist vistilSector;
        Triangular::extractSector(visible, vistilSector);
        visible.swap(vistilSector);
        cerr << "Reduced visible tiling to a sector containing "
             << visible.size() << " vertices.\n";
      }
    }
    break;

    case triangular_radprj:
    {
      Common::vec2ilist vistilSector;

      Triangular::tilingVisLocal(init, steps, tiling, visible);
      Triangular::extractSector(visible, vistilSector);
      Triangular::radialProj(vistilSector, output, mean);

      output_vertices = false;
    }
    break;

    case hexagonal_tiling:
    {
      Hexagonal::tiling(init, steps, tiling);
      visible.swap(tiling); // TODO: debug

      if (sector) {
        // TODO
      }
    }
    break;

    case hexagonal_radprj:
    {
      // TODO: implement

      output_vertices = false;
    }
    break;

    case generic_tiling:
    {
      /* The 'sector' parameter is ignored here. */
      GenericLattice::tilingVisLocal(init, lattice, steps, tiling, visible);
    }
    break;

    case generic_radprj:
    {
      GenericLattice::tilingVisLocal(init, lattice, steps, tiling, visible);
      GenericLattice::radialProj(visible, lattice, output, mean);

      output_vertices = false;
    }
    break;

    default:
      cerr << "error: unsupported processing mode selected!\n";
      return 1;
  }

  if (output_vertices) {
    cerr << "Size of visible point data is around "
         << uint(double(visible.size() * sizeof(vec2i)) / 1024.0)
         << " kilobytes.\n";

    cout << visible;
  } else {
    cerr << "mean distance " << mean
         << " during radial projection of " << (output.size() + 1)
         << " vertices.\n";

    Common::writeRawConsole(output);
  }

  return 0;
}

