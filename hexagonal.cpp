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

int main(int argc, char* argv[]) {
  using namespace Common;

  uint steps = 40;
  uint mode = 0;

  if (argc >= 2) {
    stringstream ss(argv[1]);
    ss >> steps;
  }

  if (argc >= 3) {
    stringstream ss(argv[2]);
    ss >> mode;
  }

  // hex-tiling case

  // 1st true -> only generate visible vertices
  // 2nd true -> only generate a 1/6-sector
  hexTiling tiling(steps, vec2i(0, 0), true, true); 
  const vec2ilist& in = tiling.getVertices();
  vector<vec2d> vis2d;

  // convert to 2d
  vis2d.reserve(in.size());
  for (vector<vec2i>::const_iterator i = in.begin(); i != in.end(); ++i) {
    vis2d.push_back(i->transHexTo2D());
  }

  // tri-tiling code

  /*triTiling tiling(steps, vec2i(0, 0));
  const vec2ilist& in = tiling.getVertices();
  vector<vec2d> vis2d;

  // restrict to 1/6-sector, select visible tiling points and convert to 2d
  for (vector<vec2i>::const_iterator i = in.begin(); i != in.end(); ++i) {
    if (Coprime::gcdZ(abs(i->x), abs(i->y)) == 1) {
      const vec2d t(i->transHexTo2D());

      if (t.inFirstQuadrant() &&
          t.inSectorL3()) vis2d.push_back(t);
    }
  }*/

  Common::dlist angles;
  double mean;
  Common::radialProj(vis2d, angles, mean);

  cerr << "mean distance " << mean
       << " during radial projection of " << (angles.size() + 1)
       << " vertices.\n";

  Common::writeRawConsole(angles);

  //cout << vis2d;

  return 0;
}

