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

#include "pdf_writer.h"

#include <cmath>
#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>

void importRawConsole(vector<ArithVisibility::bragg>& output,
                      vec2d& range, double& radius) {
  using namespace ArithVisibility;

  unsigned in_int;
  double in_flt;
  bragg data;
  unsigned read = 0;

  cin.read(reinterpret_cast<char*>(&in_int), sizeof(unsigned));
  if (cin.eof() && cin.fail()) goto readfail;
  if (in_int != sizeof(bragg)) goto signfail;

  cin.read(reinterpret_cast<char*>(&in_int), sizeof(unsigned));
  if (cin.eof() && cin.fail()) goto readfail;

  cin.read(reinterpret_cast<char*>(&in_flt), sizeof(double));
  if (cin.eof() && cin.fail()) goto readfail;
  range.x = in_flt;

  cin.read(reinterpret_cast<char*>(&in_flt), sizeof(double));
  if (cin.eof() && cin.fail()) goto readfail;
  range.y = in_flt;

  cin.read(reinterpret_cast<char*>(&in_flt), sizeof(double));
  if (cin.eof() && cin.fail()) goto readfail;
  radius = in_flt;

  output.reserve(output.size() + in_int);

  while (true) {
    cin.read(reinterpret_cast<char*>(&data), sizeof(bragg));
    if (cin.eof()) break;

    output.push_back(data);
    ++read;
  }

  if (read != in_int) {
    cerr << "warn: incomplete raw import (" << read << " of "
         << in_int << "elements)" << endl;
  }

  return;

signfail:
  cerr << "error: verifying signature failed\n";
  return ;

readfail:
  cerr << "error: reading signature failed\n";
  return ;
}

void braggToPDF(const vector<ArithVisibility::bragg>& input,
                const vec2d& range, const double radius,
                const string& filename) {
  using namespace ArithVisibility;

  cairo_surface_t *surface;
  cairo_t *cr;

  // Base width of the PDF and offset to the borders.
  const double basewidth = 800.0;
  const double offset = 50.0;

  if (input.empty()) return;

  const double height = basewidth * (range.y / range.x);

  surface = cairo_pdf_surface_create(filename.c_str(), basewidth, height);
  cairo_pdf_surface_restrict_to_version(surface, CAIRO_PDF_VERSION_1_5);
  cr = cairo_create(surface);

  cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
  cairo_set_line_width(cr, 1.0);
  cairo_set_fill_rule(cr, CAIRO_FILL_RULE_WINDING);

  // Translate origin to the middle of the page
  cairo_translate(cr, basewidth * 0.5, height * 0.5);

  // Apply scaling
  const double scaling = std::min(
    (basewidth - offset) / (2.0 * (range.x + radius)),
    (height - offset) / (2.0 * (range.y + radius)));
  cairo_scale(cr, scaling, scaling);

  for (vector<bragg>::const_iterator k = input.begin(); k != input.end(); ++k) {
    cairo_arc(cr, k->getPosition().x, k->getPosition().y, k->getIntensity(),
              0.0, 2.0 * M_PI);
    cairo_fill(cr);
    cairo_stroke(cr);
  }

  cairo_destroy(cr);
  cairo_surface_destroy(surface);
}

int main(int argc, char* argv[]) {
  using namespace ArithVisibility;

  string filename;

  vector<bragg> list;
  vec2d range;
  double rad;

  if (argc >= 2) {
    filename = argv[1];
  } else {
    filename = "output.pdf";
  }

  importRawConsole(list, range, rad);
  braggToPDF(list, range, rad, filename);

  return 0;
}
