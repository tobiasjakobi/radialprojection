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

int main(int argc, char* argv[]) {
  cairo_surface_t *surface;
  cairo_t *cr;
  const char *filename = "testimage.pdf";

  surface = cairo_pdf_surface_create(filename, 80.0, 80.0);
  cr = cairo_create(surface);

  cairo_set_source_rgba(cr, 1.0, 0.0, 0.0, 1.0);
  cairo_set_line_width(cr, 1.0);
  cairo_set_fill_rule(cr, CAIRO_FILL_RULE_WINDING);

  cairo_arc(cr, 20.0, 20.0, 10.0, 0.0, 2*M_PI);
  cairo_fill(cr);
  cairo_stroke(cr);

  cairo_set_source_rgba(cr, 0.0, 0.0, 1.0, 1.0);
  cairo_arc(cr, 30.0, 20.0, 10.0, 0.0, 2*M_PI);
  cairo_fill(cr);
  cairo_stroke(cr);

  cairo_set_source_rgba(cr, 0.0, 1.0, 0.0, 1.0);
  cairo_arc(cr, 25.0, 28.6, 10.0, 0.0, 2*M_PI);
  cairo_fill(cr);
  cairo_stroke(cr);

  cairo_destroy(cr);
  cairo_surface_destroy(surface);

  return 0;
}
