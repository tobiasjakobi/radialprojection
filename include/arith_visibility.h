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

#ifndef _ARITH_VISIBILITY_H_
#define _ARITH_VISIBILITY_H_

#include "common.h"
#include "visibility.h"

namespace Coprime {

  // Very primitive integer factorization (should work for small numbers)
  void factorInteger(uint i, vector<uint>& factorization);

  /* Assuming that p = +1 or -1 (mod 8), this finds the tuple *
   * (m,n) that solves the equation algnorm(m,n) = p.         */
  void findTupleZ2(const int p, vec2i& out);

  /* Assuming that p = +1 (mod 4), this finds the tuple (m,n) *
   * that solves the equation algnorm(m,n) = p.               */
  void findTupleGI(const int p, vec2i& out);

  /* Assuming that p = +1 (mod 3), this finds the tuple (m,n) *
   * that solves the equation algnorm(m,n) = p.               */
  void findTupleES(const int p, vec2i& out);

  /* Assuming that p = +1 or -1 (mod 5), this finds the tuple (m,n) *
   * that solves the equation algnorm(m,n) = p.                     */
  void findTupleGM(const int p, vec2i& out);

  /* The two conditions that apply to our visibility checks (for Z[Sqrt[2]]): *
   * p = +1 or -1 (mod 8) (first condition)                                   *
   * p = +3 or -3 (mod 8) (second condition)                                  */
  bool pCond1Z2(const int p);
  bool pCond2Z2(const int p);

  /* The two conditions (for the Gaussian Integers):  *
   * p = +1 (mod 4) (first condition)                 *
   * p = +3 (mod 4) (second condition)                */
  bool pCond1GI(const int p);
  bool pCond2GI(const int p);

  /* The two conditions (for the Eisenstein Integers): *
   * p = +1 (mod 3) (first condition)                  *
   * p = +2 (mod 3) (second condition)                 */
  bool pCond1ES(const int p);
  bool pCond2ES(const int p);

  /* The two conditions (for the Z[tau] case): *
   * p = +1 or -1 (mod 5) (first condition)    *
   * p = +2 or -2 (mod 5) (second condition)   */
  bool pCond1GM(const int p);
  bool pCond2GM(const int p);

  /* Do prime factorization in Z[Sqrt[2]]. This takes 'in' as an *
   * element of Z[Sqrt[2]] and returns a list of primes that     *
   * divide 'in' (without multiplicity).                         */
  void factorZ2(const vec2i& in, vector<vec2i>& factorization);

  // Do prime factorization in the Gaussian integers.
  void factorGI(const vec2i& in, vector<vec2i>& factorization);

  // Do prime factorization in the Eisenstein integers.
  void factorES(const vec2i& in, vector<vec2i>& factorization);

  // Do prime factorization in Z[tau] (tau being the golden mean).
  void factorGM(const vec2i& in, vector<vec2i>& factorization);
};

namespace ArithVisibility {
  /* Division test (square-free case) for the primes which *
   * satisfy the first condition.                          */
  bool divTest2Free1Z2(const vec2i& in, const int p);

  /* Division test (square-free case) for the primes which *
   * satisfy the second condition.                         */
  bool divTest2Free2Z2(const vec2i& in, const int p);

  // Check for square-free 'visibility' of an element of Z[Sqrt[2]].
  bool visibility2FreeZ2(const vec2i& in);

  // Check for square-free 'visibility' of an element of the Gaussian Integers.
  bool divTest2Free1GI(const vec2i& in, const int p);
  bool divTest2Free2GI(const vec2i& in, const int p);
  bool visibility2FreeGI(const vec2i& in);

  // Check for square-free 'visibility' of an element of the Eisenstein Integers.
  bool divTest2Free1ES(const vec2i& in, const int p);
  bool divTest2Free2ES(const vec2i& in, const int p);
  bool visibility2FreeES(const vec2i& in);

  // Check for square-free 'visibility' of an element of Z[tau].
  bool divTest2Free1GM(const vec2i& in, const int p);
  bool divTest2Free2GM(const vec2i& in, const int p);
  bool visibility2FreeGM(const vec2i& in);

  /* Let x = in / c, an element of Q(Sqrt[2]), then this computes    *
   * the denominator in the Fourier module Z[Sqrt[2]] * (Sqrt[2]/4). */
  vec2i denomZ2Fourier(const vec2i& in, const int in_c);

  /* Compute the intensity of 'denom', an element of Z[Sqrt[2]].       *
   * The input are going to be denominators of elements of Q(Sqrt[2]). */
  double intensityZ2(const vec2i& denom);

  // Computation of denominator and intensity for the Gaussian Integers.
  vec2i denomGIFourier(const vec2i& in, const int in_c);
  double intensityGI(const vec2i& denom);

  // Computation of denominator and intensity for the Eisenstein Integers.
  vec2i denomESFourier(const vec2i& in, const int in_c);
  double intensityES(const vec2i& denom);

  // Computation of denominator and intensity for Z[tau].
  vec2i denomGMFourier(const vec2i& in, const int in_c);
  double intensityGM(const vec2i& denom);

  // Division tests for the primes in the cube-free case.
  bool divTest3Free1Z2(const vec2i& in, const int p);
  bool divTest3Free2Z2(const vec2i& in, const int p);

  // Check for cube-free 'visibility' of an element of Z[Sqrt[2]].
  bool visibility3FreeZ2(const vec2i& in);

  // Check for cube-free 'visibility' of an element of the Gaussian Integers.
  bool divTest3Free1GI(const vec2i& in, const int p);
  bool divTest3Free2GI(const vec2i& in, const int p);
  bool visibility3FreeGI(const vec2i& in);

  // Check for cube-free 'visibility' of an element of the Eisenstein Integers.
  bool divTest3Free1ES(const vec2i& in, const int p);
  bool divTest3Free2ES(const vec2i& in, const int p);
  bool visibility3FreeES(const vec2i& in);

  // Check for cube-free 'visibility' of an element of Z[tau].
  bool divTest3Free1GM(const vec2i& in, const int p);
  bool divTest3Free2GM(const vec2i& in, const int p);
  bool visibility3FreeGM(const vec2i& in);

  class vec2iq {
  private:
    vec2i numerator;
    int denominator;

  public:
    vec2iq() {}
    vec2iq(int a, int b, int c) : numerator(a, b), denominator(c) {}

    bool operator==(const vec2iq& v) const {
      return ((numerator == v.numerator) && (denominator == v.denominator));
    }

    /* Lexicographic ordering, this is needed to use sorting *
     * algorithms of STL containers.                         */
    bool operator<(const vec2iq& v) const {
      if (numerator < v.numerator) {
        return true;
      }

      if (numerator == v.numerator) {
        return (denominator < v.denominator);
      }

      return false;
    }

    // Minkowski embedding for Q(Sqrt[2])
    vec2d minkowskiQ2() const {
      const double inv = 1.0 / double(denominator);

      return (numerator.minkowskiZ2() * inv);
    }

    // Minkowski embedding for the Gaussian rationals (by this we mean Q(I))
    vec2d minkowskiGR() const {
      const double inv = 1.0 / double(denominator);

      return (vec2d(double(numerator.x), double(numerator.y)) * inv);
    }

    // Minkowski embedding for the Eisenstein rationals
    vec2d minkowskiER() const {
      const double inv = 1.0 / double(denominator);

      return (numerator.minkowskiES() * inv);
    }

    // Minkowski embeeding for Q(tau) (tau being the golden mean)
    // GMR = Golden Mean Rationals
    vec2d minkowskiGMR() const {
      const double inv = 1.0 / double(denominator);

      return (numerator.minkowskiGM() * inv);
    }

    const vec2i& getNumerator() const {
      return numerator;
    }

    int getDenominator() const {
      return denominator;
    }
  };

  typedef double (&scalefunc)(double);
  typedef bool (&clipfunc)(const vec2d&);

  class bragg {
  private:
    vec2d position;
    double intensity;

  public:
    bragg() {}
    bragg(double x, double y, double i) : position(x, y), intensity(i) {}
    bragg(const vec2d& p, double i) : position(p), intensity(i) {}

    void apply(scalefunc f) {
      intensity = f(intensity);
    }

    const vec2d& getPosition() const {
      return position;
    }

    double getIntensity() const {
      return intensity;
    }

    void rotate(double degree);
  };

  const double clipeps = 0.00001;

  // Compute the diffraction of the square-free points of Z[Sqrt[2]].
  void diffractionZ2(const vector<vec2iq>& in, vector<bragg>& out,
                     clipfunc f);

  /* Box clipping test that checks if the point 'x' is contained in an *
   * arrangement of fundamental domains of the reciprocal lattice.     *
   * Currently two times two domains (around the origin) are cut.      */
  bool clipFundamentalZ2(const vec2d& x);

  // Compute the diffraction of the square-free points of the Gaussian Integers.
  void diffractionGI(const vector<vec2iq>& in, vector<bragg>& out,
                     clipfunc f);
  bool clipFundamentalGI(const vec2d& x);

  // Compute the diffraction of the square-free points of the Eisenstein Integers.
  void diffractionES(const vector<vec2iq>& in, vector<bragg>& out,
                     clipfunc f);
  bool clipFundamentalES(const vec2d& x);

  // Compute the diffraction of the square-free points of Z[tau].
  void diffractionGM(const vector<vec2iq>& in, vector<bragg>& out,
                     clipfunc f);
  bool clipFundamentalGM(const vec2d& x);

  struct VisOpZ2 {
    typedef Common::vec2ilist list_type;
    static const double epsilon;

    static double angle(const vec2i& a);
    static vec2d toR2(const vec2i& a);
    static bool rayTest(const vec2i& a, const vec2i& b);
  };

  struct VisOpGI {
    typedef Common::vec2ilist list_type;
    static const double epsilon;

    static double angle(const vec2i& a);
    static vec2d toR2(const vec2i& a);
    static bool rayTest(const vec2i& a, const vec2i& b);
  };

  struct VisOpES {
    typedef Common::vec2ilist list_type;
    static const double epsilon;

    static double angle(const vec2i& a);
    static vec2d toR2(const vec2i& a);
    static bool rayTest(const vec2i& a, const vec2i& b);
  };

  struct VisOpGM {
    typedef Common::vec2ilist list_type;
    static const double epsilon;

    static double angle(const vec2i& a);
    static vec2d toR2(const vec2i& a);
    static bool rayTest(const vec2i& a, const vec2i& b);
  };

  typedef VisTest::VisibleList<VisOpZ2> VisListZ2;
  typedef VisTest::VisibleList<VisOpGI> VisListGI;
  typedef VisTest::VisibleList<VisOpES> VisListES;
  typedef VisTest::VisibleList<VisOpGM> VisListGM;

  void visCircleZ2(const uint r, Common::vec2ilist& out, bool radialproj);
  void visCircleZ2Fast(const uint r, Common::vec2ilist& out, bool radialproj);

  void visCircleGI(const uint r, Common::vec2ilist& out, bool radialproj);
  void visCircleGIFast(const uint r, Common::vec2ilist& out, bool radialproj);

  void visCircleES(const uint r, Common::vec2ilist& out, bool radialproj);
  void visCircleESFast(const uint r, Common::vec2ilist& out, bool radialproj);

  void visCircleGM(const uint r, Common::vec2ilist& out, bool radialproj);
  void visCircleGMFast(const uint r, Common::vec2ilist& out, bool radialproj);

  void radialProjZ2(const uint r, Common::dlist& out);
  void radialProjGI(const uint r, Common::dlist& out);
  void radialProjES(const uint r, Common::dlist& out);
  void radialProjGM(const uint r, Common::dlist& out);

  /*
   * Class defining a occupation array belonging to a rectangular
   * centered (around zero) subset of the Gaussian integers.
   * We use this class to produce bitmaps in plain PBM format from
   * patches of square-free Gaussian integers.
   *
   * Re-uses code from class vismap (defined in the chiral_radial
   * header, here used for computing geometric visibility for the
   * planar chair tiling).
   */
  class OccupationArray {
  private:
    typedef vector<bool> boolvec;

    int xrange[2];
    int yrange[2];
    int rowlen;

    boolvec oarray;

    bool getbit(int x, int y) const {
      if (xrange[0] <= x && x <= xrange[1] &&
          yrange[0] <= y && y <= yrange[1]) {
        const int pos = (y - yrange[0]) * rowlen + (x - xrange[0]);
        assert(pos >= 0);

        return oarray[pos];
      } else {
        // Return 'false' for all positions outside range
        return false;
      }
    }

    void setbit(int x, int y) {
      if (xrange[0] <= x && x <= xrange[1] &&
          yrange[0] <= y && y <= yrange[1]) {
        const int pos = (y - yrange[0]) * rowlen + (x - xrange[0]);
        assert(pos >= 0);

        oarray[pos] = true;
      }
    }

  public:
    // constructor
    OccupationArray(const Common::vec2ilist& input) {
      xrange[0] = -1;
      xrange[1] = 0;
      yrange[0] = -1;
      yrange[1] = 0;

      // extract ranges from input data
      for (Common::vec2ilist::const_iterator i = input.begin();
           i != input.end(); ++i) {
        int x = i->x;
        int y = i->y;

        if (x < xrange[0])
          xrange[0] = x;
        if (x > xrange[1])
          xrange[1] = x;
        if (y < yrange[0])
          yrange[0] = y;
        if (y > yrange[1])
          yrange[1] = y;
      }

      rowlen = xrange[1] - xrange[0] + 1;
      oarray.resize((yrange[1] - yrange[0] + 1) * rowlen, false);

      for (Common::vec2ilist::const_iterator i = input.begin();
           i != input.end(); ++i)
        setbit(i->x, i->y);
    }

    // destructor
    ~OccupationArray() {}

    // copy-constructor
    OccupationArray(const OccupationArray& oa) : rowlen(oa.rowlen),
                                                 oarray(oa.oarray) {
      xrange[0] = oa.xrange[0];
      xrange[1] = oa.xrange[1];
      yrange[0] = oa.yrange[0];
      yrange[1] = oa.yrange[1];
    }

    void writePlainPBM() const {
      // write header
      cout << "P1" << endl;
      cout << rowlen << ' ' << (yrange[1] - yrange[0] + 1) << endl;

      int rowcnt = 0;
      for (boolvec::const_iterator i = oarray.begin(); i != oarray.end(); ++i) {
        cout << ((*i) ? '1' : '0') << ' ';

        rowcnt++;

        if (rowcnt == rowlen) {
          rowcnt = 0;
          cout << endl;
        }
      }
    }
  };

  /*
   * Test for determining the subgroup that leaves the square-free
   * Gaussian integers invariant.
   */
  bool powerTest(const mtx2x2i& m, uint n);
  bool multTest(const mtx2x2i& m, const vec2i& v);

};

ostream& operator<<(ostream &os, const ArithVisibility::vec2iq& v);

// Format output so that we can use it in SAGE
ostream& operator<<(ostream &os, const ArithVisibility::bragg& b);

void vTableZ2(const uint r, Common::vec2ilist& table);
void vqTableRecipZ2(const uint r, const uint s,
                    vector<ArithVisibility::vec2iq>& table);

void vTableGI(const uint r, Common::vec2ilist& table);
void vqTableRecipGI(const uint r, const uint s,
                    vector<ArithVisibility::vec2iq>& table);

void vTableES(const uint r, Common::vec2ilist& table);
void vqTableRecipES(const uint r, const uint s,
                    vector<ArithVisibility::vec2iq>& table);

void vTableGM(const uint r, Common::vec2ilist& table);
void vqTableRecipGM(const uint r, const uint s,
                    vector<ArithVisibility::vec2iq>& table);

void vCircleZ2(const uint r, Common::vec2ilist& table);
void vCircleGI(const uint r, Common::vec2ilist& table);
void vCircleES(const uint r, Common::vec2ilist& table);
void vCircleGM(const uint r, Common::vec2ilist& table);

void minmax(const vector<ArithVisibility::bragg>& input, vec2d& min,
            vec2d& max, double& radius);

// Format output into encapsulated Postscript
void toEPS(const vector<ArithVisibility::bragg>& input, bool fill);

// Different method of writing an encapsulated Postscript.
// Only filled circles are supported.
void toEPS2(const vector<ArithVisibility::bragg>& input);

// Format output into Scalable Vector Graphics (SVG)
void toSVG(const vector<ArithVisibility::bragg>& input);

/* Export the diffraction pattern (which is the list of Bragg peaks) in *
 * a raw format to the console. The data can then be used with the      *
 * PDF writer part of the program to plot the diffraction to a PDF.     */
void exportRawConsole(const vector<ArithVisibility::bragg>& input);

#endif
