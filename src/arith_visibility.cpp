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

#include "arith_visibility.h"

void Coprime::factorInteger(uint i, vector<uint>& factorization) {
  factorization.clear();

  for (ulong k = 2; k*k <= i; ++k) {
    while (i % k == 0) {
      factorization.push_back(k);
      i /= k;
    }
  }

  if (i > 1) factorization.push_back(i);
}

void Coprime::findTupleZ2(const int p, vec2i& out) {
  int x = ceil(sqrt(double(p)));
  int y;

  while (true) {
    y = lround(sqrt(double(x*x - p) / 2.0));

    if (2*y*y - x*x + p == 0) break;
    ++x;
  }

  out.set(x, y);
}

void Coprime::findTupleGI(const int p, vec2i& out) {
  int x = ceil(sqrt(double(p)));
  int y;

  while (true) {
    y = lround(sqrt(double(p - x*x)));

    if (x*x + y*y - p == 0) break;
    --x;
  }

  out.set(x, y);
}

void Coprime::findTupleES(const int p, vec2i& out) {
  int x = floor(sqrt(double(p * 4) / 3.0));
  int y;

  while (true) {
    const int t = lround(sqrt(-3*x*x + 4*p));

    if ((x + t) % 2 == 0) {
      y = (x + t) / 2;
      break;
    }

    if ((x - t) % 2 == 0) {
      y = (x - t) / 2;
      break;
    }

    --x;
  }

  out.set(x, y);
}

bool Coprime::pCond1Z2(const int p) {
  return (((p - 1) % 8 == 0) || (((p + 1) % 8 == 0)));
}

bool Coprime::pCond2Z2(const int p) {
  return (((p - 3) % 8 == 0) || (((p + 3) % 8 == 0)));
}

bool Coprime::pCond1GI(const int p) {
  return ((p - 1) % 4 == 0);
}

bool Coprime::pCond2GI(const int p) {
  return ((p - 3) % 4 == 0);
}

bool Coprime::pCond1ES(const int p) {
  return ((p - 1) % 3 == 0);
}

bool Coprime::pCond2ES(const int p) {
  return ((p - 2) % 3 == 0);
}

void Coprime::factorZ2(const vec2i& in, vector<vec2i>& factorization) {
  factorization.clear();

  const int norm = in.preNormZ2();
  if (abs(norm) == 1) return;

  vector<uint> primes;
  factorInteger(abs(norm), primes);

  for (vector<uint>::const_iterator k = primes.begin(); k != primes.end(); ++k) {
    if (*k % 2 == 0) {
      /* If we encounter the prime '2' in our factorization of the norm, then *
       * the factorization of the input contains the factor Sqrt[2].          */
      factorization.push_back(vec2i(0, 1));
      continue;
    }

    if (pCond2Z2(*k)) {
      factorization.push_back(vec2i(*k, 0));
      continue;
    }

    if (pCond1Z2(*k)) {
      vec2i t;
      Coprime::findTupleZ2(*k, t);

      if (in.isDivZ2(t)) factorization.push_back(t);
      if (in.isDivZ2(t.conjZ2())) factorization.push_back(t.conjZ2());
    }
  }
}

void Coprime::factorGI(const vec2i& in, vector<vec2i>& factorization) {
  factorization.clear();

  const int norm = in.preNormGI();
  if (abs(norm) == 1) return;

  vector<uint> primes;
  factorInteger(abs(norm), primes);

  for (vector<uint>::const_iterator k = primes.begin(); k != primes.end(); ++k) {
    if (*k % 2 == 0) {
      /* If we encounter the prime '2' in our factorization of the norm, then *
       * the factorization of the input contains the factor (1 + I).          */
      factorization.push_back(vec2i(1, 1));
      continue;
    }

    if (pCond2GI(*k)) {
      factorization.push_back(vec2i(*k, 0));
      continue;
    }

    if (pCond1GI(*k)) {
      vec2i t;
      Coprime::findTupleGI(*k, t);

      if (in.isDivGI(t)) factorization.push_back(t);
      if (in.isDivGI(t.conjGI())) factorization.push_back(t.conjGI());
    }
  }
}

void Coprime::factorES(const vec2i& in, vector<vec2i>& factorization) {
  factorization.clear();

  const int norm = in.preNormES();
  if (abs(norm) == 1) return;

  vector<uint> primes;
  factorInteger(abs(norm), primes);

  for (vector<uint>::const_iterator k = primes.begin(); k != primes.end(); ++k) {
    if (*k % 3 == 0) {
      /* If we encounter the prime '3' in our factorization of the norm, then *
       * the factorization of the input contains the factor (1 - omega).      */
      factorization.push_back(vec2i(1, -1));
      continue;
    }

    if (pCond2ES(*k)) {
      factorization.push_back(vec2i(*k, 0));
      continue;
    }

    if (pCond1ES(*k)) {
      vec2i t;
      Coprime::findTupleES(*k, t);

      if (in.isDivES(t)) factorization.push_back(t);
      if (in.isDivES(t.conjES())) factorization.push_back(t.conjES());
    }
  }
}

bool ArithVisibility::divTest2Free1Z2(const vec2i& in, const int p) {
  using namespace Coprime;

  vec2i t, a1, a2;

  Coprime::findTupleZ2(p, t);

  multZ2(in, t.squareZ2(), a1);
  multZ2(in, t.conjZ2().squareZ2(), a2);

  return (a1.isDiv(p*p) || a2.isDiv(p*p));
}

bool ArithVisibility::divTest2Free2Z2(const vec2i& in , const int p) {
  return in.isDiv(p*p);
}

bool ArithVisibility::visibility2FreeZ2(const vec2i& in) {
  using namespace Coprime;

  if (((in.x - in.y) % 2 == 0) && (in.x % 2 == 0))
    return false;

  const int norm = in.preNormZ2(); /* the algebraic norm */
  if (abs(norm) == 1) return true;

  vector<uint> primes;
  factorInteger(abs(norm), primes);

  for (vector<uint>::const_iterator k = primes.begin();
       k != primes.end(); ++k) {
    if (pCond1Z2(*k)) {
      if (divTest2Free1Z2(in, *k)) return false;
    } else {
      if (pCond2Z2(*k)) {
        if (divTest2Free2Z2(in, *k)) return false;
      }
    }
  }

  return true;
}

bool ArithVisibility::divTest2Free1GI(const vec2i& in, const int p) {
  using namespace Coprime;

  vec2i t, a1, a2;

  Coprime::findTupleGI(p, t);

  multGI(in, t.squareGI(), a1);
  multGI(in, t.conjGI().squareGI(), a2);

  return (a1.isDiv(p*p) || a2.isDiv(p*p));
}

bool ArithVisibility::divTest2Free2GI(const vec2i& in, const int p) {
  return in.isDiv(p*p);
}

bool ArithVisibility::visibility2FreeGI(const vec2i& in) {
  using namespace Coprime;

  if (((in.x - in.y) % 2 == 0) && (in.x % 2 == 0))
    return false;

  const int norm = in.preNormGI(); /* the algebraic norm */
  if (abs(norm) == 1) return true;

  vector<uint> primes;
  factorInteger(abs(norm), primes);

  for (vector<uint>::const_iterator k = primes.begin();
       k != primes.end(); ++k) {
    if (pCond1GI(*k)) {
      if (divTest2Free1GI(in, *k)) return false;
    } else {
      if (pCond2Z2(*k)) {
        if (divTest2Free2GI(in, *k)) return false;
      }
    }
  }

  return true;
}

vec2i ArithVisibility::denomZ2Fourier(const vec2i& in, const int in_c) {
  const vec2i x(in.y * 4, in.x * 2);
  const vec2i c(in_c, 0);
  const vec2i g(Coprime::gcdZ2(c, x));

  return c.divZ2(g);
}

double ArithVisibility::intensityZ2(const vec2i& denom) {
  static const double zetaZ2 = 48.0 * sqrt(2.0) / (Common::pi *
    Common::pi * Common::pi * Common::pi);

  if (denom.isZero()) return 0.0;

  vector<vec2i> primesZ2;
  Coprime::factorZ2(denom, primesZ2);

  double ret = 1.0;
  for (vector<vec2i>::const_iterator k = primesZ2.begin();
       k != primesZ2.end(); ++k) {
    const int pnorm = k->preNormZ2();
    ret *= (1.0 / (double(pnorm * pnorm) - 1.0));
  }

  ret *= (1.0 / (2.0 * sqrt(2.0) * zetaZ2));
  return ret;
}

vec2i ArithVisibility::denomGIFourier(const vec2i& in, const int in_c) {
  // This case is self-dual, so we don't need any transformation of the input.
  const vec2i c(in_c, 0);
  const vec2i g(Coprime::gcdGI(c, in));

  return c.divGI(g);
}

double ArithVisibility::intensityGI(const vec2i& denom) {
  static const double catalanC = 0.915965594177219015054603514932;

  if (denom.isZero()) return 0.0;

  vector<vec2i> primesGI;
  Coprime::factorGI(denom, primesGI);

  double ret = 1.0;
  for (vector<vec2i>::const_iterator k = primesGI.begin();
       k != primesGI.end(); ++k) {
    const int pnorm = k->preNormGI();
    ret *= (1.0 / (double(pnorm * pnorm) - 1.0));
  }

  ret *= (6.0 / (Common::pi * Common::pi * catalanC));
  return ret;
}

bool ArithVisibility::divTest3Free1Z2(const vec2i& in, const int p) {
  using namespace Coprime;

  vec2i t, a1, a2;

  Coprime::findTupleZ2(p, t);

  multZ2(in, t.cubeZ2(), a1);
  multZ2(in, t.conjZ2().cubeZ2(), a2);

  return (a1.isDiv(p*p*p) || a2.isDiv(p*p*p));
}

bool ArithVisibility::divTest3Free2Z2(const vec2i& in, const int p) {
  return in.isDiv(p*p*p);
}

bool ArithVisibility::visibility3FreeZ2(const vec2i& in) {
  using namespace Coprime;

  // Check if the input is divisible by (Sqrt[2])^3 = 2*Sqrt[2].
  if (in.isDivZ2(vec2i(0, 2))) return false;

  const int norm = in.preNormZ2(); /* the algebraic norm */
  if (abs(norm) == 1) return true;

  vector<uint> primes;
  factorInteger(abs(norm), primes);

  for (vector<uint>::const_iterator k = primes.begin();
       k != primes.end(); ++k) {
    if (pCond1Z2(*k)) {
      if (divTest3Free1Z2(in, *k)) return false;
    } else {
      if (pCond2Z2(*k)) {
        if (divTest3Free2Z2(in, *k)) return false;
      }
    }
  }

  return true;
}

bool ArithVisibility::divTest3Free1GI(const vec2i& in, const int p) {
  using namespace Coprime;

  vec2i t, a1, a2;

  Coprime::findTupleGI(p, t);

  multGI(in, t.cubeGI(), a1);
  multGI(in, t.conjGI().cubeGI(), a2);

  return (a1.isDiv(p*p*p) || a2.isDiv(p*p*p));
}

bool ArithVisibility::divTest3Free2GI(const vec2i& in, const int p) {
  return in.isDiv(p*p*p);
}

bool ArithVisibility::visibility3FreeGI(const vec2i& in) {
  using namespace Coprime;

  // Check if the input is divisible by (1-i)^3.
  if (in.isDivGI(vec2i(2, 2))) return false;

  const int norm = in.preNormGI(); /* the algebraic norm */
  if (abs(norm) == 1) return true;

  vector<uint> primes;
  factorInteger(abs(norm), primes);

  for (vector<uint>::const_iterator k = primes.begin();
       k != primes.end(); ++k) {
    if (pCond1GI(*k)) {
      if (divTest3Free1GI(in, *k)) return false;
    } else {
      if (pCond2GI(*k)) {
        if (divTest3Free2GI(in, *k)) return false;
      }
    }
  }

  return true;
}

void ArithVisibility::diffractionZ2(const vector<vec2iq>& in,
        vector<bragg>& out, clipfunc f) {
  out.clear();

  for (vector<vec2iq>::const_iterator k = in.begin(); k != in.end(); ++k) {
    const vec2d pos(k->minkowskiQ2());
    if (!f(pos)) continue;

    const vec2i denom(denomZ2Fourier(k->getNumerator(), k->getDenominator()));

    if (visibility3FreeZ2(denom)) {
      out.push_back(bragg(pos, intensityZ2(denom)));
    }
  }
}

bool ArithVisibility::clipFundamentalZ2(const vec2d& x) {
  using namespace Common;

  static const vec2d vzero(0.0, 0.0);
  static const vec2d vec1(0.5, 0.5);
  static const vec2d vec2(sqrt(2.0) / 4.0, -sqrt(2.0) / 4.0);
  static const double clipeps = 0.00001;

  if (clipeps + checkPosition(vec2 + vzero, vec2 + vec1, x) < 0) return false;
  if (clipeps + checkPosition(-vec2 + vec1, -vec2 + vzero, x) < 0) return false;
  if (clipeps + checkPosition(-vec1 + vzero, -vec1 + vec2, x) < 0) return false;
  if (clipeps + checkPosition(vec1 + vec2, vec1 + vzero, x) < 0) return false;

  return true;
}

void ArithVisibility::diffractionGI(const vector<vec2iq>& in,
              vector<bragg>& out, clipfunc f) {
  out.clear();

  for (vector<vec2iq>::const_iterator k = in.begin(); k != in.end(); ++k) {
    const vec2d pos(k->minkowskiGR());
    if (!f(pos)) continue;

    const vec2i denom(denomGIFourier(k->getNumerator(), k->getDenominator()));

    if (visibility3FreeGI(denom)) {
      out.push_back(bragg(pos, intensityGI(denom)));
    }
  }
}

bool ArithVisibility::clipFundamentalGI(const vec2d& x) {
  using namespace Common;

  static const vec2d vzero(0.0, 0.0);
  static const vec2d vec1(0.0, 1.0);
  static const vec2d vec2(1.0, 0.0);
  static const double clipeps = 0.00001;

  if (clipeps + checkPosition(vec2 + vzero, vec2 + vec1, x) < 0) return false;
  if (clipeps + checkPosition(-vec2 + vec1, -vec2 + vzero, x) < 0) return false;
  if (clipeps + checkPosition(-vec1 + vzero, -vec1 + vec2, x) < 0) return false;
  if (clipeps + checkPosition(vec1 + vec2, vec1 + vzero, x) < 0) return false;

  return true;
}

ostream& operator<<(ostream &os, const ArithVisibility::vec2iq& v) {
  os << '{' << v.getNumerator() << ',' << v.getDenominator() << '}';
  return os;
}

ostream& operator<<(ostream &os, const ArithVisibility::bragg& b) {
  os << "circle((" << b.getPosition().x << ',' << b.getPosition().y << "),"
     << b.getIntensity() << ",color='black',thickness=0.5)";
}

void vTableZ2(const uint r, Common::vec2ilist& table) {
  const int temp = int(ceil(double(r) / sqrt(2.0)));
  table.clear();
  table.reserve((r + 1) * (temp + 1));

  for (int i = -int(r); i <= int(r); ++i) {
    for (int j = -temp; j <= temp; ++j) {
      table.push_back(vec2i(i, j));
    }
  }
}

void vqTableRecipZ2(const uint r, const uint s,
                    vector<ArithVisibility::vec2iq>& table) {
  using namespace Coprime;
  using namespace ArithVisibility;

  table.clear();

  for (uint c = 1; c <= s; ++c) {
    for (int i = -int(r); i <= int(r); ++i) {
      for (int j = -int(r); j <= int(r); ++j) {
        const int g = gcdZ(uint(gcdZ(abs(2 * i), abs(j))), 4 * c);

        table.push_back(vec2iq(2 * i / g, j / g, 4 * c / g));
      }
    }
  }

  sort(table.begin(), table.end());
  table.erase(unique(table.begin(), table.end()), table.end());
}

void vTableGI(const uint r, Common::vec2ilist& table) {
  table.clear();
  table.reserve((r + 1) * (r + 1));

  for (int i = -int(r); i <= int(r); ++i) {
    for (int j = -int(r); j <= int(r); ++j) {
      table.push_back(vec2i(i, j));
    }
  }
}

void vqTableRecipGI(const uint r, const uint s,
                    vector<ArithVisibility::vec2iq>& table) {
  using namespace Coprime;
  using namespace ArithVisibility;

  table.clear();

  for (uint c = 1; c <= s; ++c) {
    for (int i = -int(r); i <= int(r); ++i) {
      for (int j = -int(r); j <= int(r); ++j) {
        const int g = gcdZ(uint(gcdZ(abs(i), abs(j))), c);

        table.push_back(vec2iq(i / g, j / g, c / g));
      }
    }
  }

  sort(table.begin(), table.end());
  table.erase(unique(table.begin(), table.end()), table.end());
}

void minmax(const vector<ArithVisibility::bragg>& input,
            vec2d& min, vec2d& max, double& radius){
  using namespace ArithVisibility;

  if (input.empty()) return;

  vector<bragg>::const_iterator k = input.begin();
  vec2d a(k->getPosition()), b(k->getPosition());
  double r = k->getIntensity();
  ++k;

  for (; k != input.end(); ++k) {
    const vec2d temp(k->getPosition());
    const double rtemp = k->getIntensity();

    if (temp.x < a.x) {
      a.x = temp.x;
    }

    if (temp.x > b.x) {
      b.x = temp.x;
    }

    if (temp.y < a.y) {
      a.y = temp.y;
    }

    if (temp.y > b.y) {
      b.y = temp.y;
    }

    if (rtemp > r) r = rtemp;
  }

  min = a; max = b;
  radius = r;
}

void toEPS(const vector<ArithVisibility::bragg>& input) {
  using namespace ArithVisibility;

  vec2d min, max;
  double radius;

  // Base width of the Postscript and offset to the borders.
  const int basewidth = 800;
  const int offset = 50;

  if (input.empty()) return;

  minmax(input, min, max, radius);

  const double xfactor = std::max(abs(min.x), abs(max.x));
  const double yfactor = std::max(abs(min.y), abs(max.y));

  const int height = basewidth * lround(yfactor / xfactor);

  // Write header with bounding box information
  cout << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
  cout << "%%BoundingBox: 0 0 " << basewidth << ' ' << height << endl;
  cout << "%%BeginProlog" << endl;
  cout << "%%EndProlog" << endl;

  // Translate origin to the middle of the page
  cout << (basewidth / 2) << ' ' << (height / 2) << " translate" << endl;

  // Apply scaling
  cout << lround(double(basewidth - offset) / (2.0 * (xfactor + radius))) << ' '
       << lround(double(height - offset) / (2.0 * (yfactor + radius)))
       << " scale" << endl;

  cout << "0.002 setlinewidth" << endl; // TODO: base width on min/max

  for (vector<bragg>::const_iterator k = input.begin(); k != input.end(); ++k) {
    cout << "newpath" << endl;
    cout << k->getPosition().x << ' ' << k->getPosition().y << ' '
         << k->getIntensity() << " 0 360 arc" << endl;
    cout << "stroke" << endl;
  }

  cout << "%%EOF" << endl;
}

double linscale(double x) {
  return x * 0.4;
}

int main(int argc, char* argv[]) {
  using namespace ArithVisibility;

  vector<vec2iq> large_table;
  vector<bragg> diffraction;

  vqTableRecipZ2(35, 29, large_table);
  diffractionZ2(large_table, diffraction, clipFundamentalZ2);

  for (vector<bragg>::iterator k = diffraction.begin(); k != diffraction.end(); ++k) {
    k->apply(linscale);
  }

  toEPS(diffraction);

  /*vector<vec2iq> large_table;
  vector<bragg> diffraction;

  vqTableRecipGI(50, 47, large_table);
  diffractionGI(large_table, diffraction, clipFundamentalGI);

  for (vector<bragg>::iterator k = diffraction.begin(); k != diffraction.end(); ++k) {
    k->apply(linscale);
  }

  toEPS(diffraction);&/

  /*Common::vec2ilist large_table;
  Common::vec2ilist sqfree_table;

  for (uint size = 50; size < 2000; size += 50) {
    vTableGI(size, large_table);
  
    for (Common::vec2ilist::const_iterator i = large_table.begin();
         i != large_table.end(); ++i) {
      if (ArithVisibility::visibility2FreeGI(*i)) sqfree_table.push_back(*i);
    }
  
    cout << size << ": ";
    cout << (double(sqfree_table.size()) / double(large_table.size())) << endl;
    sqfree_table.clear();
  }*/

  return 0;
}
