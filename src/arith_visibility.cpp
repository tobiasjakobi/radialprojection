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

#include "arith_visibility.h"

#include <sstream>

const double ArithVisibility::VisOpZ2::epsilon =
  2.0 * numeric_limits<double>::epsilon();

const double ArithVisibility::VisOpGI::epsilon =
  2.0 * numeric_limits<double>::epsilon();

const double ArithVisibility::VisOpES::epsilon =
  2.0 * numeric_limits<double>::epsilon();

const double ArithVisibility::VisOpGM::epsilon =
  2.0 * numeric_limits<double>::epsilon();

void Coprime::factorInteger(uint i, vector<uint>& factorization) {
  factorization.clear();

  for (ulong k = 2; k*k <= i; ++k) {
    if (i % k == 0) {
      factorization.push_back(k);

      while (true) {
        i /= k;
        if (i % k != 0) break;
      }
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
  assert(out.preNormZ2() == p);
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
  assert(out.preNormGI() == p);
}

void Coprime::findTupleES(const int p, vec2i& out) {
  int x = floor(sqrt(double(p * 4) / 3.0));
  int y;

  while (true) {
    const int t = lround(sqrt(-3*x*x + 4*p));

    if (t*t + 3*x*x - 4*p == 0) {
      if ((x + t) % 2 == 0) {
        y = (x + t) / 2;
        break;
      }

      if ((x - t) % 2 == 0) {
        y = (x - t) / 2;
        break;
      }
    }

    --x;
  }

  out.set(x, y);
  assert(out.preNormES() == p);
}

void Coprime::findTupleGM(const int p, vec2i& out) {
  int x;
  int y = 0;

  while (true) {
    const int t = lround(sqrt(5*y*y + 4*p));

    if (t*t - 5*y*y - 4*p == 0) {
      if ((-y + t) % 2 == 0) {
        x = (-y + t) / 2;
        break;
      }

      if ((-y - t) % 2 == 0) {
        x = (-y - t) / 2;
        break;
      }
    }

    ++y;
  }

  out.set(x, y);
  assert(out.preNormZTau() == p);
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

bool Coprime::pCond1GM(const int p) {
  return (((p - 1) % 5 == 0) || (((p + 1) % 5 == 0)));
}

bool Coprime::pCond2GM(const int p) {
  return (((p - 2) % 5 == 0) || (((p + 2) % 5 == 0)));
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

void Coprime::factorGM(const vec2i& in, vector<vec2i>& factorization) {
  factorization.clear();

  const int norm = in.preNormZTau();
  if (abs(norm) == 1) return;

  vector<uint> primes;
  factorInteger(abs(norm), primes);

  for (vector<uint>::const_iterator k = primes.begin(); k != primes.end(); ++k) {
    if (*k % 5 == 0) {
      /* If we encounter the prime '5' in our factorization of the norm, then *
       * the factorization of the input contains the factor (-1 + 2*tau).     */
      factorization.push_back(vec2i(-1, 2));
      continue;
    }

    if (pCond2GM(*k)) {
      factorization.push_back(vec2i(*k, 0));
      continue;
    }

    if (pCond1GM(*k)) {
      vec2i t;
      Coprime::findTupleGM(*k, t);

      if (in.isDivGM(t)) factorization.push_back(t);
      if (in.isDivGM(t.conjGM())) factorization.push_back(t.conjGM());
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
      if (pCond2GI(*k)) {
        if (divTest2Free2GI(in, *k)) return false;
      }
    }
  }

  return true;
}

bool ArithVisibility::divTest2Free1ES(const vec2i& in, const int p) {
  using namespace Coprime;

  vec2i t, a1, a2;

  Coprime::findTupleES(p, t);

  multES(in, t.squareES(), a1);
  multES(in, t.conjES().squareES(), a2);

  return (a1.isDiv(p*p) || a2.isDiv(p*p));
}

bool ArithVisibility::divTest2Free2ES(const vec2i& in, const int p) {
  return in.isDiv(p*p);
}

bool ArithVisibility::visibility2FreeES(const vec2i& in) {
  using namespace Coprime;

  if (((in.x - in.y) % 3 == 0) && (in.x % 3 == 0))
    return false;

  const int norm = in.preNormES(); /* the algebraic norm */
  if (abs(norm) == 1) return true;

  vector<uint> primes;
  factorInteger(abs(norm), primes);

  for (vector<uint>::const_iterator k = primes.begin();
       k != primes.end(); ++k) {
    if (pCond1ES(*k)) {
      if (divTest2Free1ES(in, *k)) return false;
    } else {
      if (pCond2ES(*k)) {
        if (divTest2Free2ES(in, *k)) return false;
      }
    }
  }

  return true;
}

bool ArithVisibility::divTest2Free1GM(const vec2i& in, const int p) {
  using namespace Coprime;

  vec2i t, a1, a2;

  Coprime::findTupleGM(p, t);

  multZTau(in, t.squareGM(), a1);
  multZTau(in, t.conjGM().squareGM(), a2);

  return (a1.isDiv(p*p) || a2.isDiv(p*p));
}

bool ArithVisibility::divTest2Free2GM(const vec2i& in, const int p) {
  return in.isDiv(p*p);
}

bool ArithVisibility::visibility2FreeGM(const vec2i& in) {
  using namespace Coprime;

  if (((in.x - in.y) % 5 == 0) && (in.x % 5 == 0))
    return false;

  const int norm = in.preNormZTau(); /* the algebraic norm */
  if (abs(norm) == 1) return true;

  vector<uint> primes;
  factorInteger(abs(norm), primes);

  for (vector<uint>::const_iterator k = primes.begin();
       k != primes.end(); ++k) {
    if (pCond1GM(*k)) {
      if (divTest2Free1GM(in, *k)) return false;
    } else {
      if (pCond2GM(*k)) {
        if (divTest2Free2GM(in, *k)) return false;
      }
    }
  }

  return true;
}

vec2i ArithVisibility::denomZ2Fourier(const vec2i& in, const int in_c) {
  const vec2i x(in.x * 2, in.y * 4);
  const vec2i c(in_c, 0);
  const vec2i g(Coprime::gcdZ2(c, x));

  return c.divZ2(g);
}

double ArithVisibility::intensityZ2(const vec2i& denom) {
  if (denom.isZero()) return 0.0;

  vector<vec2i> primesZ2;
  Coprime::factorZ2(denom, primesZ2);

  double ret = 1.0;
  for (vector<vec2i>::const_iterator k = primesZ2.begin();
       k != primesZ2.end(); ++k) {
    const int pnorm = k->preNormZ2();
    ret *= (1.0 / (double(pnorm * pnorm) - 1.0));
  }

  ret *= (1.0 / (2.0 * sqrt(2.0) * Constants::zetaZ2));
  return (ret*ret);
}

vec2i ArithVisibility::denomGIFourier(const vec2i& in, const int in_c) {
  // This case is self-dual, so we don't need any transformation of the input.
  const vec2i c(in_c, 0);
  const vec2i g(Coprime::gcdGI(c, in));

  return c.divGI(g);
}

double ArithVisibility::intensityGI(const vec2i& denom) {
  if (denom.isZero()) return 0.0;

  vector<vec2i> primesGI;
  Coprime::factorGI(denom, primesGI);

  double ret = 1.0;
  for (vector<vec2i>::const_iterator k = primesGI.begin();
       k != primesGI.end(); ++k) {
    const int pnorm = k->preNormGI();
    ret *= (1.0 / (double(pnorm * pnorm) - 1.0));
  }

  ret *= (6.0 / (Constants::pi * Constants::pi * Constants::catalanC));
  return (ret*ret);
}

vec2i ArithVisibility::denomESFourier(const vec2i& in, const int in_c) {
  const vec2i x(in.x * 2 - in.y, -in.x + in.y * 2);
  const vec2i d(in_c * 2, 0);
  const vec2i g(Coprime::gcdES(d, x));

  return d.divES(g);
}

double ArithVisibility::intensityES(const vec2i& denom) {
  if (denom.isZero()) return 0.0;

  vector<vec2i> primesES;
  Coprime::factorES(denom, primesES);

  double ret = 1.0;
  for (vector<vec2i>::const_iterator k = primesES.begin();
       k != primesES.end(); ++k) {
    const int pnorm = k->preNormES();
    ret *= (1.0 / (double(pnorm * pnorm) - 1.0));
  }

  ret *= (2.0 / (sqrt(3.0) * Constants::zetaES));
  return (ret*ret);
}

vec2i ArithVisibility::denomGMFourier(const vec2i& in, const int in_c) {
  const vec2i x(in.x * 2 + in.y, in.x + in.y * 3);
  const vec2i c(in_c, 0);
  const vec2i g(Coprime::gcdZTau(c, x));

  return c.divGM(g);
}

double ArithVisibility::intensityGM(const vec2i& denom) {
  if (denom.isZero()) return 0.0;

  vector<vec2i> primesGM;
  Coprime::factorGM(denom, primesGM);

  double ret = 1.0;
  for (vector<vec2i>::const_iterator k = primesGM.begin();
       k != primesGM.end(); ++k) {
    const int pnorm = k->preNormZTau();
    ret *= (1.0 / (double(pnorm * pnorm) - 1.0));
  }

  ret *= (1.0 / (sqrt(5.0) * Constants::zetaGM));
  return (ret*ret);
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

bool ArithVisibility::divTest3Free1ES(const vec2i& in, const int p) {
  using namespace Coprime;

  vec2i t, a1, a2;

  Coprime::findTupleES(p, t);

  multES(in, t.cubeES(), a1);
  multES(in, t.conjES().cubeES(), a2);

  return (a1.isDiv(p*p*p) || a2.isDiv(p*p*p));
}

bool ArithVisibility::divTest3Free2ES(const vec2i& in, const int p) {
  return in.isDiv(p*p*p);
}

bool ArithVisibility::visibility3FreeES(const vec2i& in) {
  using namespace Coprime;

  // Check if the input is divisible by (1-omega)^3.
  if (in.isDivES(vec2i(3, 6))) return false;

  const int norm = in.preNormES(); /* the algebraic norm */
  if (abs(norm) == 1) return true;

  vector<uint> primes;
  factorInteger(abs(norm), primes);

  for (vector<uint>::const_iterator k = primes.begin();
       k != primes.end(); ++k) {
    if (pCond1ES(*k)) {
      if (divTest3Free1ES(in, *k)) return false;
    } else {
      if (pCond2ES(*k)) {
        if (divTest3Free2ES(in, *k)) return false;
      }
    }
  }

  return true;
}

bool ArithVisibility::divTest3Free1GM(const vec2i& in, const int p) {
  using namespace Coprime;

  vec2i t, a1, a2;

  Coprime::findTupleGM(p, t);

  multZTau(in, t.cubeGM(), a1);
  multZTau(in, t.conjGM().cubeGM(), a2);

  return (a1.isDiv(p*p*p) || a2.isDiv(p*p*p));
}

bool ArithVisibility::divTest3Free2GM(const vec2i& in, const int p) {
  return in.isDiv(p*p*p);
}

bool ArithVisibility::visibility3FreeGM(const vec2i& in) {
  using namespace Coprime;

  // Check if the input is divisible by (-1+2*tau)^3.
  if (in.isDivGM(vec2i(-5, 10))) return false;

  const int norm = in.preNormZTau(); /* the algebraic norm */
  if (abs(norm) == 1) return true;

  vector<uint> primes;
  factorInteger(abs(norm), primes);

  for (vector<uint>::const_iterator k = primes.begin();
       k != primes.end(); ++k) {
    if (pCond1GM(*k)) {
      if (divTest3Free1GM(in, *k)) return false;
    } else {
      if (pCond2GM(*k)) {
        if (divTest3Free2GM(in, *k)) return false;
      }
    }
  }

  return true;
}

void ArithVisibility::bragg::rotate(double degree) {
  const double rad = 2.0 * Constants::pi * (degree / 360.0);
  const double a = sin(rad);
  const double b = cos(rad);

  const vec2d temp(position.applyRotation(a, b));
  position.set(temp[0], temp[1]);
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

  /*if (clipeps + checkPosition(vec2 + vzero, vec2 + vec1, x) < 0) return false;
  if (clipeps + checkPosition(-vec2 + vec1, -vec2 + vzero, x) < 0) return false;
  if (clipeps + checkPosition(-vec1 + vzero, -vec1 + vec2, x) < 0) return false;
  if (clipeps + checkPosition(vec1 + vec2, vec1 + vzero, x) < 0) return false;*/
  if (clipeps + checkPosition(vzero, vec1, x) < 0) return false;
  if (clipeps + checkPosition(-vec2 + vec1, -vec2 + vzero, x) < 0) return false;
  if (clipeps + checkPosition(vzero, vec2, x) < 0) return false;
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

  if (clipeps + checkPosition(vec2 + vzero, vec2 + vec1, x) < 0) return false;
  if (clipeps + checkPosition(-vec2 + vec1, -vec2 + vzero, x) < 0) return false;
  if (clipeps + checkPosition(-vec1 + vzero, -vec1 + vec2, x) < 0) return false;
  if (clipeps + checkPosition(vec1 + vec2, vec1 + vzero, x) < 0) return false;

  return true;
}

void ArithVisibility::diffractionES(const vector<vec2iq>& in,
              vector<bragg>& out, clipfunc f) {
  out.clear();

  for (vector<vec2iq>::const_iterator k = in.begin(); k != in.end(); ++k) {
    const vec2d pos(k->minkowskiER());
    if (!f(pos)) continue;

    const vec2i denom(denomESFourier(k->getNumerator(), k->getDenominator()));

    if (visibility3FreeES(denom)) {
      out.push_back(bragg(pos, intensityES(denom)));
    }
  }
}

bool ArithVisibility::clipFundamentalES(const vec2d& x) {
  using namespace Common;

  static const vec2d vzero(0.0, 0.0);
  static const vec2d vec1(0.0, 2.0 / sqrt(3.0));
  static const vec2d vec2(1.0, 1.0 / sqrt(3.0));

  if (clipeps + checkPosition(vec2 + vzero, vec2 + vec1, x) < 0) return false;
  if (clipeps + checkPosition(-vec2 + vec1, -vec2 + vzero, x) < 0) return false;
  if (clipeps + checkPosition(-vec1 + vzero, -vec1 + vec2, x) < 0) return false;
  if (clipeps + checkPosition(vec1 + vec2, vec1 + vzero, x) < 0) return false;

  return true;
}

void ArithVisibility::diffractionGM(const vector<vec2iq>& in,
              vector<bragg>& out, clipfunc f) {
  out.clear();

  for (vector<vec2iq>::const_iterator k = in.begin(); k != in.end(); ++k) {
    const vec2d pos(k->minkowskiGMR());
    if (!f(pos)) continue;

    const vec2i denom(denomGMFourier(k->getNumerator(), k->getDenominator()));

    if (visibility3FreeGM(denom)) {
      out.push_back(bragg(pos, intensityGM(denom)));
    }
  }
}

bool ArithVisibility::clipFundamentalGM(const vec2d& x) {
  using namespace Common;

  static const vec2d vzero(0.0, 0.0);
  static const vec2d vec1((5.0 - sqrt(5.0)) / 10.0, (5.0 + sqrt(5.0)) / 10.0);
  static const vec2d vec2(1.0 / sqrt(5.0), -1.0 / sqrt(5.0));

  if (clipeps + checkPosition(vec2 + vzero, vec2 + vec1, x) < 0) return false;
  if (clipeps + checkPosition(-vec2 + vec1, -vec2 + vzero, x) < 0) return false;
  if (clipeps + checkPosition(-vec1 + vzero, -vec1 + vec2, x) < 0) return false;
  if (clipeps + checkPosition(vec1 + vec2, vec1 + vzero, x) < 0) return false;

  return true;
}

double ArithVisibility::VisOpZ2::angle(const vec2i& a) {
  return a.minkowskiZ2().angle();
}

vec2d ArithVisibility::VisOpZ2::toR2(const vec2i& a) {
  return a.minkowskiZ2();
}

bool ArithVisibility::VisOpZ2::rayTest(const vec2i& a, const vec2i& b) {
  if (a.x == 0) {
    return (b.x == 0);
  }

  if (b.x == 0) {
    return (a.x == 0);
  }

  if (a.y == 0) {
    return (b.y == 0);
  }

  if (b.y == 0) {
    return (a.y == 0);
  }

  return (b.x * a.y == b.y * a.x);
}

double ArithVisibility::VisOpGI::angle(const vec2i& a) {
  return a.angle();
}

vec2d ArithVisibility::VisOpGI::toR2(const vec2i& a) {
  return vec2d(double(a.x), double(a.y));
}

bool ArithVisibility::VisOpGI::rayTest(const vec2i& a, const vec2i& b) {
  if (a.x == 0) {
    return (b.x == 0);
  }

  if (b.x == 0) {
    return (a.x == 0);
  }

  if (a.y == 0) {
    return (b.y == 0);
  }

  if (b.y == 0) {
    return (a.y == 0);
  }

  return (b.x * a.y == b.y * a.x);
}

double ArithVisibility::VisOpES::angle(const vec2i& a) {
  return a.minkowskiES().angle();
}

vec2d ArithVisibility::VisOpES::toR2(const vec2i& a) {
  return a.minkowskiES();
}

bool ArithVisibility::VisOpES::rayTest(const vec2i& a, const vec2i& b) {
  if (a.x == 0) {
    return (b.x == 0);
  }

  if (b.x == 0) {
    return (a.x == 0);
  }

  if (a.y == 0) {
    return (b.y == 0);
  }

  if (b.y == 0) {
    return (a.y == 0);
  }

  return (b.x * a.y == b.y * a.x);
}

double ArithVisibility::VisOpGM::angle(const vec2i& a) {
  return a.minkowskiGM().angle();
}

vec2d ArithVisibility::VisOpGM::toR2(const vec2i& a) {
  return a.minkowskiGM();
}

bool ArithVisibility::VisOpGM::rayTest(const vec2i& a, const vec2i& b) {
  if (a.x == 0) {
    return (b.x == 0);
  }

  if (b.x == 0) {
    return (a.x == 0);
  }

  if (a.y == 0) {
    return (b.y == 0);
  }

  if (b.y == 0) {
    return (a.y == 0);
  }

  return (b.x * a.y == b.y * a.x);
}

void ArithVisibility::visCircleZ2(const uint r, Common::vec2ilist& out,
                            bool radialproj) {
  using namespace Common;

  vec2ilist circleZ2;
  vCircleZ2(r, circleZ2);

  cerr << "info: constructed patch of the Minkowski embedding of Z[Sqrt[2]] with "
       << circleZ2.size() << " vertices.\n";

  VisListZ2* vlist = new VisListZ2;
  vlist->reserve(lround(double(circleZ2.size()) * 0.71));

  vlist->init();

  for (vec2ilist::const_iterator i = circleZ2.begin(); i != circleZ2.end(); ++i) {
    if (visibility2FreeZ2(*i)) vlist->insertSorted(*i);
  }

  cerr << "info: isolated " << vlist->size()
       << " square-free elements of the patch.\n";

  circleZ2.resize(0); // deallocate the initial vertices

  if (radialproj)
    vlist->removeInvisibleFast();
  else
    vlist->removeInvisibleProper();

  out.clear();
  out.reserve(vlist->size());
  vlist->dump(out);

  delete vlist;
}

void ArithVisibility::visCircleZ2Fast(const uint r,
                            Common::vec2ilist& out, bool radialproj) {
  using namespace Common;

  vec2ilist circleZ2;
  vCircleZ2(r, circleZ2);

  cerr << "info: constructed patch of the Minkowski embedding of Z[Sqrt[2]] with "
       << circleZ2.size() << " vertices.\n";

  vec2ielist ext;
  uint num_sqfree;

  // estimated number of square-free elements
  const uint est_sqfree = lround(double(circleZ2.size()) * 0.71);

  out.clear();
  out.reserve(est_sqfree);

  if (radialproj) {
    // write directly into the output container
    for (vec2ilist::const_iterator i = circleZ2.begin(); i != circleZ2.end(); ++i)
      if (visibility2FreeZ2(*i)) out.push_back(i->primitive());

    num_sqfree = out.size();
  } else {
    ext.reserve(est_sqfree);

    for (vec2ilist::const_iterator i = circleZ2.begin(); i != circleZ2.end(); ++i)
      if (visibility2FreeZ2(*i)) ext.push_back(*i);

    num_sqfree = ext.size();
  }

  cerr << "info: isolated " << num_sqfree
       << " square-free elements of the patch.\n";

  circleZ2.resize(0);

  if (radialproj) {
    sort(out.begin(), out.end());
    out.erase(unique(out.begin(), out.end()), out.end());
  } else {
    sort(ext.begin(), ext.end());
    normalize(ext);
    ext.erase(unique(ext.begin(), ext.end()), ext.end());

    for (vec2ielist::const_iterator i = ext.begin(); i != ext.end(); ++i)
      out.push_back(*i);
  }

  cerr << "info: " << out.size() << " elements are visible.\n";
}

void ArithVisibility::visCircleGI(const uint r, Common::vec2ilist& out,
                            bool radialproj) {
  using namespace Common;

  vec2ilist circleGI;
  vCircleGI(r, circleGI);

  cerr << "info: constructed patch of the Gaussian Integers with "
       << circleGI.size() << " vertices.\n";

  VisListGI* vlist = new VisListGI;
  vlist->reserve(lround(double(circleGI.size()) * 0.69));

  vlist->init();

  for (vec2ilist::const_iterator i = circleGI.begin(); i != circleGI.end(); ++i) {
    if (visibility2FreeGI(*i)) vlist->insertSorted(*i);
  }

  cerr << "info: isolated " << vlist->size()
       << " square-free elements of the patch.\n";

  circleGI.resize(0); // deallocate the initial vertices

  if (radialproj)
    vlist->removeInvisibleFast();
  else
    vlist->removeInvisibleProper();

  out.clear();
  out.reserve(vlist->size());
  vlist->dump(out);

  delete vlist;
}

void ArithVisibility::visCircleGIFast(const uint r,
                            Common::vec2ilist& out, bool radialproj) {
  using namespace Common;

  vec2ilist circleGI;
  vCircleGI(r, circleGI);

  cerr << "info: constructed patch of the Gaussian Integers with "
       << circleGI.size() << " vertices.\n";

  vec2ielist ext;
  uint num_sqfree;

  // estimated number of square-free elements
  const uint est_sqfree = lround(double(circleGI.size()) * 0.69);

  out.clear();
  out.reserve(est_sqfree);

  if (radialproj) {
    // write directly into the output container
    for (vec2ilist::const_iterator i = circleGI.begin(); i != circleGI.end(); ++i)
      if (visibility2FreeGI(*i)) out.push_back(i->primitive());

    num_sqfree = out.size();
  } else {
    ext.reserve(est_sqfree);

    for (vec2ilist::const_iterator i = circleGI.begin(); i != circleGI.end(); ++i)
      if (visibility2FreeGI(*i)) ext.push_back(*i);

    num_sqfree = ext.size();
  }

  cerr << "info: isolated " << num_sqfree
       << " square-free elements of the patch.\n";

  circleGI.resize(0);

  if (radialproj) {
    sort(out.begin(), out.end());
    out.erase(unique(out.begin(), out.end()), out.end());
  } else {
    sort(ext.begin(), ext.end());
    normalize(ext);
    ext.erase(unique(ext.begin(), ext.end()), ext.end());

    for (vec2ielist::const_iterator i = ext.begin(); i != ext.end(); ++i)
      out.push_back(*i);
  }

  cerr << "info: " << out.size() << " elements are visible.\n";
}

void ArithVisibility::visCircleES(const uint r, Common::vec2ilist& out,
                            bool radialproj) {
  using namespace Common;

  vec2ilist circleES;
  vCircleES(r, circleES);

  cerr << "info: constructed patch of the Eisenstein Integers with "
       << circleES.size() << " vertices.\n";

  VisListES* vlist = new VisListES;
  vlist->reserve(lround(double(circleES.size()) * 0.8));

  vlist->init();

  for (vec2ilist::const_iterator i = circleES.begin(); i != circleES.end(); ++i) {
    if (visibility2FreeES(*i)) vlist->insertSorted(*i);
  }

  cerr << "info: isolated " << vlist->size()
       << " square-free elements of the patch.\n";

  circleES.resize(0); // deallocate the initial vertices

  if (radialproj)
    vlist->removeInvisibleFast();
  else
    vlist->removeInvisibleProper();

  out.clear();
  out.reserve(vlist->size());
  vlist->dump(out);

  delete vlist;
}

void ArithVisibility::visCircleESFast(const uint r,
                            Common::vec2ilist& out, bool radialproj) {
  using namespace Common;

  vec2ilist circleES;
  vCircleES(r, circleES);

  cerr << "info: constructed patch of the Eisenstein Integers with "
       << circleES.size() << " vertices.\n";

  vec2ielist ext;
  uint num_sqfree;

  // estimated number of square-free elements
  const uint est_sqfree = lround(double(circleES.size()) * 0.8);

  out.clear();
  out.reserve(est_sqfree);

  if (radialproj) {
    // write directly into the output container
    for (vec2ilist::const_iterator i = circleES.begin(); i != circleES.end(); ++i)
      if (visibility2FreeES(*i)) out.push_back(i->primitive());

    num_sqfree = out.size();
  } else {
    ext.reserve(est_sqfree);

    for (vec2ilist::const_iterator i = circleES.begin(); i != circleES.end(); ++i)
      if (visibility2FreeES(*i)) ext.push_back(*i);

    num_sqfree = ext.size();
  }

  cerr << "info: isolated " << num_sqfree
       << " square-free elements of the patch.\n";

  circleES.resize(0);

  if (radialproj) {
    sort(out.begin(), out.end());
    out.erase(unique(out.begin(), out.end()), out.end());
  } else {
    sort(ext.begin(), ext.end());
    normalize(ext);
    ext.erase(unique(ext.begin(), ext.end()), ext.end());

    for (vec2ielist::const_iterator i = ext.begin(); i != ext.end(); ++i)
      out.push_back(*i);
  }

  cerr << "info: " << out.size() << " elements are visible.\n";
}

void ArithVisibility::visCircleGM(const uint r, Common::vec2ilist& out,
                            bool radialproj) {
  using namespace Common;

  vec2ilist circleGM;
  vCircleGM(r, circleGM);

  cerr << "info: constructed patch of the Minkowski embedding of Z[tau] "
       << "(tau the golden mean) with " << circleGM.size()
       << " vertices.\n";

  VisListGM* vlist = new VisListGM;
  vlist->reserve(lround(double(circleGM.size()) * 0.89));

  vlist->init();

  for (vec2ilist::const_iterator i = circleGM.begin(); i != circleGM.end(); ++i) {
    if (visibility2FreeGM(*i)) vlist->insertSorted(*i);
  }

  cerr << "info: isolated " << vlist->size()
       << " square-free elements of the patch.\n";

  circleGM.resize(0); // deallocate the initial vertices

  if (radialproj)
    vlist->removeInvisibleFast();
  else
    vlist->removeInvisibleProper();

  out.clear();
  out.reserve(vlist->size());
  vlist->dump(out);

  delete vlist;
}

void ArithVisibility::visCircleGMFast(const uint r,
                            Common::vec2ilist& out, bool radialproj) {
  using namespace Common;

  vec2ilist circleGM;
  vCircleGM(r, circleGM);

  cerr << "info: constructed patch of the Minkowski embedding of Z[tau] with "
       << circleGM.size() << " vertices.\n";

  vec2ielist ext;
  uint num_sqfree;

  // estimated number of square-free elements
  const uint est_sqfree = lround(double(circleGM.size()) * 0.89);

  out.clear();
  out.reserve(est_sqfree);

  if (radialproj) {
    // write directly into the output container
    for (vec2ilist::const_iterator i = circleGM.begin(); i != circleGM.end(); ++i)
      if (visibility2FreeGM(*i)) out.push_back(i->primitive());

    num_sqfree = out.size();
  } else {
    ext.reserve(est_sqfree);

    for (vec2ilist::const_iterator i = circleGM.begin(); i != circleGM.end(); ++i)
      if (visibility2FreeGM(*i)) ext.push_back(*i);

    num_sqfree = ext.size();
  }

  cerr << "info: isolated " << num_sqfree
       << " square-free elements of the patch.\n";

  circleGM.resize(0);

  if (radialproj) {
    sort(out.begin(), out.end());
    out.erase(unique(out.begin(), out.end()), out.end());
  } else {
    sort(ext.begin(), ext.end());
    normalize(ext);
    ext.erase(unique(ext.begin(), ext.end()), ext.end());

    for (vec2ielist::const_iterator i = ext.begin(); i != ext.end(); ++i)
      out.push_back(*i);
  }

  cerr << "info: " << out.size() << " elements are visible.\n";
}

void ArithVisibility::radialProjZ2(const uint r, Common::dlist& out) {
  using namespace Common;

  vec2ilist circleZ2, sqfreeZ2;
  vCircleZ2(r, circleZ2);

  cerr << "info: constructed patch of the Minkowski embedding of Z[Sqrt[2]] with "
       << circleZ2.size() << " vertices.\n";

  sqfreeZ2.reserve(lround(double(circleZ2.size()) * 0.71));

  for (vec2ilist::const_iterator i = circleZ2.begin(); i != circleZ2.end(); ++i)
    if (visibility2FreeZ2(*i)) sqfreeZ2.push_back(i->primitive());

  cerr << "info: isolated " << sqfreeZ2.size()
       << " square-free elements of the patch.\n";

  circleZ2.resize(0); // deallocate the initial vertices

  sort(sqfreeZ2.begin(), sqfreeZ2.end());
  sqfreeZ2.erase(unique(sqfreeZ2.begin(), sqfreeZ2.end()), sqfreeZ2.end());

  cerr << "info: applying radial projection to " << sqfreeZ2.size()
       << " visible elements.\n";

  dlist angles;
  double meandist;

  angles.reserve(sqfreeZ2.size());
  for (vec2ilist::const_iterator i = sqfreeZ2.begin(); i != sqfreeZ2.end(); ++i)
    angles.push_back(i->minkowskiZ2().angle());

  sort(angles.begin(), angles.end());

  out.clear();
  out.reserve(angles.size() - 1);
  neighbourDiff(angles, out, meandist);
  normalizeAngDists(out, meandist);
}

void ArithVisibility::radialProjGI(const uint r, Common::dlist& out) {
  using namespace Common;

  vec2ilist circleGI, sqfreeGI;
  vCircleGI(r, circleGI);

  cerr << "info: constructed patch of the Gaussian Integers with "
       << circleGI.size() << " vertices.\n";

  sqfreeGI.reserve(lround(double(circleGI.size()) * 0.69));

  for (vec2ilist::const_iterator i = circleGI.begin(); i != circleGI.end(); ++i)
    if (visibility2FreeGI(*i)) sqfreeGI.push_back(i->primitive());

  cerr << "info: isolated " << sqfreeGI.size()
       << " square-free elements of the patch.\n";

  circleGI.resize(0); // deallocate the initial vertices

  sort(sqfreeGI.begin(), sqfreeGI.end());
  sqfreeGI.erase(unique(sqfreeGI.begin(), sqfreeGI.end()), sqfreeGI.end());

  cerr << "info: applying radial projection to " << sqfreeGI.size()
       << " visible elements.\n";

  dlist angles;
  double meandist;

  angles.reserve(sqfreeGI.size());
  for (vec2ilist::const_iterator i = sqfreeGI.begin(); i != sqfreeGI.end(); ++i)
    angles.push_back(i->angle());

  sort(angles.begin(), angles.end());

  out.clear();
  out.reserve(angles.size() - 1);
  neighbourDiff(angles, out, meandist);
  normalizeAngDists(out, meandist);
}

void ArithVisibility::radialProjES(const uint r, Common::dlist& out) {
  using namespace Common;

  vec2ilist circleES, sqfreeES;
  vCircleES(r, circleES);

  cerr << "info: constructed patch of the Eisenstein Integers with "
       << circleES.size() << " vertices.\n";

  sqfreeES.reserve(lround(double(circleES.size()) * 0.8));

  for (vec2ilist::const_iterator i = circleES.begin(); i != circleES.end(); ++i)
    if (visibility2FreeES(*i)) sqfreeES.push_back(i->primitive());

  cerr << "info: isolated " << sqfreeES.size()
       << " square-free elements of the patch.\n";

  circleES.resize(0); // deallocate the initial vertices

  sort(sqfreeES.begin(), sqfreeES.end());
  sqfreeES.erase(unique(sqfreeES.begin(), sqfreeES.end()), sqfreeES.end());

  cerr << "info: applying radial projection to " << sqfreeES.size()
       << " visible elements.\n";

  dlist angles;
  double meandist;

  angles.reserve(sqfreeES.size());
  for (vec2ilist::const_iterator i = sqfreeES.begin(); i != sqfreeES.end(); ++i)
    angles.push_back(i->minkowskiES().angle());

  sort(angles.begin(), angles.end());

  out.clear();
  out.reserve(angles.size() - 1);
  neighbourDiff(angles, out, meandist);
  normalizeAngDists(out, meandist);
}

void ArithVisibility::radialProjGM(const uint r, Common::dlist& out) {
  using namespace Common;

  vec2ilist circleGM, sqfreeGM;
  vCircleGM(r, circleGM);

  cerr << "info: constructed patch of the Minkowski embedding of Z[tau] with "
       << circleGM.size() << " vertices.\n";

  sqfreeGM.reserve(lround(double(circleGM.size()) * 0.89));

  for (vec2ilist::const_iterator i = circleGM.begin(); i != circleGM.end(); ++i)
    if (visibility2FreeGM(*i)) sqfreeGM.push_back(i->primitive());

  cerr << "info: isolated " << sqfreeGM.size()
       << " square-free elements of the patch.\n";

  circleGM.resize(0); // deallocate the initial vertices

  sort(sqfreeGM.begin(), sqfreeGM.end());
  sqfreeGM.erase(unique(sqfreeGM.begin(), sqfreeGM.end()), sqfreeGM.end());

  cerr << "info: applying radial projection to " << sqfreeGM.size()
       << " visible elements.\n";

  dlist angles;
  double meandist;

  angles.reserve(sqfreeGM.size());
  for (vec2ilist::const_iterator i = sqfreeGM.begin(); i != sqfreeGM.end(); ++i)
    angles.push_back(i->minkowskiGM().angle());

  sort(angles.begin(), angles.end());

  out.clear();
  out.reserve(angles.size() - 1);
  neighbourDiff(angles, out, meandist);
  normalizeAngDists(out, meandist);
}

/*
 * If we have an invertible (over Z) 2x2 matrix M, then M can only leave
 * the square-free Gaussian invariant, if the following property holds:
 *
 * Let n be a natural number, then for every 1 <= k <= n, the k-th
 * power of M should have square-free rows, i.e. write
 * M^k = {{a, b}, {c, d}}, then both a + c*I and b + d*I should be
 * square-free Gaussian integers.
 */
bool ArithVisibility::powerTest(const mtx2x2i& m, uint n) {
  mtx2x2i temp(mtx2x2i::identity);

  for (uint i = 0; i <= n; ++i) {
    temp = temp * m;

    if (!visibility2FreeGI(temp.column(0)) ||
        !visibility2FreeGI(temp.column(1)))
      return false;
  }

  return true;
}

/*
 * If we have an invertible (over Z) 2x2 matrix M, then M can only leave
 * the square-free Gaussian invariant, if the following property holds:
 *
 * For each vector (m, n) where m and n are coprime and m + n*I is
 * square-free, M multipied with (m,n) has to be square-free again.
 *
 * We usually create a large region of the square-free coprime Gaussians
 * and then apply the test for each element.
 */
bool ArithVisibility::multTest(const mtx2x2i& m, const vec2i& v) {
  const vec2i mv(m * v);

  return visibility2FreeGI(mv);
}

ostream& operator<<(ostream &os, const ArithVisibility::vec2iq& v) {
  os << '{' << v.getNumerator() << ',' << v.getDenominator() << '}';
  return os;
}

ostream& operator<<(ostream &os, const ArithVisibility::bragg& b) {
  os << "circle((" << b.getPosition()[0] << ',' << b.getPosition()[1] << "),"
     << b.getIntensity() << ",color='black',thickness=0.5)";
  return os;
}

void vTableZ2(const uint r, Common::vec2ilist& table) {
  const int temp = int(ceil(double(r) / sqrt(2.0)));
  table.clear();
  table.reserve((2*r + 1) * (2*temp + 1));

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
        const int g = gcdZFast(gcdZFast(abs(2 * i), abs(j)), 4 * c);

        table.push_back(vec2iq(2 * i / g, j / g, 4 * c / g));
      }
    }
  }

  sort(table.begin(), table.end());
  table.erase(unique(table.begin(), table.end()), table.end());
}

void vTableGI(const uint r, Common::vec2ilist& table) {
  table.clear();
  table.reserve((2*r + 1) * (2*r + 1));

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
        const int g = gcdZFast(gcdZFast(abs(i), abs(j)), c);

        table.push_back(vec2iq(i / g, j / g, c / g));
      }
    }
  }

  sort(table.begin(), table.end());
  table.erase(unique(table.begin(), table.end()), table.end());
}

void vTableES(const uint r, Common::vec2ilist& table) {
  table.clear();
  table.reserve((2*r + 1) * (2*r + 1));

  for (int i = -int(r); i <= int(r); ++i) {
    for (int j = -int(r); j <= int(r); ++j) {
      table.push_back(vec2i(i, j));
    }
  }
}

void vqTableRecipES(const uint r, const uint s,
                    vector<ArithVisibility::vec2iq>& table) {
  using namespace Coprime;
  using namespace ArithVisibility;

  table.clear();

  for (uint c = 1; c <= s; ++c) {
    for (int i = -int(r); i <= int(r); ++i) {
      for (int j = -int(r); j <= int(r); ++j) {
        const int g = gcdZFast(gcdZFast(abs(-2*i + 4*j), abs(-4*i + 2*j)), 3*c);

        table.push_back(vec2iq((-2*i + 4*j) / g, (-4*i + 2*j) / g, 3*c / g));
      }
    }
  }

  sort(table.begin(), table.end());
  table.erase(unique(table.begin(), table.end()), table.end());
}

void vTableGM(const uint r, Common::vec2ilist& table) {
  table.clear();
  table.reserve((2*r + 1) * (2*r + 1));

  for (int i = -int(r); i <= int(r); ++i) {
    for (int j = -int(r); j <= int(r); ++j) {
      table.push_back(vec2i(i, j));
    }
  }
}

void vqTableRecipGM(const uint r, const uint s,
                    vector<ArithVisibility::vec2iq>& table) {
  using namespace Coprime;
  using namespace ArithVisibility;

  table.clear();

  for (uint c = 1; c <= s; ++c) {
    for (int i = -int(r); i <= int(r); ++i) {
      for (int j = -int(r); j <= int(r); ++j) {
        const int g = gcdZFast(gcdZFast(abs(-i + 2*j), abs(2*i + j)), 5*c);

        table.push_back(vec2iq((-i + 2*j) / g, (2*i + j) / g, 5*c / g));
      }
    }
  }

  sort(table.begin(), table.end());
  table.erase(unique(table.begin(), table.end()), table.end());
}

void vCircleZ2(const uint r, Common::vec2ilist& table) {
  const int temp = int(ceil(double(r) / sqrt(2.0)));
  table.clear();
  table.reserve((2*r + 1) * (2*temp + 1));

  const double cradSq = double(r * r) * 2.0;

  for (int i = -int(r); i <= int(r); ++i) {
    for (int j = -temp; j <= temp; ++j) {
      const vec2i vtx(i, j);

      if (vtx.minkowskiZ2().lengthSquared() <= cradSq)
        table.push_back(vtx);
    }
  }
}

void vCircleGI(const uint r, Common::vec2ilist& table) {
  table.clear();
  table.reserve((2*r + 1) * (2*r + 1));

  const double cradSq = double(r * r);

  for (int i = -int(r); i <= int(r); ++i) {
    for (int j = -int(r); j <= int(r); ++j) {
      if (vec2d(double(i), double(j)).lengthSquared() <= cradSq)
        table.push_back(vec2i(i, j));
    }
  }
}

void vCircleES(const uint r, Common::vec2ilist& table) {
  table.clear();
  table.reserve((2*r + 1) * (2*r + 1));

  const double cradSq = double(r * r) * 0.75;

  for (int i = -int(r); i <= int(r); ++i) {
    for (int j = -int(r); j <= int(r); ++j) {
      const vec2i vtx(i, j);

      if (vtx.minkowskiES().lengthSquared() <= cradSq)
        table.push_back(vtx);
    }
  }
}

void vCircleGM(const uint r, Common::vec2ilist& table) {
  table.clear();
  table.reserve((2*r + 1) * (2*r + 1));

  /* Computation of the cut-off radius (geometry in parallelogram): *
   * R_max = Cos[ArcTan[1, 1] + ArcTan[tau, 1 - tau]] * Sqrt[2]     *
   * ArcTan[1, 1]: slope of first basis vector.                     *
   * ArcTan[tau, 1 - tau]: slope of second basis vector.            *
   * Sqrt[2]: length of first basis vector.                         */
  const double cradSq = double(r * r) * (5.0 / 3.0);

  for (int i = -int(r); i <= int(r); ++i) {
    for (int j = -int(r); j <= int(r); ++j) {
      const vec2i vtx(i, j);

      if (vtx.minkowskiGM().lengthSquared() <= cradSq)
        table.push_back(vtx);
    }
  }
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

    if (temp[0] < a[0]) {
      a[0] = temp[0];
    }

    if (temp[0] > b[0]) {
      b[0] = temp[0];
    }

    if (temp[1] < a[1]) {
      a[1] = temp[1];
    }

    if (temp[1] > b[1]) {
      b[1] = temp[1];
    }

    if (rtemp > r) r = rtemp;
  }

  min = a; max = b;
  radius = r;
}

void toEPS(const vector<ArithVisibility::bragg>& input, bool fill) {
  using namespace ArithVisibility;

  vec2d min, max;
  double radius;

  // Base width of the Postscript and offset to the borders.
  const int basewidth = 800;
  const int offset = 50;

  if (input.empty()) return;

  minmax(input, min, max, radius);

  // Compute range and midpoint of the input data
  const vec2d range(max - min);
  const vec2d mid(min + range * 0.5);

  // Aspect ratio used to compute height from basewidth
  const double aspect = (range[0] == 0.0) ? 1.0 : range[1] / range[0];
  const int height = lround(double(basewidth) * aspect);

  // Write header with bounding box information
  cout << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
  cout << "%%BoundingBox: 0 0 " << basewidth << ' ' << height << endl;
  cout << "%%BeginProlog" << endl;
  cout << "%%EndProlog" << endl;

  // Translate origin to the middle of the page
  cout << (basewidth / 2) << ' ' << (height / 2) << " translate" << endl;

  const double scaling = std::min(
    double(basewidth - offset) / (range[0] + radius),
    double(height - offset) / (range[1] + radius));

  // Apply scaling
  cout << lround(scaling) << ' ' << lround(scaling) << " scale" << endl;

  cout << "0.002 setlinewidth" << endl; // TODO: base width on min/max

  const string mode = (fill ? "fill" : "stroke");

  for (vector<bragg>::const_iterator k = input.begin(); k != input.end(); ++k) {
    const vec2d pos(k->getPosition() - mid);

    cout << "newpath" << endl;
    cout << pos[0] << ' ' << pos[1] << ' '
         << k->getIntensity() << " 0 360 arc" << endl;
    cout << mode << endl;
  }

  cout << "showpage" << endl;
  cout << "%%EOF" << endl;
}

void toEPS2(const vector<ArithVisibility::bragg>& input) {
  using namespace ArithVisibility;

  vec2d min, max;
  double radius;

  // Base width of the Postscript and offset to the borders.
  const int basewidth = 800;
  const int offset = 50;

  if (input.empty()) return;

  minmax(input, min, max, radius);

  // Compute range and midpoint of the input data
  const vec2d range(max - min);
  const vec2d mid(min + range * 0.5);

  // Aspect ratio used to compute height from basewidth
  const double aspect = (range[0] == 0.0) ? 1.0 : range[1] / range[0];
  const int height = lround(double(basewidth) * aspect);

  // Write header with bounding box information
  cout << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
  cout << "%%BoundingBox: 0 0 " << basewidth << ' ' << height << endl;
  cout << "%%BeginProlog" << endl;
  cout << "%%EndProlog" << endl;

  cout << "1 setlinecap" << endl
       << "1 setlinejoin" << endl;

  // circle hack: draw a line width linewidth set to radius and length zero
  cout << "/circlehack {2 mul setlinewidth moveto 0 0 rlineto stroke} def" << endl;

  // Translate origin to the middle of the page
  cout << (basewidth / 2) << ' ' << (height / 2) << " translate" << endl;

  const double scaling = std::min(
    double(basewidth - offset) / (range[0] + radius),
    double(height - offset) / (range[1] + radius));

  // Apply scaling
  cout << lround(scaling) << ' ' << lround(scaling) << " scale" << endl;

  for (vector<bragg>::const_iterator k = input.begin(); k != input.end(); ++k) {
    const vec2d pos(k->getPosition() - mid);

    cout << pos[0] << ' ' << pos[1] << ' '
         << k->getIntensity() << " circlehack" << endl;
  }

  cout << "showpage" << endl;
  cout << "%%EOF" << endl;
}

void toSVG(const vector<ArithVisibility::bragg>& input, bool fill) {
  using namespace ArithVisibility;

  vec2d min, max;
  double radius;

  // Base width of the SVG and offset to the borders.
  const int basewidth = 800;
  const int offset = 50;

  if (input.empty()) return;

  minmax(input, min, max, radius);

  // Compute range and midpoint of the input data
  const vec2d range(max - min);
  const vec2d mid(min + range * 0.5);

  // Aspect ratio used to compute height from basewidth
  const double aspect = (range[0] == 0.0) ? 1.0 : range[1] / range[0];
  const int height = lround(double(basewidth) * aspect);

  // Write header with geometry/dimension information
  cout << "<svg width=\"" << basewidth << "\" height=\"" << height
       << "\" xmlns=\"http://www.w3.org/2000/svg\" "
       << "version=\"1.1\">" << endl;

  // Translate origin to the middle of the page
  cout << "\t" << "<g transform=\"translate(" << (basewidth / 2) << ","
       << (height / 2) << ")\">" << endl;

  const double scaling = std::min(
    double(basewidth - offset) / (range[0] + radius),
    double(height - offset) / (range[1] + radius));

  // Apply scaling
  cout << "\t\t" << "<g transform=\"scale(" << lround(scaling) << ")\">" << endl;

  const string strwidth = "0.0004"; // TODO: base width on min/max
  const string fopacity = (fill ? "1.0" : "0.0");
  const string fcolor = (fill ? "255, 255, 255" : "0, 0, 0");

  for (vector<bragg>::const_iterator k = input.begin(); k != input.end(); ++k) {
    const vec2d pos(k->getPosition() - mid);

    cout << "\t\t\t" << "<circle cx=\"" << pos[0] << "\" cy=\""
         << pos[1] << "\" r=\"" << k->getIntensity()
         << "\" fill=\"rgb(" << fcolor << ")\" fill-opacity=\""
         << fopacity << "\" stroke-width=\"" << strwidth
         << "\" stroke=\"rgb(0, 0, 0)\" />" << endl;
  }

  cout << "\t\t" << "</g>" << endl
       << "\t" << "</g>" << endl;

  cout << "</svg>" << endl;
}

void exportRawConsole(const vector<ArithVisibility::bragg>& input) {
  using namespace ArithVisibility;

  vec2d min, max;
  double radius;

  minmax(input, min, max, radius);

  const double range[2] = {std::max(abs(min[0]), abs(max[0])),
                           std::max(abs(min[1]), abs(max[1]))};

  unsigned out;
  const char* data;

  /* write signature first */

  out = sizeof(bragg); /* element size */
  cout.write(reinterpret_cast<const char*>(&out), sizeof(unsigned) * 1);

  out = input.size();
  cout.write(reinterpret_cast<const char*>(&out), sizeof(unsigned) * 1);

  /* write range and radius */
  cout.write(reinterpret_cast<const char*>(range), sizeof(double) * 2);
  cout.write(reinterpret_cast<const char*>(&radius), sizeof(double) * 1);

  data = reinterpret_cast<const char*>(&(*input.begin()));
  cout.write(data, sizeof(bragg) * input.size());
}

double linscale(double x) {
  return x * 0.4;
}

double rootscale(double x) {
  return sqrt(x) * 0.1;
}

double powscale(double x) {
  return pow(x, 0.25) * 0.05;
}

typedef void (*tablefunc)(const uint, Common::vec2ilist&);
typedef bool (*visfunc)(const vec2i&);

int main_diffraction(int argc, char* argv[]) {
  using namespace ArithVisibility;

  cerr << "info: diffraction main mode selected." << endl;

  stringstream parser;

  uint mode = 0;
  uint num_range = 30;
  uint denom_range = 27;
  uint scale_idx = 0;

  vector<vec2iq> raster;
  vector<bragg> diffraction;

  tablefunc tfunc = NULL;
  visfunc vfunc = NULL;

  // Parse parameters passed to the application
  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> mode;

    if (argc >= 3) {
      parser.str(argv[2]);
      parser.clear();
      parser >> num_range;

      if (argc >= 4) {
        parser.str(argv[3]);
        parser.clear();
        parser >> denom_range;

        if (argc >= 5) {
          parser.str(argv[4]);
          parser.clear();
          parser >> scale_idx;
        }
      }
    }
  }

  // Clamp scaling function index and select the appropriate scaler
  if (scale_idx > 2) scale_idx = 2;
  const scalefunc sfunc = (scale_idx == 0 ? linscale :
                          (scale_idx == 1 ? rootscale : powscale));

  if (mode > 7) {
    cerr << "error: unknown mode (" << mode <<  ") selected.\n";
    return 1;
  }

  switch (mode) {
    case 0:
      vqTableRecipZ2(num_range, denom_range, raster);
      diffractionZ2(raster, diffraction, clipFundamentalZ2);
    break;

    case 1:
      tfunc = vTableZ2;
      vfunc = visibility2FreeZ2;
    break;

    case 2:
      vqTableRecipGI(num_range, denom_range, raster);
      diffractionGI(raster, diffraction, clipFundamentalGI);
    break;

    case 3:
      tfunc = vTableGI;
      vfunc = visibility2FreeGI;
    break;

    case 4:
      vqTableRecipES(num_range, denom_range, raster);
      diffractionES(raster, diffraction, clipFundamentalES);
    break;

    case 5:
      tfunc = vTableES;
      vfunc = visibility2FreeES;
    break;

    case 6:
      vqTableRecipGM(num_range, denom_range, raster);
      diffractionGM(raster, diffraction, clipFundamentalGM);
    break;

    case 7:
      tfunc = vTableGM;
      vfunc = visibility2FreeGM;
    break;

    default:
      assert(false);
    break;
  }

  /* Even modes are computing diffraction, so apply rescaling here *
   * and then export the pattern to the console.                   */
  if (mode % 2 == 0) {
    for (vector<bragg>::iterator k = diffraction.begin();
         k != diffraction.end(); ++k) {
      k->apply(sfunc);
    }

    //exportRawConsole(diffraction);
    toEPS2(diffraction);
  } else {
    // Density computation mode
    using namespace Common;

    vec2ilist table, sqfree;

    for (uint size = 50; size < 1000; size += 50) {
      tfunc(size, table);

      for (vec2ilist::const_iterator i = table.begin(); i != table.end(); ++i) {
        if (vfunc(*i)) sqfree.push_back(*i);
      }

      cout << size << ": ";
      cout << (double(sqfree.size()) / double(table.size())) << endl;
      sqfree.clear();
    }
  }

  return 0;
}

int main_radialproj(int argc, char* argv[]) {
  using namespace ArithVisibility;

  cerr << "info: radial projection main mode selected." << endl;

  stringstream parser;

  uint mode = 0;
  uint range = 30;

  Common::vec2ilist patch;
  Common::dlist spacings;

  // Parse parameters passed to the application
  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> mode;

    if (argc >= 3) {
      parser.str(argv[2]);
      parser.clear();
      parser >> range;
    }
  }

  if (mode > 7) {
    cerr << "error: unknown mode (" << mode <<  ") selected.\n";
    return 1;
  }

  switch (mode) {
    case 0:
      visCircleZ2Fast(range, patch, false);
    break;

    case 1:
      radialProjZ2(range, spacings);
    break;

    case 2:
      visCircleGIFast(range, patch, false);
    break;

    case 3:
      radialProjGI(range, spacings);
    break;

    case 4:
      visCircleESFast(range, patch, false);
    break;

    case 5:
      radialProjES(range, spacings);
    break;

    case 6:
      visCircleGMFast(range, patch, false);
    break;

    case 7:
      radialProjGM(range, spacings);
    break;

    default:
      assert(false);
    break;
  }

  if (mode % 2 == 0) {
    cout << patch << endl;
  } else {
    Common::writeRawConsole(spacings);
  }

  return 0;
}

int main_sqfree_grid(int argc, char* argv[]) {
  using namespace ArithVisibility;

  cerr << "info: square-free grid/bitmap main mode selected." << endl;

  stringstream parser;

  uint mode = 0;
  uint range = 30;

  Common::vec2ilist table, sqfree;

  tablefunc tfunc = NULL;
  visfunc vfunc = NULL;

  // Parse parameters passed to the application
  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> mode;

    if (argc >= 3) {
      parser.str(argv[2]);
      parser.clear();
      parser >> range;
    }
  }

  if (mode > 7) {
    cerr << "error: unknown mode (" << mode <<  ") selected.\n";
    return 1;
  }

  /*
   * Presenting the square-free integers in grid form (i.e. as bitmap)
   * might not make so such sense for the cases other than the
   * Gaussian integers.
   * We still allow the other cases here though.
   */
  switch (mode) {
    case 0:
    case 1:
      tfunc = vTableZ2;
      vfunc = visibility2FreeZ2;
    break;

    case 2:
    case 3:
      tfunc = vTableGI;
      vfunc = visibility2FreeGI;
    break;

    case 4:
    case 5:
      tfunc = vTableES;
      vfunc = visibility2FreeES;
    break;

    case 6:
    case 7:
      tfunc = vTableGM;
      vfunc = visibility2FreeGM;
    break;

    default:
      assert(false);
    break;
  }

  tfunc(range, table);

  for (Common::vec2ilist::const_iterator i = table.begin();
       i != table.end(); ++i) {
    if (vfunc(*i))
      sqfree.push_back(*i);
  }

  const OccupationArray occu(sqfree);
  occu.writePlainPBM();

  return 0;
}

/*
 * We know that the visible lattice points of Z^2 are invariant under
 * GL_2(Z). Now consider the square-free integers of Z[i] as a subset
 * of Z^2. Denote the subset as S. We can ask the same question here.
 * Which matrices leave S invariant? S is certainly going to be a
 * subgroup of GL_2(Z).
 *
 * To explore this question we check invertible matrices with components
 * in a certain range (here given by 'inv_range'). We apply the following
 * two tests:
 * (1) power test (see ArithVisibility::powerTest() for docu)
 * (2) multiplication test (see ArithVisibility::multTest() for docu)
 *
 * The powers used in (1) are configured through 'pow_range', the number
 * of elements used for (2) is configured through 'sqf_range'.
 */
int main_invariant_subgroup(int argc, char* argv[]) {
  using namespace ArithVisibility;

  cerr << "info: invariant subgroup main mode selected." << endl;

  stringstream parser;

  uint inv_range = 30, sqf_range = 40, pow_range = 20;

  Common::vec2ilist table, sampleCoprimeSqFree;
  vector<mtx2x2i> sampleGL2Z, output;

  // Parse parameters passed to the application
  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> inv_range;

    if (argc >= 3) {
      parser.str(argv[2]);
      parser.clear();
      parser >> sqf_range;

      if (argc >= 4) {
        parser.str(argv[3]);
        parser.clear();
        parser >> pow_range;
      }
    }
  }

  cerr << "info: input parameters = {" << inv_range << ", " << sqf_range
       << ", " << pow_range << '}' << endl;

  // Sample GL_2(Z)
  mtx2x2i::invertibleList(sampleGL2Z, inv_range);

  // Sample coprime square-free Gaussian integers
  vTableGI(sqf_range, table);
  for (Common::vec2ilist::const_iterator i = table.begin();
       i != table.end(); ++i) {
    if (visibility2FreeGI(*i) && i->coprime())
      sampleCoprimeSqFree.push_back(*i);
  }

  cerr << "info: number of matrices from GL_2(Z) = "
       << sampleGL2Z.size() << endl;
  cerr << "info: number of coprime square-free Gaussians = "
       << sampleCoprimeSqFree.size() << endl;

  for (vector<mtx2x2i>::const_iterator i = sampleGL2Z.begin();
       i != sampleGL2Z.end(); ++i) {
    bool failed = false;

    if (!ArithVisibility::powerTest(*i, pow_range))
      continue;

    for (Common::vec2ilist::const_iterator j = sampleCoprimeSqFree.begin();
         j != sampleCoprimeSqFree.end(); ++j) {
      if (!ArithVisibility::multTest(*i, *j)) {
        failed = true;
        break;
      }
    }

    if (!failed)
      output.push_back(*i);
  }

  cerr << "info: number of elements found = " << output.size() << endl;
  cout << "info: elements = " << output << endl;

  return 0;
}

int main_invariant_subgroup_chuck(int argc, char* argv[]) {
  using namespace ArithVisibility;

  cerr << "info: invariant subgroup main mode selected." << endl;

  stringstream parser;

  uint inv_range = 30, sqf_range = 40;

  Common::vec2ilist table, sampleSqFree;
  vector<mtx2x2i> sampleGL2Z, output;

  // Parse parameters passed to the application
  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> inv_range;

    if (argc >= 3) {
      parser.str(argv[2]);
      parser.clear();
      parser >> sqf_range;
    }
  }

  cerr << "info: input parameters = {" << inv_range << ", " << sqf_range
       << '}' << endl;

  // Sample GL_2(Z)
  mtx2x2i::invertibleList(sampleGL2Z, inv_range);

  // Sample coprime square-free Gaussian integers
  vTableGI(sqf_range, table);
  for (Common::vec2ilist::const_iterator i = table.begin();
       i != table.end(); ++i) {
    if (visibility2FreeGI(*i))
      sampleSqFree.push_back(*i);
  }

  cerr << "info: number of matrices from GL_2(Z) = "
       << sampleGL2Z.size() << endl;
  cerr << "info: number of square-free Gaussians = "
       << sampleSqFree.size() << endl;

  const vec2i splittings[] = {
    vec2i(-3,  4),
    vec2i(-3, -4)
  };
  const uint num_splitting = sizeof(splittings) / sizeof(splittings[0]);

  for (vector<mtx2x2i>::const_iterator i = sampleGL2Z.begin();
       i != sampleGL2Z.end(); ++i) {

    bool is_divisible = false;

    for (Common::vec2ilist::const_iterator j = sampleSqFree.begin();
         j != sampleSqFree.end(); ++j) {

      const vec2i Av((*i) * (*j));

      for (unsigned k = 0; k < num_splitting; ++k) {
        if (Av.isDivGI(splittings[k])) {
          is_divisible = true;
          break;
        }
      }

      if (is_divisible)
        break;
    }

    if (is_divisible)
      continue;

    output.push_back(*i);
  }

  cerr << "info: number of elements found = " << output.size() << endl;
  cout << "info: elements = " << output << endl;

  return 0;
}

void print_usage() {
  cerr << "arithmetic visibility: usage:" << endl;

  cerr << "arith_visibility --diffraction: selects diffraction main mode" << endl;
  cerr << "\t\tparameter 1: mode (even = diffraction; odd = density)" << endl;
    cerr << "\t\t\t" << "0/1 = Z[Sqrt[2]]; 2/3 = Gaussian Integers;" << endl;
    cerr << "\t\t\t" << "4/5 = Eisenstein Integers; 6/7 = Z[tau] (golden mean)" << endl;
  cerr << "\t\tparameter 2: numerator range" << endl;
  cerr << "\t\tparameter 3: denominator range" << endl;
  cerr << "\t\tparameter 4: scaling function "
       << "(0 = linear; 1 = square root; 2 = power)" << endl;

  cerr << "arith_visibility --radialproj: selects radial projection main mode" << endl;
  cerr << "\t\tparameter 1: mode (even = visible square-free; odd = radial projection)" << endl;
  cerr << "\t\tparameter 2: range (number of vertices depends on mode)" << endl;

  cerr << "arith_visibility --sqfree-grid: selects square-free grid/bitmap main mode" << endl;
  cerr << "\t\tparameter 1: mode (see diffraction main mode)" << endl;
  cerr << "\t\t\t(both even and odd map to the same mode)" << endl;
  cerr << "\t\tparameter 2: range (number of vertices depends on mode)" << endl;

  cerr << "arith_visibility --invariant-subgroup: selects invariant subgroup main mode" << endl;
  cerr << "\t\tparameter 1: sampling range of GL_2(Z) (default = 30)" << endl;
  cerr << "\t\tparameter 2: sampling range of square-free Gaussians (default = 40)" << endl;
  cerr << "\t\tparameter 3: range of powers applied to the matrices (default = 20)" << endl;

  cerr << "arith_visibility --invariant-subgroup-chuck: selects invariant subgroup (Christian) main mode" << endl;
  cerr << "\t\tparameter 1: sampling range of GL_2(Z) (default = 30)" << endl;
  cerr << "\t\tparameter 2: sampling range of square-free Gaussians (default = 40)" << endl;
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

  if (main_mode == "--diffraction")
    ret = main_diffraction(argc - 1, argv + 1);
  else if (main_mode == "--radialproj")
    ret = main_radialproj(argc - 1, argv + 1);
  else if (main_mode == "--sqfree-grid")
    ret = main_sqfree_grid(argc - 1, argv + 1);
  else if (main_mode == "--invariant-subgroup")
    ret = main_invariant_subgroup(argc - 1, argv + 1);
  else if (main_mode == "--invariant-subgroup-chuck")
    ret = main_invariant_subgroup_chuck(argc - 1, argv + 1);
  else
    print_usage();

  return ret;
}
