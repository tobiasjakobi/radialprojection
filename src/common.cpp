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

#include "common.h"

#include <fstream>
#include <algorithm>

vec2d vec4i::shift(0.0, 0.0);

const double vec2d::sectorL3 = sqrt(3.0); /* = tan(2*pi/6) */
const double vec2d::sectorL5 = sqrt(5.0 + 2.0*sqrt(5.0)); /* = tan(2*pi/5) */
const double vec2d::sectorL12 = sqrt(3.0); /* = tan(2*pi/6) */

// Used for the higher-order cyclotomic cases:
const double vec2d::sectorL7 = tan(2.0 * Common::pi / 7.0);

double Common::RadiusSelector::radiusSq = 0.0;

vec2i vec2i::multUnitGM(int k) const {
  int a = x;
  int b = y;

  if (k != 0) {
    if (k > 0) {
      while (k != 0) {
        const int temp = a;
        a = b;
        b = temp + b;
        --k;
      }
    } else {
      while (k != 0) {
        const int temp = a;
        a = b - a;
        b = temp;
        ++k;
      }
    }
  }

  return vec2i(a, b);
}

vec2i vec2i::reduceGM(int& k) const {
  static const double tau = (1.0 + sqrt(5.0)) * 0.5;

  // Nothing to reduce if the element is already zero.
  if (x == 0 && y == 0) {
    k = 0;
    return *this;
  }

  int a, b, t;
  bool sign = false;
  double z = double(x) + double(y) * tau;

  if (z < 0.0) {
    a = -x;
    b = -y;
    z = -z;
    sign = true;
  } else {
    a = x;
    b = y;
  }

  /* If the absolute value of the algebraic norm is one, then we're already        *
   * dealing with a unit. In this case the log-fraction is already an integer.     *
   * Use 'lround' instead of 'floor' in this case, to avoid numerical instability. */
  if (this->normZTau() == 1)
    t = lround(log(z) / log(tau));
  else
    t = floor(log(z) / log(tau));

  assert(double(t) * log(tau) <= log(z) + 0.0001);
  assert(log(z) < double(t+1) * log(tau));

  k = -t;

  /* Move the element into the interval [1, tau) through *
   * multiplication with either tau or tau^{-1}.         */
  if (t != 0) {
    if (t < 0) {
      // Multiply with tau
      while (t != 0) {
        const int temp = a;
        a = b;
        b = temp + b;
        ++t;
      }
    } else {
      // Divide by tau
      while (t != 0) {
        const int temp = a;
        a = b - a;
        b = temp;
        --t;
      }
    }
  }

  assert(1.0 <= double(a) + double(b) * tau);
  assert(double(a) + double(b) * tau < tau);

  if (sign) {
    a = -a;
    b = -b;
  }

  return vec2i(a, b);
}

vec2d vec2i::transTriToR2() const {
  static const double factor = sqrt(3.0) * 0.5;

  return vec2d(
    double(x) + double(y) * 0.5,
    double(y) * factor);
}

vec2d vec2i::transGenericToR2(const vec2d& v) const {
  return vec2d(
    double(x) + v.x * double(y),
    v.y * double(y));
}

vec2d vec2i::minkowskiZ2() const {
  return vec2d(
    double(x) + double(y) * sqrt(2.0),
    double(x) - double(y) * sqrt(2.0)
  );
}

vec2d vec2i::minkowskiES() const {
  static const double factor = sqrt(3.0) * 0.5;

  return vec2d(
    double(x) - double(y) * 0.5,
    double(y) * factor);
}

vec2d vec2i::minkowskiGM() const {
  static const double tau = (1.0 + sqrt(5.0)) * 0.5;

  return vec2d(
    double(x) + double(y) * tau,
    double(x) + double(y) * (1.0 - tau)
  );
}

bool vec2i::coprime() const {
  return (Coprime::gcdZFast(abs(x), abs(y)) == 1);
}

vec2i vec2i::primitive() const {
  const int g = Coprime::gcdZFast(abs(x), abs(y));

  return vec2i(x / g, y / g);
}

vec4s vec4s::directL10ToUnique() const {
  if (this->isZero()) return *this;

  const vec2i x(a[0], a[1]);
  const vec2i y(a[2], a[3]);

  // Special case handling
  if (x.isZero()) return vec4s(0, 0, y.isPositiveGM() ? 1 : -1, 0);
  if (y.isZero()) return vec4s(x.isPositiveGM() ? 1 : -1, 0, 0, 0);

  // Compute a positive GCD in Z[tau]
  const vec2i g(Coprime::gcdZTau(x, y).positiveGM());

  int k;

  /* Divide the coordinates by g (making them primitive) *
   * and then apply reduction, making (rx, ry) unique.   */
  const vec2i rx(x.divGM(g).reduceGM(k));
  const vec2i ry(y.divGM(g).multUnitGM(k));

  return vec4s(rx, ry);
}

vec2iExt::vec2iExt(const vec2i& in) {
  e = Coprime::gcdZFast(abs(in.x), abs(in.y));
  v.x = in.x / e;
  v.y = in.y / e;
}

vec2iExt::vec2iExt(const vec2iExt& c) {
  v = c.v;
  e = c.e;
}

bool vec2iExt::operator==(const vec2iExt& rhs) const {
  return (v == rhs.v);
}

bool vec2iExt::operator<(const vec2iExt& rhs) const {
  return (v < rhs.v);
}

vec2iExt::operator vec2i() const {
  return vec2i(v.x * e, v.y * e);
}

vec2iExt::operator int() const {
  return e;
}

uint Coprime::gcdZFast(uint u, uint v) {
  uint shift;

  // GCD(0,v) == v; GCD(u,0) == u, GCD(0,0) == 0
  if (u == 0) return v;
  if (v == 0) return u;

  /* Let shift := lg K, where K is the greatest power of 2
     dividing both u and v. */
  for (shift = 0; ((u | v) & 1) == 0; ++shift) {
    u >>= 1;
    v >>= 1;
  }

  while ((u & 1) == 0)
    u >>= 1;

  // From here on, u is always odd.
  do {
    /* remove all factors of 2 in v -- they are not common
       note: v is not zero, so while will terminate */
    while ((v & 1) == 0)  /* Loop X */
      v >>= 1;
 
    /* Now u and v are both odd. Swap if necessary so u <= v,
       then set v = v - u (which is even). For bignums, the
       swapping is just pointer movement, and the subtraction
       can be done in-place. */
    if (u > v) {
      // Swap u and v.
      uint t = v;
      v = u;
      u = t;
    }

    v = v - u; // Here v >= u.
  } while (v != 0);
 
  // restore common factors of 2
  return u << shift;
}

int Coprime::gcdZ(int a, int b) {
  int x = a, y = b, z;

  while (y != 0) {
    z = y;
    y = x % y;
    x = z;
  }

  return x;
}

void Coprime::gcdZTest(uint count) {
  srand(time(NULL));

  for (uint i = 0; i < count; ++i) {
    const uint test[2] = {abs(rand() % 65536), abs(rand() % 65536)};

    const uint ref = abs(gcdZ(int(test[0]), int(test[1])));
    const uint fast = gcdZFast(test[0], test[1]);

    if (ref != fast) {
      cerr << "error: gcd test in Z failed:" << endl;
      cerr << "reference = " << ref << "; fast = " << fast << endl;
      break;
    }
  }
}

vec2i Coprime::moduloZ2(const vec2i& a, const vec2i& b) {
  const double N = 1.0 / double(b.preNormZ2());

  const int alpha = lround(double(a.x * b.x - 2 * a.y * b.y) * N);
  const int beta = lround(double(a.y * b.x - a.x * b.y) * N);

  return vec2i(a.x - (alpha * b.x + 2 * beta * b.y),
               a.y - (alpha * b.y + beta * b.x));
}

vec2i Coprime::gcdZ2(const vec2i& a, const vec2i& b) {
  vec2i x(a);
  vec2i y(b);
  vec2i z;

  while (!y.isZero()) {
    z = y;
    y = moduloZ2(x, y);
    x = z;
  }

  return x;
}

/* a = a.x + a.y * sqrt[2]                              *
 * b = b.x + b.y * sqrt[2]                              *
 * This computes c=a*b in terms of c.x + c.y * sqrt[2]  *
 * and stores the result in out.                        */
void Coprime::multZ2(const vec2i& a, const vec2i& b, vec2i& out) {
  out.set(a.x * b.x + 2 * a.y * b.y,
          a.x * b.y + a.y * b.x);
}

vec2i Coprime::moduloZTau(const vec2i& a, const vec2i& b) {
  const double N = 1.0 / double(b.preNormZTau());

  const int alpha = lround(double(a.x*b.x - a.y*b.y + a.x*b.y) * N);
  const int beta = lround(double(a.y*b.x - a.x*b.y) * N);

  return vec2i(a.x - (alpha*b.x + beta*b.y),
               a.y - (alpha*b.y + beta * (b.x + b.y)));
}

vec2i Coprime::gcdZTau(const vec2i& a, const vec2i& b) {
  vec2i x(a);
  vec2i y(b);
  vec2i z;

  while (!y.isZero()) {
    z = y;
    y = moduloZTau(x, y);
    x = z;
  }

  return x;
}

/* a = a.x + a.y * tau                              *
 * b = b.x + b.y * tau (with tau the golden mean)   *
 * This computes c=a*b in terms of c.x + c.y * tau  *
 * and stores the result in out.                    */
void Coprime::multZTau(const vec2i& a, const vec2i& b, vec2i& out) {
  out.set(a.x * b.x + a.y * b.y,
          a.y * (b.x + b.y) + a.x * b.y);
}

vec2i Coprime::moduloZ3(const vec2i& a, const vec2i& b) {
  const double N = 1.0 / double(b.preNormZ3());

  const int alpha = lround(double(a.x * b.x - 3 * a.y * b.y) * N);
  const int beta = lround(double(a.y * b.x - a.x * b.y) * N);

  return vec2i(a.x - (alpha * b.x + 3 * beta * b.y),
               a.y - (alpha * b.y + beta * b.x));
}

vec2i Coprime::gcdZ3(const vec2i& a, const vec2i& b) {
  vec2i x(a);
  vec2i y(b);
  vec2i z;

  while (!y.isZero()) {
    z = y;
    y = moduloZ3(x, y);
    x = z;
  }

  return x;
}

/* a = a.x + a.y * sqrt[3]                              *
 * b = b.x + b.y * sqrt[3]                              *
 * this computes c=a*b in terms of c.x + c.y * sqrt[3]  *
 * and stores the result in out.                        */
void Coprime::multZ3(const vec2i& a, const vec2i& b, vec2i& out) {
  out.set(a.x * b.x + 3 * a.y * b.y,
          a.x * b.y + a.y * b.x);
}

vec2i Coprime::moduloGI(const vec2i& a, const vec2i& b) {
  const double N = 1.0 / double(b.preNormGI());

  const int alpha = lround(double(a.x * b.x + a.y * b.y) * N);
  const int beta = lround(double(a.y * b.x - a.x * b.y) * N);

  return vec2i(a.x - (alpha * b.x - beta * b.y),
               a.y - (alpha * b.y + beta * b.x));
}


vec2i Coprime::gcdGI(const vec2i& a, const vec2i& b) {
  vec2i x(a);
  vec2i y(b);
  vec2i z;

  while (!y.isZero()) {
    z = y;
    y = moduloGI(x, y);
    x = z;
  }

  return x;
}

void Coprime::multGI(const vec2i& a, const vec2i& b, vec2i& out) {
  out.set(a.x * b.x - a.y * b.y,
          a.x * b.y + a.y * b.x);
}

vec2i Coprime::moduloES(const vec2i& a, const vec2i& b) {
  const double N = 1.0 / double(b.preNormES());

  const int alpha = lround(double(a.x * b.x + a.y * b.y - a.x * b.y) * N);
  const int beta = lround(double(a.y * b.x - a.x * b.y) * N);

  return vec2i(a.x - (alpha * b.x - beta * b.y),
               a.y - (alpha * b.y + beta * (b.x - b.y)));
}

vec2i Coprime::gcdES(const vec2i& a, const vec2i& b) {
  vec2i x(a);
  vec2i y(b);
  vec2i z;

  while (!y.isZero()) {
    z = y;
    y = moduloES(x, y);
    x = z;
  }

  return x;
}

void Coprime::multES(const vec2i& a, const vec2i& b, vec2i& out) {
  out.set(a.x * b.x - a.y * b.y,
          a.x * b.y + a.y * (b.x - b.y));
}

void Common::normalize(vec2ielist& in) {
  if (in.size() < 2) return;

  vec2ielist::iterator i = in.begin();
  vec2ielist::iterator j = i;
  ++i;

  bool seq = false;
  int min = *j;

  for (; i != in.end(); ++i) {
    if (*i == *j) {
      seq = true;
      if (*i < min) min = *i;

      continue;
    }

    if (seq) {
      while (true) {
        j->set(min);

        ++j;
        if (j == i) break;
      }
    }

    j = i;
    min = *j;
    seq = false;
  }
}

double Common::checkPosition(const vec2d& a, const vec2d& b, const vec2d& c) {
  const double pos = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  return pos;
}

bool Common::circularCheck(double radSquared, double xSquared) {
  if (radSquared - xSquared > eps) {
    return true;
  } else {
    if (radSquared - xSquared < -eps) return false;
  }

  cerr << "Warning: Insufficient accuracy in function circularCheck.\n";
  return false;
}

void Common::srandExt() {
  uint seed;

  ifstream urandom("/dev/urandom", ios::in|ios::binary);
  urandom.read(reinterpret_cast<char*>(&seed), sizeof(uint));
  urandom.close();

  srand(seed);
}

// Create n random double floats in the range [0.0, 1.0]
// Warning: This call is slow and not suited for large amount of random numbers.
//          Better use srand and rand in this context!
void Common::random(uint n, double out[]) {
  if (n == 0) return;

  uint* mem = new uint[n];

  // Use /dev/urandom as RNG (better than the default C/C++ RNG)
  ifstream urandom("/dev/urandom", ios::in|ios::binary);
  urandom.read(reinterpret_cast<char*>(mem), n * sizeof(uint));
  urandom.close();  

  for (uint i = 0; i < n; ++i) {
    out[i] = double(mem[i]) / double(numeric_limits<uint>::max());
  }
  delete [] mem;
}

// Same as random, but creates data in the range [0.0, range]
void Common::random(uint n, double range, double out[]) {
  if (n == 0) return;

  uint* mem = new uint[n];
  const double scaler = range / double(numeric_limits<uint>::max());

  ifstream urandom("/dev/urandom", ios::in|ios::binary);
  urandom.read(reinterpret_cast<char*>(mem), n * sizeof(uint));
  urandom.close();  

  for (uint i = 0; i < n; ++i) {
    out[i] = double(mem[i]) * scaler;
  }
  delete [] mem;
}

// VERY naive implementation of a prime number test
bool Common::isprime(uint x) {
  if (x == 0 || x == 1) return false;

  bool prime = true;
  for (uint i = 2; i < x; ++i) {
    if (x % i == 0) {
      prime = false;
      break;
    }
  }

  return prime;
}

uint Common::eulerPhi(uint i) {
  if (i == 1) return 1;

  uint res = i;

  // Check for divisibility by every prime number below the square root.
  // Start with 2.
  if (i % 2 == 0) {
    res -= (res / 2);
    do i /= 2;
    while (i % 2 == 0);
  }

  // Since this doesn't use a list of primes, check every odd number. Ideally, skip past composite numbers.
  for (uint j = 3; j*j <= i; j += 2) {
    if (i % j == 0) {
      res -= (res / j);
      do i /= j;
      while (i % j == 0);
    }
  }

  // If i>1, then it's the last factor at this point.
  if (i > 1) res -= (res / i);

  return res; 
}

double Common::power(double x, uint y) {
  if (y == 0) return 1.0;

  double temp = 1.0;
  for (uint i = 0; i < y; ++i) {
    temp *= x;
  }

  return temp;
}

int Common::ipower(int x, uint y) {
  if (y == 0) return 1;

  int temp = 1;
  for (uint i = 0; i < y; ++i) {
    temp *= x;
  }

  return temp;
}

void Common::readRawConsole(dlist& output) {
  double temp;
  while (true) {
    cin.read(reinterpret_cast<char*>(&temp), sizeof(double));
    if (cin.eof()) break;

    output.push_back(temp);
  }
}

void Common::writeRawConsole(const vec4ilist& input) {
  unsigned out;
  const char* data;

  /* write signature first */

  out = sizeof(int); /* element size */
  cout.write(reinterpret_cast<const char*>(&out), sizeof(unsigned) * 1);

  out = 4; /* number of elements per entry */
  cout.write(reinterpret_cast<const char*>(&out), sizeof(unsigned) * 1);

  out = input.size();
  cout.write(reinterpret_cast<const char*>(&out), sizeof(unsigned) * 1);

  data = reinterpret_cast<const char*>(&(*input.begin()));
  cout.write(data, sizeof(vec4i) * input.size());
}

void Common::readRawConsole(vec4ilist& output) {
  unsigned in;
  vec4i data;

  cin.read(reinterpret_cast<char*>(&in), sizeof(unsigned));
  if (cin.eof() && cin.fail()) goto readfail;
  if (in != sizeof(int)) goto signfail;

  cin.read(reinterpret_cast<char*>(&in), sizeof(unsigned));
  if (cin.eof() && cin.fail()) goto readfail;
  if (in != 4) goto signfail;

  cin.read(reinterpret_cast<char*>(&in), sizeof(unsigned));
  if (cin.eof() && cin.fail()) goto readfail;
    
  output.reserve(output.size() + in);

  while (true) {
    cin.read(reinterpret_cast<char*>(&data), sizeof(vec4i));
    if (cin.eof()) break;

    output.push_back(data);
  }

  return;
signfail:
  cerr << "error: verifying signature failed\n";
  return ;

readfail:
  cerr << "error: reading signature failed\n";
  return ;
}

void Common::minmax(const dlist& input, double& min, double& max) {
  double a, b;

  dlist::const_iterator i = input.begin();
  a = b = *i;

  ++i;
  for (; i != input.end(); ++i) {
    if (*i < a) a = *i;
    if (*i > b) b = *i;
  }

  min = a;
  max = b;
}

void Common::binmax(const BinningData& input, uint& index) {
  uint idx = 0;

  const vector<uint>& data = input.data;

  for (uint i = 0; i < input.data.size(); ++i) {
    if (data[i] > data[idx]) idx = i;
  }

  index = idx;
}

void Common::printstats(const dlist& data, const BinningData& bin) {
  double min, max;
  minmax(data, min, max);

  cerr << "info: minimum input value = " << min << endl;
  cerr << "info: maximum input value = " << max << endl;

  uint index;
  binmax(bin, index);

  cerr << "info: largest bin with " << bin.data[index] << " elements has index "
       << index << endl;
  cerr << "info: bin = [" << bin.range[0] + bin.step * double(index)
       << ", " << bin.range[0] + bin.step * double(index + 1)
       << "] (midpoint = " << bin.range[0] + bin.step * (double(index) + 0.5)
       << ")\n";
}

void Common::histogramBinning(const dlist& input, BinningData& output) {
  if (input.empty()) return;

  uint catched = 0;
  double input_max;

  const double a = output.range[0];
  const double step = output.step;

  if (output.tail) {
    // find maximum value of input list
    dlist::const_iterator i = input.begin();
    input_max = *i; ++i;
    for (; i != input.end(); ++i) {
      if (*i > input_max) input_max = *i;
    }

    // check if we don't overflow the (maximum) bin count
    if ((input_max + 1.0 - a) / step > double(numeric_limits<uint>::max()))
      return;

    cerr << "info: computing histogram binning for data tail "
         << "(maximum input value = " << input_max << ")\n";
  }

  const double b = (output.tail ? input_max + 1.0 : output.range[1]);

  vector<uint>& data = output.data;
  data.clear();
  data.resize(uint((b - a) / step), 0);

  for (dlist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const double cur = *i;

    if (cur < a || cur >= b) continue;

    ++data[uint((cur - a) / step)];
    ++catched;
  }

  cerr << "info: computed histogram with " << data.size() << " bins (interval = ["
       << a << ',' << b << "); step width = " << step << ")\n";
  cerr << "statistics: " << catched << " data points (from " << input.size()
       << ") fall into the binning area\n";

  output.catched = catched;
}

template <typename T>
void Common::histogramScale(const BinningData& input,
                            vector<T>& output, T scale) {
  output.clear();
  output.reserve(input.data.size());

  for (vector<uint>::const_iterator i = input.data.begin();
       i != input.data.end(); ++i) {
    output.push_back(T(*i) * scale);
  }
}

/* Creates "envelope" data for given histogram input:                 *
 * Can be used for ListPlot to visualize the distributions coming     *
 * from numerical simulations of the radial projection.               */
void Common::histogramEnvelope(double a, double b, double step, bool stats) {
  // fetch input data from console
  dlist inputData;
  readRawConsole(inputData);

  // compute histogram binning
  BinningData binData;
  binData.range[0] = a;
  binData.range[1] = b;
  binData.step = step;
  binData.tail = false;
  histogramBinning(inputData, binData);

  // convert to envelope by applying scaling
  dlist envelopeData;
  histogramScale(binData, envelopeData, 1.0 / (double(inputData.size()) * step));

  if (stats) printstats(inputData, binData);

  // output to console
  writeRawConsole(envelopeData);
}

// Same as histogramEnvelope but processes the "tail" of the data.
void Common::histoTailEnvelope(double a, double step, bool stats) {
  dlist inputData;
  readRawConsole(inputData);

  BinningData binData;
  binData.range[0] = a;
  binData.step = step;
  binData.tail = true;
  histogramBinning(inputData, binData);

  if (binData.data.empty()) {
    cerr << "error: binning failed (range or precision issues?)\n";
    return;
  }

  // Since the tail is usually quite long, we need extended precision here
  eflist envelopeData;
  histogramScale(binData, envelopeData,
    1.0L / (long double)(double(inputData.size()) * step));

  if (stats) printstats(inputData, binData);

  writeRawConsole(envelopeData);
}

void Common::histogramEnvelopeLD(double a, double b, double step) {
  const uint num_bin = uint((b - a) / step);

  ullong ndata = 0;
  ullong in_bin = 0;
  ullong* bins = new ullong[num_bin];

  for (uint i = 0; i < num_bin; ++i)
    bins[i] = 0;

  while (true) {
    double data;

    cin.read(reinterpret_cast<char*>(&data), sizeof(double));
    ++ndata;

    if (data >= a && data < b) {
      ++bins[uint((data - a) / step)];
      ++in_bin;
    }

    if (cin.eof()) break;
  }

  dlist envelopeData;

  cerr << "Computing histogram with " << num_bin << " bins (interval = ["
       << a << ',' << b << "); step width = " << step << ")\n";
  cerr << "statistics: " << in_bin << " data points (from " << ndata
       << ") fall into the binning area\n";

  const __float128 scaler = __float128(1.0) / (__float128(ndata) * __float128(step));
  for (uint i = 0; i < num_bin; ++i) {
    const __float128 t = __float128(bins[i]) * scaler;

    envelopeData.push_back(t);
  }
  delete [] bins;

  writeRawConsole(envelopeData);
}

void Common::histoTailEnvelopeLD(double a, double step) {
  ullong ndata = 0;
  ullong in_bin = 0;
  uint lindex = 0;  

  vector<ullong> bins;
  bins.resize(1);

  while (true) {
    double data;

    cin.read(reinterpret_cast<char*>(&data), sizeof(double));
    ++ndata;

    if (data >= a) {
      const uint bindex = uint((data - a) / step);

      if (bindex > lindex) {
        bins.resize(bindex + 1);
        lindex = bindex;
      }

      ++bins[bindex];
      ++in_bin;
    }

    if (cin.eof()) break;
  }

  eflist envelopeData;

  cerr << "Computing histogram with " << (lindex + 1) << " bins (tail beginning from "
       << a << "; step width = " << step << ")\n";
  cerr << "statistics: " << in_bin << " data points (from " << ndata
       << ") fall into the binning area\n";
  

  const __float128 scaler = __float128(1.0) / (__float128(ndata) * __float128(step));
  for (uint i = 0; i < lindex + 1; ++i) {
    const __float128 t = __float128(bins[i]) * scaler;

    envelopeData.push_back(t);
  }

  writeRawConsole(envelopeData);
}

void Common::binningCopySettings(const BinningStats& input, BinningData& output) {
  output.range[0] = input.range[0];
  output.range[1] = input.range[1];
  output.step = input.step;
  output.tail = input.tail;
}

void Common::histogramStatistics(const dlist& input, BinningStats& output) {
  BinningData binData;
  binningCopySettings(output, binData);

  histogramBinning(input, binData);

  dlist envelopeData;
  histogramScale(binData, envelopeData,
    1.0 / (double(input.size()) * binData.step));

  minmax(input, output.min, output.max);
  binmax(binData, output.maxbin_index);

  output.maxbin_position = binData.range[0] +
    binData.step * (double(output.maxbin_index) + 0.5);
}

void Common::neighbourDiff(const dlist& input, dlist& output, double& mean) {
  assert(input.size() != 0 && input.size() != 1);

  dlist::const_iterator j = input.begin();
  mean = *j;

  dlist::const_iterator k = j;
  ++k;
  while (k != input.end()) {
    output.push_back(*k - *j);
    ++j;
    ++k;
  }

  --j;
  mean = (*j - mean) / double(input.size() - 1);
}

void Common::normalizeAngDists(dlist& input, double mean) {
  const double normalizer = 1.0 / mean;

  for (dlist::iterator i = input.begin(); i != input.end(); ++i) {
    (*i) *= normalizer;
  }
}

void Common::radialProj(const vec2dlist& input,
                dlist& output, double& meandist) {
  output.clear();
  output.reserve(input.size());

  dlist angles;
  angles.reserve(input.size());

  for (vec2dlist::const_iterator i = input.begin(); i != input.end(); ++i)
    angles.push_back(i->angle());

  sort(angles.begin(), angles.end());
  neighbourDiff(angles, output, meandist);
  normalizeAngDists(output, meandist);
}

ostream& operator<<(ostream &os, const vec2i& v) {
  os << '{' << v.x << ',' << v.y << '}';
  return os;
}

ostream& operator<<(ostream &os, const vec2d& v) {
  os << '{' << fixed << setprecision(3)
     << v.x << ',' << v.y << '}';
  return os;
}

ostream& operator<<(ostream &os, const vec4i& v) {
  os << '{' << v[0] << ',' << v[1] << ',' << v[2] << ',' << v[3] << '}';
  return os;
}

ostream& operator<<(ostream &os, const vec8s& v) {
  os << '{' << int(v[0]);

  for (uint i = 1; i < 8; ++i) {
    os << ',' << int(v[i]);
  }

  os << '}';

  return os;
}

ostream& operator<<(ostream &os, const vec4s& v) {
  os << '{' << int(v[0]) << ',' << int(v[1])
     << ',' << int(v[2]) << ',' << int(v[3])
     << '}';

  return os;
}

ostream& operator<<(ostream &os, const vec2s& v) {
  os << '{' << int(v[0])
     << ',' << int(v[1]) << '}';

  return os;
}

ostream& operator<<(ostream &os, const vec2iExt& rhs) {
  os << '{' << rhs.get() << ',' << (int)rhs << '}';

  return os;
}

ostream& operator<<(ostream &os, const tilingEdge& e) {
  os << '{' << e[0]
     << ',' << e[1] << '}';

  return os;
}

