#include "griddual.h"

#include <sstream>
#include <algorithm>

namespace GridDualizing {

  void createGammaRandom(uint n, Common::dlist& output);
  void createGammaRandom(uint n, double sum, Common::dlist& output);
  void createGamma(uint n, double entry, Common::dlist& output);

  template <typename T>
  void generate(vec2dlist& r, vec2dlist& o);

  template <typename T>
  void generate(vector<T>& s);

  template <typename T>
  void fill(const vec2d& i, const vec2dlist& r,
            const Common::dlist& g, T& out);

  // Generate tiling by grid dualization method
  template <typename T>
  void dualize(const T& min, const T& max,
               const Common::dlist& gamma, GridTiling<T>& output);

  // Generate only the tiling vertices
  template <typename T>
  void dualize(const T& min, const T& max,
               const Common::dlist& gamma, vector<T>& output);

  template <typename T>
  void createVertices(GridTiling<T>& gt);

  template <typename T>
  void createLines(GridTiling<T>& gt);

  template <typename T>
  void dump2D(const GridTiling<T>& gt, vec2dlist& output);

  template <typename T>
  void dump2D(const vector<T>& verts, vec2dlist& output);

  template <typename T>
  double outerRadius(const GridTiling<T>& gt);

  template <typename T>
  void cutVertices(const GridTiling<T>& gt, vector<T>& output);

};

void GridDualizing::createGammaRandom(uint n, Common::dlist& output) {
  if (n == 0) return;

  double* values = new double[n];
  Common::random(n, values);

  output.clear();
  output.reserve(n);

  for (uint i = 0; i < n; ++i) {
    const double entry = values[i] - 0.5;

    output.push_back(entry);
  }
  delete [] values;
}

void GridDualizing::createGammaRandom(uint n, double sum, Common::dlist& output){
  if (n == 0) return;

  double* values = new double[n];
  Common::random(n, values);

  output.clear();
  output.reserve(n);

  double temp = 0.0;
  for (uint i = 0; i < n-1; ++i) {
    const double entry = values[i] - 0.5;

    temp += entry;
    output.push_back(entry);
  }
  delete [] values;

  output.push_back(sum - temp);
}

void GridDualizing::createGamma(uint n, double entry, Common::dlist& output){
  output.clear();
  output.reserve(n);
  for (uint i = 0; i < n; ++i) {
    output.push_back(entry);
  }
}

template <typename T>
void GridDualizing::generate(vec2dlist& r, vec2dlist& o) {
  using namespace Common;

  r.resize(T::size);
  o.resize(T::size);

  for (uint i = 0; i < T::size; ++i) {
    const double x = cos(2.0 * pi * double(i) / double(T::realsize));
    const double y = sin(2.0 * pi * double(i) / double(T::realsize));

    r[i].set(x, y);
    o[i].set(-y, x);
  }
}

template <typename T>
void GridDualizing::generate(vector<T>& s) {
  s.resize(T::size);

  for (uint i = 0; i < T::size; ++i) {
    for (uint j = 0; j < T::size; ++j) {
      s[i][j] = (j == i) ? 1 : 0;
    }
  }
}

template <typename T>
void GridDualizing::fill(const vec2d& i, const vec2dlist& r,
               const Common::dlist& g, T& out) {
  for (uint k = 0; k < T::size; ++k) {
     out[k] = ceil(i.dot(r[k]) + g[k]); // K(j) = ceil(Re(z \times xi^{-j}) + gamma_j)
  }
}

template <typename T>
void GridDualizing::dualize(const T& min, const T& max,
               const Common::dlist& gamma, GridTiling<T>& output) {
  if (T::size == 0) return;
  if (gamma.size() != T::size) return;

  // generate roots and orthogonal roots
  vec2dlist roots, rootsorth;
  generate<T>(roots, rootsorth);

  vector<T> step;
  generate(step);

  output.tiles.clear();
  // Reserving memory here is complicated since the amount of produced tiles
  // is kinda unpredictable (depends on mix, max, gamma and of course on T itself)

  T tij;

  cerr << "Regularity check not yet available!" << endl
       << "Note: No check on regularity of grid will be performed." << endl
       << "      Result is unpredictable if there are points where" << endl
       << "      more than two lines intersect!" << endl;

  for (uint i = 0; i < T::size-1; ++i) {
    const double gi = gamma[i]; // entry from translation vector (for i)
    const vec2d vi(roots[i]); // root entry (for i)
    const vec2d vior(rootsorth[i]); // orthogonal root entry

    const int ri[2] = {min[i], max[i]}; // range for i

    for (uint j = i+1; j < T::size; ++j) {
      const double gj = gamma[j]; // entry from translation vector (for j)
      const vec2d vj(roots[j]); // root entry (for j)

      const int rj[2] = {min[j], max[j]}; // range for i

      for (int ki = ri[0]; ki <= ri[1]; ++ki) {
        for (int kj = rj[0]; kj <= rj[1]; ++kj) {
          const vec2d intersec(vi * (double(ki) - gi) +
                               vior * (((double(kj) - gj) - (double(ki) - gi) * vi.dot(vj)) / vior.dot(vj)));

          // TODO: if (checkregulariy(intersec)) {cerr << "warning: grid is not regular!" << endl;}

          fill(intersec, roots, gamma, tij);
          tij[i] = ki;
          tij[j] = kj;

          /* Grouping the vertices to a tile is only important when we're interested
           * in the adjacency information (the tile lines).
           * For just the tiling vertices the grouping isn't needed and maybe it's best
           * to already start "filtering" vertices at this point (reduction to unique
           * integer representation, unique insert, etc.)
           */
          GridTile<T> newtile(tij, tij + step[i],
                              tij + step[i] + step[j], tij + step[j]);
          newtile.clip(min, max);
          output.tiles.push_back(newtile);
        }
      }
    }
  }
}

template <typename T>
void GridDualizing::dualize(const T& min, const T& max,
             const Common::dlist& gamma, vector<T>& output) {

  // TODO: implement!
}

template <typename T>
void GridDualizing::createVertices(GridTiling<T>& gt) {
  gt.vertices.clear();

  for (typename vector<typename GridTiling<T>::ttype>::const_iterator i = gt.tiles.begin();
       i != gt.tiles.end(); ++i) {
    for (uint j = 0; j < 4; ++j) {
      const typename GridTiling<T>::vtype current(i->getVertex(j));
      if (find(gt.vertices.rbegin(), gt.vertices.rend(), current) == gt.vertices.rend()) {
        gt.vertices.push_back(current);
      }
    }
  }
}

template <typename T>
void GridDualizing::createLines(GridTiling<T>& gt) {
  gt.lines.clear();

  typename GridTiling<T>::ltype lines[4];
  for (typename vector<typename GridTiling<T>::ttype>::const_iterator i = gt.tiles.begin();
       i != gt.tiles.end(); ++i) {
    i->getLines(lines);

    for (uint j = 0; j < 4; ++j) {
      typename GridTiling<T>::ltype& current = lines[j];
      if (find(gt.lines.rbegin(), gt.lines.rend(), current) == gt.lines.rend()) {
        gt.lines.push_back(current);
      }
    }
  }
}

template <typename T>
void GridDualizing::dump2D(const GridTiling<T>& gt, vec2dlist& output) {
  if (gt.vertices.size() == 0) return;

  output.clear();
  output.reserve(gt.vertices.size());

  for (typename vector<typename GridTiling<T>::vtype>::const_iterator i = gt.vertices.begin();
       i != gt.vertices.end(); ++i) {
    output.push_back(i->to2D());
  }
}

template <typename T>
void GridDualizing::dump2D(const vector<T>& verts, vec2dlist& output) {
  if (verts.size() == 0) return;

  output.clear();
  output.reserve(verts.size());

  for (typename vector<T>::const_iterator i = verts.begin();
       i != verts.end(); ++i) {
    output.push_back(i->to2D());
  }
}

template <typename T>
double GridDualizing::outerRadius(const GridTiling<T>& gt) {
  if (gt.vertices.size() == 0) return 0.0;
  double res = 0.0;

  for (typename vector<typename GridTiling<T>::vtype>::const_iterator i = gt.vertices.begin();
       i != gt.vertices.end(); ++i) {
    const double lsq = i->to2D().lengthSquared();
    if (lsq > res) res = lsq;
  }

  return sqrt(res);
}

/* A tiling generated by grid dualization needs trimming, since vertices   *
 * at the boundary of the patch might not lie in the intersection of       *
 * all grid "stripes".                                                     */
template <typename T>
void GridDualizing::cutVertices(const GridTiling<T>& gt, vector<T>& output) {
  using namespace Common;

  if (gt.vertices.size() == 0) return;

  // The 0.8 factor was derived by experimentation with lower symmetries
  // It might be completly off for higher ones...
  output.clear();
  output.reserve(double(gt.vertices.size()) * 0.8);

  // INFO: cos(pi/n) is the ratio between outer and inner radius of a regular n-gon

  /* cutoff = r_out * cos(pi / 2n)^4 [important: 2n]                        *
   * The cutoff is based on the observation that the "useable" area of the  *
   * dualized grid is derived from the original 2n-gon boundary by nesting  *
   * rotated 2n-gons into it. This can be seen even with n=5 and rather     *
   * small grids.                                                           */
  const double cutoff = outerRadius(gt) * pow(cos(pi / double(2 * T::realsize)), 4.0);
  const double cosq = cutoff*cutoff;

  for (typename vector<typename GridTiling<T>::vtype>::const_iterator i = gt.vertices.begin();
       i != gt.vertices.end(); ++i) {

    if (i->to2D().lengthSquared() > cosq) continue;
    output.push_back(*i);
  }
}

/* TODO: Visibility test                                                     *
 * For the ray visibility test we need to implement both conjugation and     *
 * multiplication for the GridVertex class. The rationale is the following:  *
 * Let x1, x2 in Z[xi_n] (x1 neq x2 and both non-zero).                      *
 * Then x2 is occluded by x1 iff                                             *
 * there exists t > 1 : t x1 = x2 <=> z := x2 x1^{-1} in Real_{>1} <=>       *
 * Im(z) = 0 and Re(z) > 1.                                                  *
 * The first condition reduces to x2 conj(x1) = conj(x2) x1                 */
int main(int argc, char* argv[]) {
  int steps = 40;
  uint mode = 0;
  bool sector = false;

  using namespace Common;
  
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

  const uint symmetry = 11;

  // TODO: All algorithm currently only work properly for prime numbers
  assert(Meta::isprime<symmetry>());

  typedef GridDualizing::GridVertex<symmetry> grVtx;
  grVtx::initRoots();
  double gm[symmetry];
  for (uint j = 0; j < symmetry; ++j) {gm[j] = 0.5;}
  GridDualizing::GridTiling<grVtx> output;


  GridDualizing::dualize<grVtx>(grVtx(-steps), grVtx(steps),
                                dlist(gm, gm + symmetry), output);

  cerr << "stats: generated " << output.tiles.size() << " tiles" << endl;

  GridDualizing::createVertices(output);

  // the representation isn't unique, so this doesn't produce unique vertices
  cerr << "stats: tiling has " << output.vertices.size() << " non-unique vertices" << endl;

  cerr << "stats: outer radius of tiling is " << GridDualizing::outerRadius(output) << endl;

  vector<grVtx> cutverts;
  GridDualizing::cutVertices(output, cutverts);

  cerr << "stats: after cutoff " << cutverts.size() << " vertices remain" << endl;

  cerr << "stats: " << 100.0 * double(cutverts.size()) / double(output.vertices.size())
       << " percent of the vertices survived" << endl;

  vec2dlist output2d;
  GridDualizing::dump2D(output, output2d);
  //GridDualizing::dump2D(cutverts, output2d);
  cout << output2d;

  /*GridDualizing::createLines(output);

  cerr << "stats: tiling has " << output.lines.size() << " unique lines" << endl;

  cout << output;*/

  return 0;
}

