// #include"../../std.hpp"
namespace ffd {

template <class F>
struct vector1d {
  std::size_t nx;
  std::vector<F> vector;
  vector1d() : nx(0) {}
  vector1d(std::size_t j) : nx(j), vector(j) {}
  template <class T>
  vector1d(std::size_t j, T x) : nx(j), vector(j, x) {}

  inline F operator[](std::size_t j) const { return vector[j]; }
  inline F& operator[](std::size_t j) { return vector[j]; }

  template <class int0_t>
  inline F operator()(int0_t x) const {
    assert((x >= 0 && x < vector.size()));
    return vector[x];
  }
  template <class int0_t>
  inline F& operator()(int0_t x) {
    assert((x >= 0 && x < vector.size()));
    return vector[x];
  }

  auto indexes(std::size_t j) const {
    std::array<std::size_t, 1> ret{j};
    return ret;
  }

  template <class T>
  void fill(T x) {
    std::fill(begin(vector), end(vector), x);
  }

  auto size() const { return vector.size(); }
};

template <class F>
auto size(vector1d<F> const& x) {
  return x.vector.size();
}

template <class F>
struct vector2d {
  std::size_t nx;
  std::vector<F> vector;
  vector2d() : nx(0) {}
  vector2d(std::size_t j) : nx(j), vector(j) {}
  vector2d(std::size_t j0, std::size_t j1) : nx(j0), vector(j0 * j1) {}
  template <class T>
  vector2d(std::size_t j0, std::size_t j1, T x) : nx(j0), vector(j0 * j1, x) {}

  inline F operator[](std::size_t j) const { return vector[j]; }
  inline F& operator[](std::size_t j) { return vector[j]; }

  template <class int0_t, class int1_t>
  inline F operator()(int0_t x, int1_t y) const {
    assert((x >= 0 && x < nx && y >= 0 && y < vector.size() / nx));
    return vector[x + y * nx];
  }
  template <class int0_t, class int1_t>
  inline F& operator()(int0_t x, int1_t y) {
    assert((x >= 0 && x < nx && y >= 0 && y < vector.size() / nx));
    return vector[x + y * nx];
  }

  auto indexes(std::size_t j) const {
    std::array<std::size_t, 2> ret;
    ret[0] = j % nx;
    ret[1] = j / nx;
    return ret;
  }

  template <class T>
  void fill(T x) {
    std::fill(begin(vector), end(vector), x);
  }

  auto size() const { return vector.size(); }
};

template <class F>
auto size(vector2d<F> const& x) {
  return x.vector.size();
}

template <class F, std::size_t nx, std::size_t ny>
auto sizes(vector2d<F> const& x) {
  return std::array<std::size_t, 2>(nx, x.vector.size() / nx);
}

template <class F>
auto size_x(vector2d<F> const& x) {
  return x.nx;
}

template <class F>
auto size_y(vector2d<F> const& x) {
  return x.vector.size() / x.nx;
}

template <class F>
struct vector3d {
  std::size_t nx, ny;
  std::vector<F> vector;
  vector3d() : nx(0), ny(0) {}
  vector3d(std::size_t j) : nx(j), ny(1), vector(j) {}
  vector3d(std::size_t j0, std::size_t j1) : nx(j0), ny(j1), vector(j0 * j1) {}
  vector3d(std::size_t j0, std::size_t j1, std::size_t j2)
      : nx(j0), ny(j1), vector(j0 * j1 * j2) {}
  template <class T>
  vector3d(std::size_t j0, std::size_t j1, std::size_t j2, T x)
      : nx(j0), ny(j1), vector(j0 * j1 * j2, x) {}

  inline F operator[](std::size_t j) const { return vector[j]; }
  inline F& operator[](std::size_t j) { return vector[j]; }

  template <class int0_t, class int1_t, class int2_t>
  inline F operator()(int0_t x, int1_t y, int2_t z) const {
    assert((x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 &&
            z < vector.size() / (nx * ny)));
    return vector[x + nx * (y + ny * z)];
  }
  template <class int0_t, class int1_t, class int2_t>
  inline F& operator()(int0_t x, int1_t y, int2_t z) {
    assert((x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 &&
            z < vector.size() / (nx * ny)));
    return vector[x + nx * (y + ny * z)];
  }

  auto indexes(std::size_t j) const {
    std::array<std::size_t, 3> ret;
    ret[0] = j % nx;
    ret[1] = (j / nx) % ny;
    ret[2] = j / (nx * ny);
    return ret;
  }

  template <class T>
  void fill(T x) {
    std::fill(begin(vector), end(vector), x);
  }

  std::size_t size() const { return vector.size(); }
};

template <class F>
auto size(vector3d<F> const& x) {
  return x.vector.size();
}

template <class F>
auto sizes(vector3d<F> const& x) {
  return std::array<std::size_t, 3>(
      x.nx, x.ny,
      (x.nx != 0 && x.ny != 0) ? x.vector.size() / (x.nx * x.ny) : 0);
}

template <class F>
auto size_x(vector3d<F> const& x) {
  return x.nx;
}

template <class F>
auto size_y(vector3d<F> const& x) {
  return x.ny;
}

template <class F>
std::size_t size_z(vector3d<F> const& x) {
  return (x.nx != 0 && x.ny != 0) ? x.vector.size() / (x.nx * x.ny)
                                  : std::size_t(0);
}

template <class F>
struct vector4d {
  std::size_t n0, n1, n2;
  std::vector<F> vector;
  vector4d() : n0(0), n1(0), n2(0) {}
  vector4d(std::size_t j) : n0(j), n1(1), n2(1), vector(j) {}
  vector4d(std::size_t j0, std::size_t j1)
      : n0(j0), n1(j1), n2(1), vector(j0 * j1) {}
  vector4d(std::size_t j0, std::size_t j1, std::size_t j2)
      : n0(j0), n1(j1), n2(j2), vector(j0 * j1 * j2) {}
  vector4d(std::size_t j0, std::size_t j1, std::size_t j2, std::size_t j3)
      : n0(j0), n1(j1), n2(j2), vector(j0 * j1 * j2 * j3) {}
  template <class T>
  vector4d(std::size_t j0, std::size_t j1, std::size_t j2, std::size_t j3, T x)
      : n0(j0), n1(j1), n2(j2), vector(j0 * j1 * j2 * j3, x) {}

  inline F operator[](std::size_t j) const { return vector[j]; }
  inline F& operator[](std::size_t j) { return vector[j]; }

  template <class int0_t, class int1_t, class int2_t, class int3_t>
  inline F const& operator()(int0_t x, int1_t y, int2_t z, int3_t w) const {
    assert(!(x < 0 || x >= n0 || y < 0 || y >= n1 || z < 0 || z >= n2 ||
             w < 0 || w >= vector.size() / (n0 * n1 * n2)));
    return vector[x + n0 * (y + n1 * (z + w * n2))];
  }
  template <class int0_t, class int1_t, class int2_t, class int3_t>
  inline F& operator()(int0_t x, int1_t y, int2_t z, int3_t w) {
    assert(!(x < 0 || x >= n0 || y < 0 || y >= n1 || z < 0 || z >= n2 ||
             w < 0 || w >= vector.size() / (n0 * n1 * n2)));
    return vector[x + n0 * (y + n1 * (z + w * n2))];
  }

  auto indexes(std::size_t j) const {
    std::array<std::size_t, 4> ret;
    ret[0] = j % n0;
    ret[1] = (j / n0) % n1;
    ret[2] = (j / (n0 * n1)) % n2;
    ret[3] = j / (n0 * n1 * n2);
    return ret;
  }

  template <class T>
  void fill(T x) {
    std::fill(begin(vector), end(vector), x);
  }

  std::size_t size() const { return vector.size(); }
};

template <class F>
auto size(vector4d<F> const& x) {
  return x.vector.size();
}

template <class F>
auto sizes(vector4d<F> const& x) {
  return std::array<std::size_t, 4>(x.n0, x.n1, x.n2,
                                    (x.n0 != 0 && x.n1 != 0 && x.n2 != 0)
                                        ? x.vector.size() / (x.n0 * x.n0)
                                        : 0);
}

template <class F>
auto size_0(vector4d<F> const& x) {
  return x.n0;
}

template <class F>
auto size_1(vector4d<F> const& x) {
  return x.n1;
}

template <class F>
auto size_2(vector4d<F> const& x) {
  return x.n2;
}

template <class F>
std::size_t size_3(vector4d<F> const& x) {
  return (x.n0 != 0 && x.n1 != 0 && x.n2 != 0)
             ? x.vector.size() / (x.n0 * x.n1 * x.n2)
             : std::size_t(0);
}

template <auto n, class F>
std::size_t size_(vector4d<F> const& x) {
  static_assert(n >= 0 && n < 4);
  if constexpr (n == 0) {
    return size_0(x);
  }
  if constexpr (n == 1) {
    return size_1(x);
  }
  if constexpr (n == 2) {
    return size_2(x);
  }
  return size_3(x);
}

template <class F>
struct vector5d {
  std::size_t n0, n1, n2, n3;
  std::vector<F> vector;
  vector5d() : n0(0), n1(0), n2(0), n3(0) {}
  vector5d(std::size_t j) : n0(j), n1(1), n2(1), n3(1), vector(j) {}
  vector5d(std::size_t j0,
           std::size_t j1,
           std::size_t j2,
           std::size_t j3,
           std::size_t j4)
      : n0(j0), n1(j1), n2(j2), n3(j3), vector(j0 * j1 * j2 * j3 * j4) {}
  template <class T>
  vector5d(std::size_t j0,
           std::size_t j1,
           std::size_t j2,
           std::size_t j3,
           std::size_t j4,
           T x)
      : n0(j0), n1(j1), n2(j2), n3(j3), vector(j0 * j1 * j2 * j3 * j4, x) {}

  inline F const& operator[](std::size_t j) const { return vector[j]; }
  inline F& operator[](std::size_t j) { return vector[j]; }

  template <class Int0_t,
            class Int1_t,
            class Int2_t,
            class Int3_t,
            class Int4_t>
  inline F const& operator()(Int0_t x0,
                             Int1_t x1,
                             Int2_t x2,
                             Int3_t x3,
                             Int4_t x4) const {
    assert(!(x0 < 0 || x0 >= n0 || x1 < 0 || x1 >= n1 || x2 < 0 || x2 >= n2 ||
             x3 < 0 || x3 >= n3 || x4 < 0 ||
             x4 >= vector.size() / (n0 * n1 * n2 * n3)));
    return vector[x0 + n0 * (x1 + n1 * (x2 + n2 * (x3 + n3 * x4)))];
  }

  template <class Int0_t,
            class Int1_t,
            class Int2_t,
            class Int3_t,
            class Int4_t>
  inline F& operator()(Int0_t x0, Int1_t x1, Int2_t x2, Int3_t x3, Int4_t x4) {
    assert(!(x0 < 0 || x0 >= n0 || x1 < 0 || x1 >= n1 || x2 < 0 || x2 >= n2 ||
             x3 < 0 || x3 >= n3 || x4 < 0 ||
             x4 >= vector.size() / (n0 * n1 * n2 * n3)));
    return vector[x0 + n0 * (x1 + n1 * (x2 + n2 * (x3 + n3 * x4)))];
  }

  auto indexes(std::size_t j) const {
    std::array<std::size_t, 5> ret;
    ret[0] = j % n0;
    ret[1] = (j / n0) % n1;
    ret[2] = (j / (n0 * n1)) % n2;
    ret[3] = (j / (n0 * n1 * n2)) % n3;
    ret[4] = j / (n0 * n1 * n2 * n3);
    return ret;
  }

  template <class T>
  void fill(T x) {
    std::fill(begin(vector), end(vector), x);
  }

  std::size_t size() const { return vector.size(); }
};

template <class F>
auto size(vector5d<F> const& x) {
  return x.vector.size();
}

template <class F>
auto size_0(vector5d<F> const& x) {
  return x.n0;
}

template <class F>
auto size_1(vector5d<F> const& x) {
  return x.n1;
}

template <class F>
auto size_2(vector5d<F> const& x) {
  return x.n2;
}

template <class F>
std::size_t size_3(vector5d<F> const& x) {
  return x.n3;
}

template <class F>
std::size_t size_4(vector5d<F> const& x) {
  return (x.n0 != 0 && x.n1 != 0 && x.n2 != 0 && x.n3 != 0)
             ? x.vector.size() / (x.n0 * x.n1 * x.n2 * x.n3)
             : std::size_t(0);
}

template <auto n, class F>
std::size_t size_(vector5d<F> const& x) {
  static_assert(n >= 0 && n < 5);
  if constexpr (n == 0) {
    return size_0(x);
  }
  if constexpr (n == 1) {
    return size_1(x);
  }
  if constexpr (n == 2) {
    return size_2(x);
  }
  if constexpr (n == 3) {
    return size_3(x);
  }
  return size_4(x);
}

}  // namespace ffd
