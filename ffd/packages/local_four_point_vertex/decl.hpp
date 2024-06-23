namespace ffd::local_four_point_vertex {

class V {
 public:
  Real const Beta;
  ffd::vector3d<Real> const grid;

  V(ffd::vector3d<Real> grid_, Real Beta_) : grid(grid_), Beta(Beta_) {
    assert((size_x(grid) == size_y(grid)));
    assert((size_x(grid) == size_z(grid)));
  }

  Real operator()(std::array<Real, 4> tau) const {
    Real sign = 1;
    std::array<Real, 3> tau_;
    for (ulong j = 0; j < 3; ++j) {
      tau_[j] = tau[j] - tau[3];
      while (tau_[j] < 0) {
        tau_[j] += Beta;
        sign = -sign;
      }
      while (tau_[j] > Beta) {
        tau_[j] -= Beta;
        sign = -sign;
      }
      assert((tau_[j] != 0));
      assert((tau_[j] != Beta));
    }  // for j in range(0, 3)
    std::array<Real, 3> dtau;
    dtau.fill(Beta / size_x(grid));
    return sign *
           ffd::linear_cubic_interpolation::Interpolate(grid, tau_, dtau);
  }
};

}  // namespace ffd::local_four_point_vertex
