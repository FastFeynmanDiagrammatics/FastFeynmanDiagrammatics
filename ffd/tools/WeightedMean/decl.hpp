namespace ffd::user_space {

template <class field>
struct WeightedMean {
  ENumber<field> mean;
  bool has_value = false;

  WeightedMean() {}
  WeightedMean(field value_, Real error_)
      : mean(ENumber<field>(value_, error_)), has_value(true) {}
  WeightedMean(std::pair<field, Real> x)
      : mean(ENumber<field>(x)), has_value(true) {}

  auto eval() const { return mean.eval(); }

  WeightedMean& operator&=(WeightedMean const& x) {
    if (!has_value) {
      *this = x;
    } else {
      auto const [val0, err0] = this->eval();
      if (std::abs(err0) < 10 * std::numeric_limits<Real>::min())
        return *this;
      auto const w0 = 1 / (err0 * err0);
      auto const [val1, err1] = x.eval();
      if (std::abs(err1) < 10 * std::numeric_limits<Real>::min())
        return *this = x;
      auto const w1 = 1 / (err1 * err1);
      auto const one_w01 = 1 / (w0 + w1);
      field const val_ret = (val0 * w0 + val1 * w1) * one_w01;
      Real const err_ret = std::sqrt(one_w01);
      mean = ENumber<field>(val_ret, err_ret);
    }
    return *this;
  }
};

}  // namespace ffd::user_space
