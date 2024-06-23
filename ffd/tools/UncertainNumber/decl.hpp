namespace ffd::user_space{

  struct UncertainNumber{
    Real value;
    Real error;

    UncertainNumber(std::pair<Real, Real> x): value(x.first), error(x.second) {}
    UncertainNumber(std::array<Real, 2> x): value(x[0]), error(x[1]) {}
    UncertainNumber(Real value_, Real error_): value(value_), error(error_) {}
    UncertainNumber(Real value_): value(value_), error(0.) {}

    UncertainNumber operator+(UncertainNumber const& x) const{
      return UncertainNumber(value+x.value, error+x.error);}
    UncertainNumber& operator+=(UncertainNumber const& x){return *this = *this + x;}
    UncertainNumber operator-(UncertainNumber const& x) const{
      return UncertainNumber(value-x.value, error-x.error);}
    UncertainNumber operator-=(UncertainNumber const& x){return *this = *this - x;}
    UncertainNumber operator*(Real x) const{return UncertainNumber(value*x, error*std::abs(x));}
    UncertainNumber& operator*=(Real x){return *this = *this * x;}
    UncertainNumber operator/(Real x) const{return UncertainNumber(value/x, error/std::abs(x));}
    UncertainNumber& operator/=(Real x){return *this = *this / x;}
    UncertainNumber operator*(UncertainNumber const& x) const{
      return UncertainNumber(value*x.value, error*std::abs(x.value)+std::abs(value)*x.error);}
    UncertainNumber& operator*=(UncertainNumber const& x){return *this = *this * x;}
    UncertainNumber operator/(UncertainNumber const& x) const{
      return UncertainNumber(value/x.value,
			     std::abs(error/x.value)+
			     std::abs(value*x.error/(x.value*x.value)));}
    UncertainNumber& operator/=(UncertainNumber const& x){return *this = *this / x;}

  };


  struct UncertainNumber2{
    Real value;
    Real error;

    UncertainNumber2(std::pair<Real, Real> x): value(x.first), error(x.second) {}
    UncertainNumber2(std::array<Real, 2> x): value(x[0]), error(x[1]) {}
    UncertainNumber2(Real value_, Real error_): value(value_), error(error_) {}
    UncertainNumber2(Real value_): value(value_), error(0.) {}

    UncertainNumber2 operator+(UncertainNumber2 const& x) const{
      return UncertainNumber2(value+x.value, std::hypot(error, x.error));}
    UncertainNumber2& operator+=(UncertainNumber2 const& x){return *this = *this + x;}
    UncertainNumber2 operator-(UncertainNumber2 const& x) const{
      return UncertainNumber2(value-x.value, std::hypot(error, x.error));}
    UncertainNumber2 operator-=(UncertainNumber2 const& x){return *this = *this - x;}
    UncertainNumber2 operator*(Real x) const{return UncertainNumber2(value*x, error*std::abs(x));}
    UncertainNumber2& operator*=(Real x){return *this = *this * x;}
    UncertainNumber2 operator/(Real x) const{return UncertainNumber2(value/x, error/std::abs(x));}
    UncertainNumber2& operator/=(Real x){return *this = *this / x;}
    UncertainNumber2 operator*(UncertainNumber2 const& x) const{
      return UncertainNumber2(value*x.value, std::hypot(error*x.value, value*x.error));}
    UncertainNumber2& operator*=(UncertainNumber2 const& x){return *this = *this * x;}
    UncertainNumber2 operator/(UncertainNumber2 const& x) const{
      return UncertainNumber2(value/x.value,
			      std::hypot(error/x.value,
					 value*x.error/(x.value*x.value)));}
    UncertainNumber2& operator/=(UncertainNumber2 const& x){return *this = *this / x;}

  };


}//namespace
