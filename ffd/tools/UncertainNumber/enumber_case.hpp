namespace ffd::user_space{

  template < class Field >
  struct ENumber {
    Field value=0;
    Real error=0;

    ENumber(std::pair<Field, Real> x): value(x.first), error(x.second) {}
    ENumber(Field value_, Real error_): value(value_), error(error_) {}
    ENumber(Field value_): value(value_), error(0.) {}

    ENumber(){}

    ENumber operator+(ENumber const& x) const{
      return ENumber(value+x.value, std::hypot(error, x.error));}
    ENumber& operator+=(ENumber const& x){return *this = *this + x;}
    ENumber operator-(ENumber const& x) const{
      return ENumber(value-x.value, std::hypot(error, x.error));}
    ENumber operator-=(ENumber const& x){return *this = *this - x;}
    ENumber operator*(Field x) const{return ENumber(value*x, error*std::abs(x));}
    ENumber& operator*=(Field x){return *this = *this * x;}
    ENumber operator/(Field x) const{return ENumber(value/x, error/std::abs(x));}
    ENumber& operator/=(Field x){return *this = *this / x;}
    ENumber operator*(ENumber const& x) const{
      return ENumber(value*x.value, std::hypot(error*std::abs(x.value), std::abs(value)*x.error));}
    ENumber& operator*=(ENumber const& x){return *this = *this * x;}
    ENumber operator/(ENumber const& x) const{
      return ENumber(value/x.value,
		     std::hypot(error/std::abs(x.value),
				std::abs(value*x.error/(x.value*x.value))));}
    ENumber& operator/=(ENumber const& x){return *this = *this / x;}

    auto eval() const { return std::make_pair(value, error);}
  };



  template < class F >
  auto
  to_string(ENumber<F> const& x){
    std::stringstream ss;
    ss << x.value << "+/-" << x.error;
    return ss.str();
  }

  
}//namespace
