namespace ffd::nilpotent_polynomial{

  template<typename Field, SetT SetType>
  NilpotentPolynomial<Field, SetType>
  abs(const NilpotentPolynomial<Field, SetType>& Poly){
    using std::abs;
    NilpotentPolynomial<Field, SetType> Poly_abs(Poly.linear_size());
    for (unsigned int set=0; set < (unsigned int)Poly.size(); ++set){
      Poly_abs[set] = abs(Poly[set]);
    }
    return Poly_abs;
  }
}//namespace ffd::nilpotent_polynomial
