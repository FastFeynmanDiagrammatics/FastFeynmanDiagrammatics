namespace ffd::imaginary_time{

  template<bool b>
  inline Real Beta(ImaginaryTime<b> const& I){
    return I.LowerUpperBound[1] - I.LowerUpperBound[0];
  }

}//namespace
