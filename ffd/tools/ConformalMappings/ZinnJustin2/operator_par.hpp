

namespace ffd::conformal_mappings{

  template<typename Field>
  Complex ZinnJustin2<Field>::operator()(Complex z) const{
    Complex w = (z + 2. - 2.*std::sqrt( 1. + z ) ) / z;
    Complex ret = 0., w_pow = 1.;
    for(unsigned long j=0; j < std::size(Coef_w); ++j, w_pow *= w){
      ret += w_pow*Coef_w[j];
    }
    return ret;
  }

}//namespace
