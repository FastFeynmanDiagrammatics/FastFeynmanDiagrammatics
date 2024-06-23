namespace ffd::matsubara_frequency{
	
  template<bool is_fermion=true>
  class MatsubaraFrequency{
    public: 
     int n_freq;
     Real Beta;
     MatsubaraFrequency(Real Beta_, int n_freq_=0): n_freq(n_freq_), Beta(Beta_){};
  };

  
  template<bool is_fermion>
  Real frequency(const MatsubaraFrequency<is_fermion> & omega){
    using ffd::core_math::Pi;
    return (2*omega.n_freq + is_fermion)*Pi/omega.Beta;
  }

	
}//namespace
