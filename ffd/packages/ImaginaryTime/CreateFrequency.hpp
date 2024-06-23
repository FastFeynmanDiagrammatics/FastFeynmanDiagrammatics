namespace ffd::imaginary_time{

  template<bool NotAntiperiodic = true>
  auto CreateFrequency(ffd::imaginary_time::ImaginaryTime<NotAntiperiodic> const& ImagTime_, int n_freq = 0){
    return ffd::matsubara_frequency::MatsubaraFrequency<!NotAntiperiodic>(Beta(ImagTime_),n_freq);
}

}//namespace
