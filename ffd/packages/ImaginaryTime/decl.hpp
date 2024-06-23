

namespace ffd::imaginary_time{

  template<bool NotAntiperiodic = true>
  class ImaginaryTime{
  public:
    std::array<Real, 2> LowerUpperBound;

    ImaginaryTime(Real Beta):
     LowerUpperBound{0, Beta}{}
    
  };

  
  
  using     PeriodicImaginaryTime = ImaginaryTime<true>;
  using AntiperiodicImaginaryTime = ImaginaryTime<false>;
  
  
}//namespace ffd::imaginary_time
