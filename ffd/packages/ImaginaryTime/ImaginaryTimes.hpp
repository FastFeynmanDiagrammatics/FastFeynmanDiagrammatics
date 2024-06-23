
namespace ffd::user_space{

auto ImaginaryTimes(Real Beta){
   return std::make_pair(ffd::imaginary_time::PeriodicImaginaryTime(Beta), ffd::imaginary_time::AntiperiodicImaginaryTime(Beta));
}

}//namespace
