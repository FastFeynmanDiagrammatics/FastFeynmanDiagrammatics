

namespace ffd::user_space{

  class SetOstreamFloatFlags{
  public:
    int precision = 10;
    SetOstreamFloatFlags() {}
    SetOstreamFloatFlags(int x): precision(x) {}
  };

  
  std::ostream& operator<<(std::ostream& out, SetOstreamFloatFlags O_){
    return out<<std::showpos<<std::scientific<<std::setprecision(O_.precision);
  }

    
}
