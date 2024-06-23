namespace ffd::SetT_wrapper{
  enum SetT{Full, Even, Odd};
}

namespace ffd::ranked_zeta{
  using namespace ffd::SetT_wrapper;
  
  template<typename Field, SetT SetType>
    class RankedZeta;

}


