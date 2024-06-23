

namespace ffd::user_space::imaginary_time_lattice_propagator{

  template<typename T>

  auto

  SwapArrayElements(std::array<T, 2> x){
    auto& [x0, x1] = x;

    
    std::swap(x0, x1);

    
    return x;
  }

  
}//namespace
