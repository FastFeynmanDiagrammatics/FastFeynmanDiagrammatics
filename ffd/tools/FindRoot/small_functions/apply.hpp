

namespace ffd::find_root{

  template<typename function_t,
	   typename vector_t>

  auto

  apply_vec(function_t F, vector_t vec){
    auto F_vec = vec;

    for( std::size_t j=0; j < size(vec); ++j){
      F_vec[j] = F(vec[j]);
    }

    return F_vec;
  }
  
  


}//namespace
