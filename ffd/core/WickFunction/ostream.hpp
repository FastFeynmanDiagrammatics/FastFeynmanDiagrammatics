

namespace ffd::user_space{

  std::ostream& operator<<(std::ostream& out, const ffd::wick_function::WickFunction& W_){
    using std::size;
    if(W_.Sign == -1){
      out<<"-";
    }
    out<<"< ";
    for(std::size_t k=0; k < size(W_); ++k){
      auto [Block0_1, Block1] = W_[k];
      if(size(Block1) != 0){
	for(std::size_t j=size(Block0_1)-1; ; --j){
	  out<<Block0_1[j];
	  if(j==0)
	    break;
	}
	for(std::size_t j=0; j < size(Block1); j++){
	  out<<Block1[j];
	}
      }else{
	for(std::size_t j=0; j < size(Block0_1); j++){
	  out<<Block0_1[j];
	}
      }
      if(k != size(W_)-1){
	out<<" ";
      }
    }
    return out<<" >";
  }


}//namespace ffd::user_space
