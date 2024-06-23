

namespace ffd::user_space{

  std::ostream& operator<<(std::ostream& out, const QuantumFieldVertex& V_){
    for(std::size_t j=0; j < size(V_); ++j){
      out<<V_[j];
      if(j != size(V_) - 1 ){
	out<<" * ";
      }
    }
    return out;
  }

}//namespace ffd::user_space
