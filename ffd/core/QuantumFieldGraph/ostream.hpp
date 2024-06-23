


namespace ffd::user_space{

  std::ostream& operator<<(std::ostream& out, const QuantumFieldGraph& G_){
    for(std::size_t j=0; j < size(G_); ++j){
      out<<G_[j];
      if(j != size(G_)-1){
	out<<"  |  ";
      }
    }
    return out;
  }

}//namespace ffd::user_space
