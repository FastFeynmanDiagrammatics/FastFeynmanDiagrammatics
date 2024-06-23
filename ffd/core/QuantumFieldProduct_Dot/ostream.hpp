

namespace ffd::user_space{

  using ffd::qf_dot::QuantumFieldDot;
  
  std::ostream& operator<<(std::ostream& out, const QuantumFieldDot& Q_){
    using std::size;
    if(size(Q_)==0){
      return out;
    }else if(size(Q_)==1){
      out<<Q_[0];
    }else{
      out<<"(";
      for(std::size_t j=0; j < size(Q_); ++j){
	out<<Q_[j];
	if(j != size(Q_) - 1 ){
	  out<<" ";
	}
      }
      out<<")";
    }
    return out<<"("<<Q_.Position.type().name()<<")";;
  }

}//namespace ffd::user_space
