namespace ffd::user_space{

  template<typename Field>
  std::ostream& operator<<(std::ostream& out_, QuadraticAction<Field> const& A_){
    using std::size;
    for(int j=0; j < size(A_); ++j){
      if(j != 0){
	out_ << " + ";
      }
      out_ << "{" << A_[j][0] << A_[j][1] << "}*(" << A_[j].Value<<")";
    }
    return out_;
  }

}//namespace ffd::quadratic_action
