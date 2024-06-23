namespace ffd::user_space{

  template<typename Field>
  typename std::enable_if<std::is_same<Field, Real>::value ||
  			  std::is_same<Field, Complex>::value,
  			  QuadraticActionTerm<Field>>::type
  operator*(ffd::qf_vertex::QuantumFieldVertex const& V_,
  	    Field t_){
    QuadraticActionTerm<Field> ret;
    assert(std::size(V_) == 2);
    for(auto j: {0, 1}){
      assert(std::size(V_[j]) == 1);
      ret[j] = V_[j];
    }
    ret.Value = t_;
    return ret;
  }


  template<typename Field>
  typename std::enable_if<std::is_same<Field, Real>::value ||
  			  std::is_same<Field, Complex>::value,
			  QuadraticActionTerm<Field>>::type
  operator*(Field t_,
	    ffd::qf_vertex::QuantumFieldVertex const& V_){
    return V_*t_;
  }
	    


}//namespace ffd::quadratic_action
