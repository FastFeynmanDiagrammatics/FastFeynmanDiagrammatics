

namespace ffd::user_space{

  template<typename Field>
  QuadraticAction<Field> operator+(QuadraticAction<Field> const& A1_,
				   QuadraticAction<Field> const& A2_){
    QuadraticAction<Field> ret;
    for(auto const& a: A1_){
      ret.push_back(a);
    }
    for(auto const& a: A2_){
      ret.push_back(a);
    }
    return ret;
  }

  
  template<typename Field>
  QuadraticAction<Field> operator+(QuadraticActionTerm<Field> const& A1_,
				   QuadraticActionTerm<Field> const& A2_){
    return QuadraticAction<Field>(A1_) + QuadraticAction<Field>(A2_);
  }

  
  template<typename Field>
  QuadraticAction<Field>& QuadraticAction<Field>::operator+=(QuadraticAction<Field> const& A_){
    for(auto const& a: A_){
      (*this).push_back(a);
    }
    return *this;
  }

}//namespace ffd::user_space
  

