namespace ffd::user_space{

  template<typename Field>
  class QuadraticAction: public std::vector<QuadraticActionTerm<Field>>{
  public:

    QuadraticAction(){}
    
    QuadraticAction(QuadraticActionTerm<Field> const& T_){
      (*this).resize(1);
      (*this)[0] = T_;
    }

    QuadraticAction& operator+=(QuadraticAction const&);

  };

}//namespace ffd::quadratic_action
