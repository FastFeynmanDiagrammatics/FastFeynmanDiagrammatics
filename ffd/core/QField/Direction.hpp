namespace ffd::user_space{

  ffd::phys::direction
  Direction(QField x){
    using namespace q_field;
    if(!value_bit(x, direction_bit)) return ffd::phys::in;
    else return ffd::phys::ou;
  }

}//namespace
