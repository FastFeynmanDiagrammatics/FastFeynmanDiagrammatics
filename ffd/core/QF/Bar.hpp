namespace ffd::user_space{

  QF
  Bar(QF const& qf){
    auto ret = qf;
    if(ret.direction == phys::in)
      ret.direction = phys::ou;
    else
      ret.direction = phys::in;
    return ret;
  }

}//namespace
