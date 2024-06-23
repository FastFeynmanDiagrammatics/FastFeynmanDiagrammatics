namespace ffd::user_space{

  auto
  Nambu(QF const& qf){
    auto ret = qf;
    ret.is_nambu = !ret.is_nambu;
    if(ret.component < 0){
      ret.direction =
	((ret.direction == phys::in) ?
	 phys::ou :
	 phys::in);
    }
    return ret;
  }


}//namespace
