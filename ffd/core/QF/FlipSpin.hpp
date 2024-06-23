namespace ffd::user_space{

  QF
  FlipSpin(QF const& qf){
    auto ret = qf;
    ret.component = -ret.component;
    return ret;
  }

}//namespace
