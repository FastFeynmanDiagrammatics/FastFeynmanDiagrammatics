namespace ffd::user_space{

  template<typename... T>
  QFSum<nothing_d, Real>
  psi(T... args){
    auto psif = psi_f(args...);
    QFSum ret;
    ret.fields.push_back(std::make_pair(psif, nothing_d()));
    ret.plus_pos[1] += 1;
    return ret;
  }
  
}//namespace
