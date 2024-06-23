namespace ffd::user_space{

  QFS<void_t, Real>
  Psi(int component = 0){
    QFS<void_t, Real> ret;
    QF qfield;
    qfield.component = component;
    qfield.statistics = phys::fermi;
    std::vector<std::pair<QF, void_t>> vector_qfield;
    vector_qfield.push_back(std::make_pair(qfield, void_t()));
    ret.fields.push_back(vector_qfield);
    ret.coef.push_back(1.);
    return ret;
  }
  
  QFS<void_t, Real>
  Psi(char name_component){
    if( name_component == 'u')
      return Psi(1);
    if( name_component == 'd')
      return Psi(-1);
    return Psi(0);
  }

}//namespace
