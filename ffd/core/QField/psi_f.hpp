namespace ffd::user_space{

  QField
  psi_f(int component = 0){
    using namespace q_field;
    assert(( std::abs(component) < (1<<c_bits) ));
    QField ret = 0u;
    ret = SetComponent(ret, component);
    return ret;
  }
  
  QField
  psi_f(char component){
    if(component == 'u') return psi_f(1);
    if(component == 'd') return psi_f(-1);
    return psi_f(int(component));
  }

}//namespace
