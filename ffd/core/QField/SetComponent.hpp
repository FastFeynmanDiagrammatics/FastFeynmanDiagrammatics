namespace ffd::user_space{
  
  QField
  SetComponent(QField x,
	       int component){
    using namespace q_field;
    auto ret = x;
    int sign = component >= 0 ? 1 : -1;
    if(sign < 0) ret += (1<<(sign_bit));
    if(sign < 0) component = -component;
    ret += (component<<1);
    return ret;
  }

}//namespace
