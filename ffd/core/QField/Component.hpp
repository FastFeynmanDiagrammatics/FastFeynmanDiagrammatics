namespace ffd::user_space{

  int
  Component(QField const& q){
    using namespace q_field;
    int a_cmp = (q>>1)%(1<<(c_bits-1));
    if(value_bit(q, sign_bit)) a_cmp = -a_cmp;
    return a_cmp;
  }

}//namespace
