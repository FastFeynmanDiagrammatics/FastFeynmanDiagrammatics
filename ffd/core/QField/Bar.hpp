namespace ffd::user_space{

  QField
  Bar(QField const& q){
    using namespace q_field;
    auto ret = q;
    if(!value_bit(ret, scalar_bit)) ret = toggle_bit(ret, direction_bit);
    return ret;
  }

}//namespace
