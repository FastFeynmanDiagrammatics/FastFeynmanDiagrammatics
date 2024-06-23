namespace ffd::user_space{

  QField
  Nambu(QField const& q){
    using namespace q_field;
    auto ret = toggle_bit(q, nambu_bit);
    if ( value_bit(q, sign_bit) ){
      ret = toggle_bit(ret, direction_bit);
    }
    return ret;
  }

}//namespace
