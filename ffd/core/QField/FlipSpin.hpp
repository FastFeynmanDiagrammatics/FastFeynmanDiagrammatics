namespace ffd::user_space{

  QField
  FlipSpin(QField const& q){
    using namespace q_field;
    auto ret = q;
    if(Component(ret) != 0 ) ret = toggle_bit(ret, sign_bit);
    return ret;
  }

}//namespace
