namespace ffd::user_space::q_field{

  template<class int1_t, class int2_t>
  bool
  value_bit(int1_t x, int2_t j){
    int1_t one1=(1<<j);
    return ((x&one1) == one1);
  }

}//namespace
