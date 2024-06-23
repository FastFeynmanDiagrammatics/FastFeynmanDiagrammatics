namespace ffd::user_space::q_field{

  template<class int1_t, class int2_t>
  int1_t
  toggle_bit(int1_t x, int2_t j){
    int1_t const two_j = (1<<j);
    if((x & two_j) != 0) return x - two_j;
    else return x + two_j;
  }

}//namespace
