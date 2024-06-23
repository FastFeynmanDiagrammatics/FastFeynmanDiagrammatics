namespace ffd::user_space{

  using QField = uint_fast16_t;

  namespace q_field{
    inline constexpr uint sign_bit = 0;
    inline constexpr uint c_bits = 8;
    inline constexpr uint scalar_bit = c_bits+1;
    inline constexpr uint direction_bit = scalar_bit+1;
    inline constexpr uint statistics_bit = direction_bit+1;
    inline constexpr uint nambu_bit = statistics_bit+1;
  }
  
}//namespace
