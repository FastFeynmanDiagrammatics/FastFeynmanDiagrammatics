namespace ffd::user_space::parameters{
  auto constexpr t_hopping = 1.;
  auto constexpr tp_hopping = -.3;
  auto constexpr U = 5.6;
  auto constexpr alpha_shift = 0.;
  auto constexpr a_sh2_BCS = 0.;
  auto constexpr a_sh_BCS = 0.;
  auto constexpr BCSdlt = 0.35;
  auto constexpr PARABdlt = 0.25;
  auto constexpr YRZdlt = 0.4;
  auto constexpr BCDdlt = 0.2;
  auto constexpr mu = 1.9;
  auto constexpr Beta = 5.;
  auto constexpr BCSBeta = 0.;
  auto constexpr BCDBeta = 1.;
  auto constexpr YRZBeta = 5.;
  auto constexpr mix_BCS = .5;
  auto constexpr mix_BCD = .25;
  auto constexpr mix_YRZ = .25;
  auto constexpr N_Chebyshev_input = 60;
  auto constexpr tR = 0.8;
  auto constexpr dmuR = 0.;
  auto constexpr LambdaR_omg_flat = 2.;
  auto constexpr n_freq_mts = 100000;
  auto constexpr h = 0.;
  std::size_t constexpr l_x = 6;
  std::size_t constexpr l_y = 6;
  std::size_t constexpr n_t = 60;


  std::size_t constexpr iter_YRZ =
    // 40;
    0;
  Real constexpr a_iter_YRZ = .7;


  std::size_t constexpr iter_BCS =
    0;
  // 0;
  Real constexpr a_iter_BCS = .7;


  std::size_t constexpr iter_BCD =
    // 40;
    0;
  Real constexpr a_iter_BCD = .7;


  std::size_t constexpr iter_mix =
    0;
  // 0;
  Real constexpr a_iter_mix = .85;


  std::size_t constexpr iter_AFB =
    // 40;
    0;
  Real constexpr a_AFB_mix = .7;


  std::size_t constexpr iter_FLAT =
    // 40;
    100;
  Real constexpr a_iter_FLAT = .7;

  
}
