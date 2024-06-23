namespace ffd::inverse_matrix_gauss::unit_test {

void UnitTest() {
  one_by_one();

  two_by_two();

  three_by_three();
  n_by_n<3>();
  n_by_n<3>();
  n_by_n<3>();

  four_by_four();

  n_by_n<5>();

  n_by_n<6>();

  n_by_n<7>();

  n_by_n<8>();

  n_by_n<9>();

  n_by_n<10>();

  n_by_n<11>();

  n_by_n<12>();

  n_by_n<13>();

  n_by_n<14>();

  //    n_by_n_zerodiag<1>();

  n_by_n_zerodiag<2>();
  n_by_n_zerodiag<3>();
  n_by_n_zerodiag<4>();
  n_by_n_zerodiag<5>();

  n_by_n_zerodiag<6>();

  n_by_n_zerodiag<6>();

  n_by_n_zerodiag<5>();

  n_by_n_zerodiag<6>();

  n_by_n_zerodiag<7>();

  n_by_n_zerodiag<8>();

  n_by_n_zerodiag<9>();

  n_by_n_zerodiag<10>();

  n_by_n_zerodiag<11>();

  n_by_n_zerodiag<12>();

  n_by_n_zerodiag<13>();

  n_by_n_zerodiag<14>();

  // n_by_n<15>();

  // n_by_n<16>();

  {
    auto [is_ok, err, err_det] = n_by_n_poly<2, 0>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }

  {
    auto [is_ok, err, err_det] = n_by_n_even_poly<2, 0>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }
  {
    auto [is_ok, err, err_det] = n_by_n_poly<2, 1>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }

  {
    auto [is_ok, err, err_det] = n_by_n_even_poly<2, 1>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }
  {
    auto [is_ok, err, err_det] = n_by_n_poly<3, 0>();
    //    std::cerr << err << " " << err_det << std::endl;
    assert((is_ok));
  }

  {
    auto [is_ok, err, err_det] = n_by_n_even_poly<3, 0>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }
  {
    auto [is_ok, err, err_det] = n_by_n_poly<3, 1>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }

  {
    auto [is_ok, err, err_det] = n_by_n_even_poly<3, 1>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }
  {
    auto [is_ok, err, err_det] = n_by_n_poly<3, 2>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }

  {
    auto [is_ok, err, err_det] = n_by_n_even_poly<3, 2>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }
  {
    auto [is_ok, err, err_det] = n_by_n_poly<3, 3>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }

  {
    auto [is_ok, err, err_det] = n_by_n_even_poly<3, 3>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }
  {
    auto [is_ok, err, err_det] = n_by_n_poly<4, 0>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }

  {
    auto [is_ok, err, err_det] = n_by_n_even_poly<4, 0>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }
  {
    auto [is_ok, err, err_det] = n_by_n_poly<4, 1>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }

  {
    auto [is_ok, err, err_det] = n_by_n_even_poly<4, 2>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }
  {
    auto [is_ok, err, err_det] = n_by_n_poly<4, 2>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }

  {
    auto [is_ok, err, err_det] = n_by_n_even_poly<4, 3>();
    // std::cerr<<err<<" "<<err_det<<std::endl;
    assert((is_ok));
  }

  // { auto [is_ok, err, err_det] = n_by_n_poly<2, 1>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<2, 2>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<3, 0>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<3, 1>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<3, 2>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<3, 3>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<3, 4>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<3, 5>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<4, 0>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<4, 1>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<4, 2>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<4, 3>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<4, 4>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<4, 5>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<5, 0>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<5, 1>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<5, 2>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<5, 3>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<5, 4>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  // { auto [is_ok, err, err_det] = n_by_n_poly<5, 5>();  std::cerr<<err<<"
  // "<<err_det<<std::endl;assert((is_ok));}

  taylor_poly<3, 3, true>();
  taylor_poly<3, 3, false>();

  taylor_poly<5, 5, true>();
  taylor_poly<5, 5, false>();
}

}  // namespace ffd::inverse_matrix_gauss::unit_test
