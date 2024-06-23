namespace ffd::heat_bath_mc{

  template<typename vector_u_t,
  typename vector_s_t>

  std::tuple<std::vector<Real>,
  std::vector<Real>,
  std::size_t>

  heat_bath_step(vector_u_t const& v_u,
    vector_s_t const& v_s){
      std::size_t const n_u = size(v_u);
      std::size_t const n_s = size(v_s);


      std::vector<Real> t_u(n_u);
      std::vector<Real> p_u(n_u);
      std::vector<Real> t_s(n_s);


      Real C_u = 0.;
      for( std::size_t j = 0; j < n_u; ++j ){
        p_u[j] = std::abs( v_u[j] );
        C_u += p_u[j];
      }


      Real const one_C_u = 1./C_u;
      Real proba = ffd::user_space::Proba();
      Real proba_C = 0.;
      std::size_t chosen_u = 0;
      bool not_chosen = true;
      for( std::size_t j = 0; j < n_u; ++j ){
        Real const p_u_C = p_u[j]*one_C_u;
        t_u[j] = std::copysign(p_u_C, v_u[j]);
        proba_C += p_u_C;
        if(proba_C > proba && not_chosen){
          chosen_u = j;
          not_chosen = false;
        }
      }


      for( std::size_t j = 0; j < n_s; ++j ){
        t_s[j] = v_s[j]*one_C_u;
      }


      return std::make_tuple(t_u, t_s, chosen_u);
    }

  }//namespace
