namespace ffd::bravais_propagator{

  template<int Lx, int Ly, int Lz>
  
  template<typename coord_t>
  
  [[nodiscard]] Real
  BravaisPropagator<Lx, Ly, Lz>::
  operator()(coord_t const& Xfin, coord_t const& Xini, int spin) const{
    std::size_t constexpr dim = std::tuple_size<coord_t>::value-1;


    std::array<int, 3> dr;
    dr.fill(0);
    if constexpr( dim > 0 ){
	auto x_f = std::get<1>(Xfin);
	auto x_i = std::get<1>(Xini);
	dr[0] = x_f - x_i;
	while( dr[0] < 0 ){
	  dr[0] += Lx;
	}
	while( dr[0] >= Lx ){
	  dr[0] -= Lx;
	}
      }
    
    if constexpr( dim > 1 ){
	auto x_f = std::get<2>(Xfin);
	auto x_i = std::get<2>(Xini);
	dr[1] = x_f - x_i;
	while( dr[1] < 0 ){
	  dr[1] += Ly;
	}
	while( dr[1] >= Ly ){
	  dr[1] -= Ly;
	}
      
      }

    if constexpr( dim > 2 ){
	auto x_f = std::get<3>(Xfin);
	auto x_i = std::get<3>(Xini);
	dr[2] = x_f - x_i;
	while( dr[2] < 0 ){
	  dr[2] += Lz;
	}
	while( dr[2] >= Lz ){
	  dr[2] -= Lz;
	}
      
      }


    Real dt = std::get<0>(Xfin) - std::get<0>(Xini);
    Real fermion_sign = 1.;
    while( dt < 0 ){
      fermion_sign = - fermion_sign;
      dt += beta_;
    }
    while( dt >= beta_ ){
      fermion_sign = - fermion_sign;
      dt -= beta_;
    }


    return fermion_sign*G0[spin+2*(dr[2]+Lz*(dr[1]+Ly*dr[0]))](dt);
  }

}//namespace
