template<typename prop_t>

auto
kt_propagator(prop_t const& prop_kt){
  fermi_G_t<l_x, l_y, n_t> ret;
  for(std::size_t nx=0; nx<size_x; ++nx){
    for(std::size_t ny=0; ny<size_y; ++ny){
      Real const kx = 2.*nx*Pi/size_x;
      Real const ky = 2.*ny*Pi/size_y;
      std::size_t const index = ny+size_y*nx;
      ret[index] = fermi_cheby_s<n_t>([=](Real t){return prop_kt(kx, ky, t);},
					    {0., Beta});
    }
  }

  return ret;
}
