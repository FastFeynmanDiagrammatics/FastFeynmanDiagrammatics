namespace ffd::user_space{

  template<
    typename coord_t,
    typename field>
  auto
  create_chis(QuantumFieldSum<coord_t, field> V_int){
    std::vector<QuantumFieldPair<coord_t>> ret;
    std::size_t const size_chis = size(V_int.fields[0])/2;
    
    
    for(uint j=0; j<size_chis; ++j){
      ret.push_back(fill_g0_matrix::pop_front_qfpair(V_int));
    }


    return ret;
  }


}//namespace
