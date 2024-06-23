namespace ffd::index_lattice{

  template<int d, int n_atoms = 1,
	   class sp_d = ffd::lattice_g::coord_d<d, n_atoms>,
	   class oracle_f = int>
  auto
  g(std::array<int, d> const& L,
    oracle_f const& b_oracle){
    using namespace ffd::user_space;
    auto blocks = b_oracle();
    
    std::vector<std::size_t> ab2_cumul;
    ab2_cumul.push_back(0);
    for(uint j=1; j<size(blocks.split); ++j){
      std::size_t const dc = (blocks.split[j]-blocks.split[j-1])*n_atoms;
      ab2_cumul.push_back(ab2_cumul[j-1] + dc*dc);
    }
    std::size_t tot_size = ab2_cumul[size(ab2_cumul)-1];
    for(uint j=0; j<d; ++j) tot_size *= L[j];
    
    auto Recenter = ffd::lattice_g::Recenter_g<d, n_atoms>(L);
    
    auto index_f =
      [ab2_cumul, b_oracle, blocks, L, Recenter]
      (std::pair<QField, sp_d> const& x0,
       std::pair<QField, sp_d> const& x1)
      -> std::size_t
      {
	using namespace ffd::user_space;
	auto sp01 = x0.second - x1.second;
	sp01 = Recenter(sp01);
	std::size_t index_brillouin = 0;
	std::size_t L_size = 1;
	for(int j=d-1; j>=0; --j){
	  index_brillouin += sp01[j]*L_size;
	  L_size *= L[j];
	}
	auto [block0, b_pos0] = b_oracle(x0.first);
	auto [block1, b_pos1] = b_oracle(x1.first);
	assert((block0 == block1));
	std::size_t const size_block = blocks.split[block0+1]-blocks.split[block0];
	std::size_t index_ab = ab2_cumul[block0];
	index_ab += b_pos0+size_block*n_atoms*b_pos1;
	if constexpr(n_atoms > 1){
	    index_ab += size_block*(x0.second[d]+size_block*n_atoms*x1.second[d]);
	  }
	return index_brillouin + L_size*index_ab;
      };

    return std::make_pair(index_f, tot_size);
  }

}//namespace
