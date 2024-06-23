namespace ffd::diagmc_sampling{
  
  template<int Order, int OrderMinus,
	   typename proposer_t,
	   typename origin_t,
	   typename vector_t>

  auto
  create_new_vertices(proposer_t proposer_,
		      origin_t O_,
		      vector_t X_,
		      int vertex_to_move_){
    auto X_ret = X_;
    

    auto [vtx_mv0, vtx_mv1] = ffd::math_tools::DivMod(vertex_to_move_, Order);

    
    if constexpr(OrderMinus == 2){
	if( ffd::user_space::Proba() < .5 ){
	  std::swap(vtx_mv0, vtx_mv1);
	}


	int random_vertex;
	do{
	  random_vertex = ffd::random_distributions::RandomInRange(Order);
	}while( random_vertex == vtx_mv1 );


	auto X_start = ( (random_vertex != vtx_mv0) ? X_ret[random_vertex] : O_);
	X_ret[vtx_mv0] = proposer_(X_start);
      }
    
    
    int random_vertex = ffd::random_distributions::RandomInRange(Order);
    auto X_start = ( (random_vertex != vtx_mv1) ? X_ret[random_vertex] : O_);
    X_ret[vtx_mv1] = proposer_( X_start );


    return X_ret;
  }


}//namespace
