

namespace ffd::find_root{
  
  template<typename function_t>
  
  std::optional<Real>
  
  FindRealRoot(function_t F,
	       Real x0,
	       Real AbsolutePrecision,
	       std::array<std::optional<Real>, 2> domain){

    for(int j: {0, 1} ){
      if( domain[j].has_value() ){
	if( (domain[j].value()-x0)*(1-2*j) > 0){
	  assert(( false ));
	  std::optional<Real> void_ret;
	  return void_ret;
	}
      }
    }

    
    if( domain[0].has_value() && domain[1].has_value() ){
      auto [xL, xR] = domain;
      assert(( xL.value() < x0 && x0 < xR.value()  ));
      if( F(xL.value()) * F(x0) > 0 &&
	  F(xR.value()) * F(x0) > 0 ){
	return {};
      }
    }

    
    Real delta_x = parameters::delta_x0;
    std::array<Real, 2> xLR, FLR;

    
    unsigned long iter = 0;
    for(; iter < parameters::num_iter_max; ++iter){
      int LR = (delta_x > 0);
      if( domain[LR].has_value() ){
	if( (domain[LR].value() -(x0+delta_x))*(1-2*LR) > 0){
	  delta_x = domain[LR].value() - x0;
	}
      }
      xLR = order_bracket({x0, x0+delta_x});
      FLR = apply_vec(F, xLR);

      
      if( FLR[0]*FLR[1] > 0 ){
	if(
	   (  std::abs( FLR[0] ) - std::abs( FLR[1] )  )*delta_x < 0
	   ){
	  delta_x /= -parameters::factor_bracket_search;
	}else{
	  delta_x *= parameters::factor_bracket_search;
	}
      }else{
	break;
      }
    }
    if(iter == parameters::num_iter_max){
      return {};
    }else{
      iter = 0;
    }


    return FindRealRoot(F, xLR, AbsolutePrecision);
    
  }
    



  template<typename function_t>
  
  std::optional<Real>

  FindRealRoot(function_t F,
	       std::array<Real, 2> xLR,
	       Real AbsolutePrecision){
    std::array<Real, 2> FLR;
    for( int j: {0, 1} ){
      FLR[j] = F(xLR[j]);
    }
    

    unsigned long iter = 0ul;
    Real Fm, xm;
    do{
      xm = .5*(xLR[0] + xLR[1]);
      Fm = F(xm);
      if( FLR[0]*Fm < 0){
	xLR[1] = xm;
	FLR[1] = Fm;
      }else{
	xLR[0] = xm;
	FLR[0] = Fm;
      }
      ++iter;
    }while( std::abs(Fm) > AbsolutePrecision &&
	    iter < parameters::num_iter_max);


    if(iter == parameters::num_iter_max){
      std::optional<Real> void_ret;
      return void_ret;
    }else{
      return xm;
    }
  }

}//namespace
