

namespace ffd::local_markov_chain{

  template<typename Field>
  Real
  MixingProbabilityDistributions(std::array<Field, 2> proba_dist,
				 Real lambda){
    assert(( std::abs(lambda) <= 1. + std::numeric_limits<Real>::epsilon() ));
    Real ret = 0;
    for( int j: {0, 1} ){
      ret += .5*(1+(1-2*j)*lambda) * std::abs(proba_dist[j]);
    }
    
    return ret;
  }

}//namespace
