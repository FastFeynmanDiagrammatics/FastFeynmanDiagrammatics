#ifndef FFD_NILPOTENTPOLY_LOGARITHM_HEADER
#define FFD_NILPOTENTPOLY_LOGARITHM_HEADER

namespace ffd::nilpotent_polynomial{

  template<int degree, typename Field>
  NilpotentPolynomial<degree, Field> Logarithm(NilpotentPolynomial<degree, Field> argument){
    NilpotentPolynomial<degree, Field> ret;
    if(argument.ClassNotInitialized)
      argument.InitializeClass();
    ret.set = argument.set;
    if(ret.set){
      ret[0] = std::log(argument[0]);
      argument /= argument[0];
      for(long V=1; V<ret.NumberCoef; ++V){
	ret[V] = argument[V]*ffd::set_theory::cardinality_set(V);
      }
      for(long V=0; V<ret.NumberSubsets - ret.NumberCoef; ++V){
	ret[ret.Sets[V]] -= ret[ret.Subsets[V]]*argument[ret.Sets[V]-ret.Subsets[V]];
      }
      for(long V=1; V<ret.NumberCoef; ++V){
	ret[V] /= ffd::set_theory::cardinality_set(V);
      }
    }
    return ret;
  }

}  

#endif
