namespace ffd::blocking_method{

  template <ullong blk_res, ullong init_blk_size, typename Field>
  
  std::pair<Field, Real>
  
  block_t<blk_res, init_blk_size, Field>::
  mean_error() const{
    using std::abs;
    using namespace ffd::user_space;
    
    
    std::array<Field, blk_res> means;
    means.fill(0.);
    std::array<ullong, blk_res> weights;
    weights.fill(0ul);
    unsigned long long weights_normalization = 0;
    for( int j: Range(blk_pos) ){
      means[j] = blk[j]/Real(blk_size);
      weights[j] = blk_size;
      weights_normalization += weights[j];
    }
    if( val_pos > 0 ){
      means[blk_pos] = blk[blk_pos]/Real(val_pos);
      weights[blk_pos] = (unsigned)val_pos;
      weights_normalization += weights[blk_pos];
    }
    

    std::array<Real, blk_res> normalized_weights;
    for( auto j: Range(blk_res) ){
      normalized_weights[j] = weights[j]*1./weights_normalization;
    }
    

    Field mean = 0.;
    for( auto j: Range(blk_res) ){
      mean += normalized_weights[j]*means[j];
    }
    

    Real error_square = 0.;
    for( int j: Range(blk_res) ){
      error_square += std::pow(normalized_weights[j]*abs(means[j]-mean), 2);
    }
    

    return std::make_pair(mean, std::sqrt(error_square));
  }

}//namespace
