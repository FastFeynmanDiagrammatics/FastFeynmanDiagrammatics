namespace ffd::blocking_method{

  template <ullong blk_res, ullong init_blk_size, typename Field>
  
  std::tuple<Field, Real, Real, Real>
  
  block_t<blk_res, init_blk_size, Field>::
  statistics() const{
    using std::abs;
    using namespace ffd::user_space;
    
    
    std::array<Field, blk_res> means;
    means.fill(0.);
    std::array<Field, blk_res> variances;
    variances.fill(0.);
    std::array<ullong, blk_res> weights;
    weights.fill(0ul);
    unsigned long long weights_normalization = 0;
    for( int j: Range(blk_pos) ){
      means[j] = blk[j]/blk_size;
      variances[j] = blk_var[j]/blk_size;
      weights[j] = blk_size;
      weights_normalization += weights[j];
    }
    if( val_pos > 0 ){
      means[blk_pos] = blk[blk_pos]/val_pos;
      variances[blk_pos] = blk_var[blk_pos]/val_pos;
      weights[blk_pos] = (unsigned)val_pos;
      weights_normalization += weights[blk_pos];
    }
    

    std::array<Real, blk_res> normalized_weights;
    for( auto j: Range(blk_res) ){
      normalized_weights[j] = weights[j]*1./weights_normalization;
    }
    

    Field mean = 0.;
    Real variance = 0.;
    for( auto j: Range(blk_res) ){
      mean += normalized_weights[j]*means[j];
      variance += normalized_weights[j]*variances[j];
    }
    Real const abs_mean = std::abs(mean);
    variance -= abs_mean*abs_mean;
    

    Real error_square = 0.;
    for( int j: Range(blk_res) ){
      error_square += std::pow(normalized_weights[j]*abs(means[j]-mean), 2);
    }
    Real const autocorr_time = .5*(tot_pos*error_square/variance-1.);
    

    return std::make_tuple(mean, std::sqrt(error_square), autocorr_time, std::sqrt(variance));
  }



}//namespace
