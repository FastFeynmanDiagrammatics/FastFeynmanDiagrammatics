namespace ffd::blocking_method{

  template <ullong blk_res, ullong init_blk_size, typename Field>
  
  ullong
  size(block_t<blk_res, init_blk_size, Field> const& acc){
    return acc.tot_pos;
  }


}//namespace
