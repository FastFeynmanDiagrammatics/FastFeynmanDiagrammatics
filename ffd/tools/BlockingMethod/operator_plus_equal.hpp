namespace ffd::blocking_method{

  template <ullong blk_res, ullong init_blk_size, typename Field>

  block_t<blk_res, init_blk_size, Field>
  
  block_t<blk_res, init_blk_size, Field>::
  operator+=(Field val) noexcept{
    this->add(val);
    return *this;
  }


}//namespace
