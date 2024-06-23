namespace ffd::blocking_method{

  template <ullong blk_res, ullong init_blk_size, typename Field>
  
  void
  
  block_t<blk_res, init_blk_size, Field>::
  deferred_add(Field val)
  {
    deferred_add_acc += val;
    ++deferred_add_counter;
  }


  template <ullong blk_res, ullong init_blk_size, typename Field>
  
  void
  
  block_t<blk_res, init_blk_size, Field>::
  exec_deferred_add () {
    if ( deferred_add_counter != 0 )
      this->add(deferred_add_acc/Real(deferred_add_counter));
    deferred_add_acc = Real(0);
    deferred_add_counter = 0;
  }
  

}//namespace
