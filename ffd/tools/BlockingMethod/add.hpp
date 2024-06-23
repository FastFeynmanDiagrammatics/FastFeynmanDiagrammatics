namespace ffd::blocking_method{

  template <ullong blk_res, ullong init_blk_size, typename Field>
  
  void
  
  block_t<blk_res, init_blk_size, Field>::
  add(Field val)
  {
    blk[blk_pos] += val;
    Real const abs_value = std::abs(val);
    blk_var[blk_pos] += abs_value*abs_value;
    tot_pos++;
    val_pos++;
    if (val_pos == blk_size)
      {
	val_pos = 0;
	blk_pos++;
	if (blk_pos == blk_res)
	  {
	    blk_pos = blk_res/2;
	    blk_size *= 2;
	    for (uint i = 0; i < blk_res/2; i++)
	      {
		blk[i] = blk[2*i] + blk[2*i+1];
		blk_var[i] = blk_var[2*i] + blk_var[2*i+1];
	      }
	    for (uint i = blk_res/2; i < blk_res; i++)
	      {
		blk[i] = 0;
		blk_var[i] = 0;
	      }
	  }
      }
  }	


}//namespace
