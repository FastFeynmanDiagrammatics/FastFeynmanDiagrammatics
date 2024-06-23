namespace ffd::user_space::correlated_blocks{

  template<typename coord_t,
	   typename field_t>
  
  ffd::relation_sets::RelationSets<QF>
  
  NormalCorrelatedBlocks(QuantumFieldSum<coord_t, field_t> const& qfs){
    ffd::relation_sets::RelationSets<QF> ret;
    for(uint j=0; j<size(qfs.fields); ++j){
      auto prod = qfs.fields[j];
      assert((size(prod) == 2));
      assert((prod[0].first.direction !=
	      prod[1].first.direction));
      int which_one_to_flip = 0;
      if(prod[1].first.direction ==
	 phys::ou){
	which_one_to_flip = 1;
      }
      prod[which_one_to_flip].first =
	Bar(prod[which_one_to_flip].first);
      ffd::relation_sets::
	AddRelation1ToRelationSets2({prod[0].first,
				     prod[1].first},
	  ret);
    }
    return ret;
  }


}//namespace
