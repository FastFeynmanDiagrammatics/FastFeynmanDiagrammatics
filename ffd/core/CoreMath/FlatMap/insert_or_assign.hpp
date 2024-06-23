namespace ffd::flat_map{

  template<typename KeyType, typename ValueType>
  void FlatMap<KeyType, ValueType>::insert_or_assign(KeyType const& key_, ValueType const& value_){
    auto it_low = std::lower_bound(std::begin(Keys), std::end(Keys), key_);
    auto it_upp = std::upper_bound(std::begin(Keys), std::end(Keys), key_);
    auto position = std::distance(std::begin(Keys), it_low);
    if(it_low == it_upp){
      Keys.resize(size(Keys)+1);
      Values.resize(size(Keys));
      for(int k=size(Keys)-1; k>position; --k){
	Keys[k] = Keys[k-1];
	Values[k] = Values[k-1];
      }
    }
    Keys[position] = key_;
    Values[position] = value_;
  }

}//namespace ffd::flat_map
