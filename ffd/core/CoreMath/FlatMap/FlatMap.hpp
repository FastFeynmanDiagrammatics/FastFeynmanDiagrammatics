namespace ffd::flat_map{

  //name taken from the boost container ``flat_map''
  template<typename KeyType, typename ValueType>
  class FlatMap{
  public:
    std::vector<KeyType> Keys;
    std::vector<ValueType> Values;

    void insert_or_assign(KeyType const&, ValueType const&);

    ValueType const& at(KeyType const& key_) const{
      assert(std::binary_search(std::begin(Keys), std::end(Keys), key_));
      return Values[std::distance(begin(Keys),
				  std::lower_bound(begin(Keys),
						   end(Keys), key_))];
    }
    
    int count(KeyType const&) const;
    
  };

  
  template<typename T, typename F>
  auto size(FlatMap<T, F> const& M_){
    return size(M_.Keys);
  }

}//namespace ffd::flat_map
