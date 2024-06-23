

namespace ffd::flat_map{

  template<typename KeyType, typename ValueType>
  int FlatMap<KeyType, ValueType>::count(KeyType const& key_) const{
    auto it_low = std::lower_bound(std::begin(Keys), std::end(Keys), key_);
    auto it_upp = std::upper_bound(std::begin(Keys), std::end(Keys), key_);
    return it_low != it_upp;
  }

}//namespace ffd::flat_map
