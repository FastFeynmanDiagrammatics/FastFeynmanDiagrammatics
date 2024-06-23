// namespace ffd::flat_map{

//   template<typename KeyType, typename ValueType>
//   ValueType const& FlatMap<KeyType, ValueType>::at(KeyType const& key_) const{
//     assert(std::lower_bound(std::begin(Keys), std::end(Keys), key_)
// 	   != std::upper_bound(std::begin(Keys), std::end(Keys), key_));
//     return Values[std::distance(std::begin(Keys),
// 				std::lower_bound(std::begin(Keys),
// 						 std::end(Keys), key_))];
//   }

// }//namespace ffd::flat_map
