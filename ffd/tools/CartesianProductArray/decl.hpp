namespace ffd::cartesian_product_array{

  template<typename T, int d>
  struct CartesianProductArray{
    std::array<std::vector<T>, d> Sets;
    std::array<T, d> value;
    std::array<std::size_t, d> indexes;
      
    
    CartesianProductArray(std::array<std::vector<T>, d> Sets,
			  std::array<std::size_t, d> indexes):
      Sets(Sets),
      indexes(indexes){
      for(std::size_t j=0; j<d; ++j){
	value[j] = Sets[j][indexes[j]];
      }
    }

    CartesianProductArray(std::array<std::vector<T>, d> Sets):
      Sets(Sets){
      for(std::size_t j=0; j < d; ++j){
	indexes[j] = 0;
	value[j] = Sets[j][indexes[j]];
      }
    }
    
    
  private:      
    void increment(std::size_t j){
      if( indexes[j] != std::size(Sets[j]) - 1 || j == d-1){
	++indexes[j];
      }else{
	indexes[j] = 0;
	increment(j+1);
      }
    }

  public:
    
    CartesianProductArray operator++(){
      increment(0);
      for(std::size_t j=0; j < d; ++j){
	value[j] = Sets[j][indexes[j]];
      }
      return *this;
    }

    bool operator!=(CartesianProductArray const& it) const{
      return indexes != it.indexes;
    }

    std::array<T, d>& operator*(){
      return value;
    }
    
    CartesianProductArray begin() const{
      return CartesianProductArray(Sets);
    }

    CartesianProductArray end() const{
      CartesianProductArray ret(Sets);
      for(std::size_t j=0; j < d; ++j){
	ret.indexes[j] = std::size(Sets[j]) - 1;
      }
      ret.increment(0);
      return ret;
    }

  };
	
	
}//namespace
