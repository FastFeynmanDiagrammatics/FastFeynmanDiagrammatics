namespace ffd::user_space{

  template<class T>
  struct VectorContainer{
    std::vector<T> container;
    std::vector<std::size_t> split;
    VectorContainer(){split.push_back(0);}
    void push_back(std::vector<T> const& v){
      container.insert(end(container),
		       cbegin(v),
		       cend(v));
      split.push_back(size(v)+split[size(split)-1]);
    }
    std::vector<T> operator[](std::size_t j) const{
      if(size(*this)==0) return std::vector<T>();
      return std::vector<T>(begin(container)+split[j],
			    begin(container)+split[j+1]);
    }
    std::array<std::size_t, 2> start_size(std::size_t j) const{
      return {split[j], split[j+1]-split[j]};
    }
  };

  template<class T>
  std::size_t size(VectorContainer<T> const& x){return size(x.split)-1;}
  
}//namespace
