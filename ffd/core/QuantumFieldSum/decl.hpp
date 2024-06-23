namespace ffd::user_space{

  template<typename coord_t = void_t,
	   typename field = Real>
  class QuantumFieldSum{
  public:
    std::vector<std::vector<std::pair<QF, coord_t>>> fields;
    std::vector<field> coef;

    QuantumFieldSum(){}

    QuantumFieldSum
    operator()(int comp) const{
      QuantumFieldSum<coord_t, field> ret = *this;
      for(uint j=0; j<size(ret); ++j){
	for(uint k=0; k<size(ret.fields[j]); ++k){
	  ret.fields[j][k].first.component = comp;
	}
      }
      return ret;
    }
    
    
    QuantumFieldSum
    operator()(char name_component) const{
      if( name_component == 'u')
	return (*this)(1);
      if( name_component == 'd')
	return (*this)(-1);
      return (*this)(0);
    }


    template<typename coord_t2>
    typename std::enable_if_t<std::is_same_v<coord_t, void_t>,
			      QuantumFieldSum<coord_t2, field>>
    operator()(coord_t2 const&) const;

    
    std::vector<std::pair<QF, coord_t>>& operator[](int j){return fields[j];}
    std::vector<std::pair<QF, coord_t>> const& operator[](int j)const{return fields[j];}
    
  };
  
  template<typename coord_t,
	   typename field>
  using QFS = QuantumFieldSum<coord_t, field>;
  
  template<typename coord_t, typename field>
  inline std::size_t size(QFS<coord_t, field> const& x){return x.fields.size();}
  
  template<typename c_1,
	   typename c_2,
	   typename fld1,
	   typename fld2>
  typename std::enable_if_t<std::is_same_v<c_1, void_t> ||
			    std::is_same_v<c_1, c_2>,
			    QFS<c_2, decltype(std::declval<fld1>()
					      *std::declval<fld2>())>>
  operator*(QFS<c_1, fld1> const &,
	    QFS<c_2, fld2> const &);
  
  
  template<typename c_1,
	   typename fld1,
	   typename fld2>
  QFS<c_1, decltype(std::declval<fld1>()
		    *std::declval<fld2>())>
  operator*(fld2,
	    QFS<c_1, fld1> const &);


  template<typename c_1,
	   typename fld1,
	   typename fld2>
  QFS<c_1, decltype(std::declval<fld1>()
		    *std::declval<fld2>())> 
  operator*(QFS<c_1, fld1> const &,
	    fld2);


  template<typename c_1,
	   typename fld1,
	   typename fld2>
  QFS<c_1, decltype(std::declval<fld1>()
		    +std::declval<fld2>())>
  operator+(QFS<c_1, fld1> const &,
	    QFS<c_1, fld2> const &);
  
  
  template<typename c_1,
	   typename fld1,
	   typename fld2>
  QFS<c_1, decltype(std::declval<fld1>()
		    +std::declval<fld2>())>
  operator-(QFS<c_1, fld1> const &,
	    QFS<c_1, fld2> const &);


  template<typename c_1,
	   typename fld1>
  QFS<c_1, fld1>
  operator-(QFS<c_1, fld1> const &);


  
}//namespace
