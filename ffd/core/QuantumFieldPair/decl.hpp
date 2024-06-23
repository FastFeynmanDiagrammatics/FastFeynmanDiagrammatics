namespace ffd::user_space{
	
  template<typename coord_t>
  class QuantumFieldPair{
  public:
    std::array<std::pair<QF, coord_t>, 2> fields;
    QF& operator[](std::size_t j){return fields[j];}
    QF operator[](std::size_t j)const{return fields[j];}
    int sign = 1;

    
    QuantumFieldPair(){}

    
    QuantumFieldPair(std::pair<QF, coord_t> const& f1,
		     std::pair<QF, coord_t> const& f2
		     ){
      fields[0] = f1;
      fields[1] = f2;
      if(fields[0].first.direction ==
	 phys::ou){
	std::swap(fields[0], fields[1]);
	if(fields[0].first.statistics ==
	   phys::fermi){
	  sign = - sign;
	}
      }
    }


    template<typename field>
    QuantumFieldPair(QuantumFieldSum<coord_t, field> const& qfs){
      assert((size(qfs) == 1));
      assert((size(qfs.fields[0]) == 2));
      assert((qfs.fields[0][0].first.direction !=
	      qfs.fields[0][1].first.direction));
      assert((qfs.fields[0][0].first.statistics ==
	      qfs.fields[0][1].first.statistics));
      *this = QuantumFieldPair(qfs.fields[0][0],
			       qfs.fields[0][1]);
    }

  };
	
}//namespace
