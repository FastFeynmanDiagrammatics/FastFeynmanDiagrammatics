namespace ffd::user_space{

  template<typename coord_t,
	   typename field>
  QFS<coord_t, field>
  HC(QFS<coord_t, field> const& qfs){
    QFS<coord_t, field> ret = Bar(qfs);


    if constexpr(std::is_same_v<field, Complex>){
    for(uint j=0; j<size(qfs); ++j){
      ret.coef[j] = std::conj(ret.coef[j]);
    }

      }
    return ret;
  }


  
  template<typename coord_t,
	   typename field>
  QFS<coord_t, field>
  HermitianConjugate(QFS<coord_t, field> const& qfs){
    return HC(qfs);
  }

}//namespace
