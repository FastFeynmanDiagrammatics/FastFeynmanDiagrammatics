

namespace ffd::user_space::imaginary_time_lattice_propagator{

  auto
  
  IsDagger_set(int dagger){
    std::vector<bool> is_dagger_set;

    
    bool is_dagger_statement = (dagger == -1);

    
    is_dagger_set.push_back( is_dagger_statement );

    
    if(dagger == 0){
      is_dagger_set.push_back( !is_dagger_statement );
    }

    
    return is_dagger_set;

    
  }
  

}//namespace
