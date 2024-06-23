namespace ffd::user_space{

  ffd::phys::statistics
  Statistics(QField x){
    using namespace q_field;
    if(!value_bit(x, statistics_bit)) return ffd::phys::fermi;
    else return ffd::phys::bose;
  }

}//namespace
