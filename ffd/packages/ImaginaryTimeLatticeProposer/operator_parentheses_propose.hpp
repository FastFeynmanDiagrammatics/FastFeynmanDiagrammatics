namespace ffd::imaginary_time_lattice_proposer{

  template<typename imaginary_time_t,
	   typename lattice_t>

  typename Proposer<imaginary_time_t, lattice_t>::coordinates_t

  Proposer<imaginary_time_t, lattice_t>::
  operator()(typename Proposer<imaginary_time_t, lattice_t>::coordinates_t X1) const{
    return this->propose(X1);
  }

}//namespace
