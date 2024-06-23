namespace ffd::user_space{
  
  struct void_t{};
  
  class QF{
  public:
    int component = 0;
    phys::statistics statistics = phys::fermi;
    phys::direction direction = phys::in;
    bool is_nambu = false;
  };

  
  bool
  operator==(QF const& q1,
	     QF const& q2);
	
}//namespace
