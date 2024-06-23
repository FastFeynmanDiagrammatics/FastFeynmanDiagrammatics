namespace ffd::blocking_method{

  using ullong = unsigned long long;

  
// this defines the block class for blocking
  template <ullong blk_res = 1024, ullong init_blk_size = 16, typename field = Real>
  
  class block_t
  {
  public:
    ullong tot_pos = 0;
    ullong val_pos = 0;
    ullong blk_size = init_blk_size;
    int blk_pos = 0;
    std::array<field, blk_res> blk;
    std::array<Real, blk_res> blk_var;
    field deferred_add_acc = Real(0);
    ullong deferred_add_counter = 0;
    
    
    block_t(){ blk.fill(0.); blk_var.fill(0.); }
    
    
    void add(field val);
    void deferred_add(field val);
    void exec_deferred_add();
    block_t operator+=(field val) noexcept;
    block_t operator&=(field val) noexcept;
    

    std::pair<field, Real> mean_error() const;
    
    std::tuple<field, Real, Real, Real> statistics() const;
    
  };

	
}//namespace
