namespace ffd::user_space{

  namespace run_tools{
    inline bool file_exists(std::string file_name) {
      std::ifstream f(file_name);
      return f.good();
    }
  } // namespace
  

  inline bool
  kill_signal_is_present(){
    return static_cast<bool>(std::ifstream("info/kill_signal"));
  }

}//namespace
