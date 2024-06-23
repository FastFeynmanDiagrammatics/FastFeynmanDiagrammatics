namespace ffd::user_space{
  
  namespace ns_generate_thread{
    char generate_random_char(){
      int const choice = RandomInRange(62);
      if(choice < 10){
	return '0' + char(RandomInRange(10));
      }else if (choice < 10+26){
	return 'a' + char(RandomInRange(26));
      }else{ //if (choice < 10+26+26){
	return 'A' + char(RandomInRange(26));
      }
    }
  }

  std::string
  generate_thread_name(int length = 6){
    std::stringstream process_name_sstream;
    for (int i = 0; i < length; ++i) {
      process_name_sstream << ns_generate_thread::generate_random_char();
    }
    return process_name_sstream.str();
  }

}//namespace
