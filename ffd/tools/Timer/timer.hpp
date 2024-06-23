namespace ffd::user_space{
  //This is used to time the code
  //you call ini() to start the timer
  //you call elapsed() if you want to know how much time
  //has passed since the last time you called ini()
  //you call fin() to stop the timer
  //you call operator() to know how many seconds
  //have passed in total (cumulative)
  struct Timer{
    std::chrono::high_resolution_clock::time_point t_ini;
    std::chrono::high_resolution_clock::time_point t_fin;
    std::chrono::high_resolution_clock::time_point t_ela;
  
    double seconds_total = 0;

  
    inline void ini(){
      t_ini = std::chrono::high_resolution_clock::now();
    }

  
    inline double elapsed(){
      t_ela = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_ela - t_ini);
      return time_span.count();
    }
    

    inline void fin(){
      t_fin = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_fin - t_ini);
      seconds_total += time_span.count();
    }
    
    
    double operator()() const{
      return seconds_total;
    }

    
    inline void reset(){
      ini();
      seconds_total = 0.;
    }
    
    
    inline double elapsed_and_reset(){
      double const ret = elapsed();
      reset();
      return ret;
    }

  };

}//namespace ffd::timer_space

