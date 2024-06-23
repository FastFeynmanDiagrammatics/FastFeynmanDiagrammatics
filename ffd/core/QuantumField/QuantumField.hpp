namespace ffd::quantum_field{
  //A QuantumField is defined by a component (e.g. for spin1/2 fermions, one can use 1 and -1)
  //statistics (encoded in "IsFermion")
  //and the fact of not having a bar (Dagger==1)
  //impossibility of having a bar (Dagger==0)
  //and the fact of having a bar (Dagger==-1)
  //not being an antiparticle/hole (encoded in "NotAntiParticle")
  //We have implemented the code to be able to answer:
  //"Are two QuantumFields equal?"
  //and also "Is QuantumField A `smaller` than QuantumField B?"
  
  
  class QuantumField: public std::tuple<char, char, bool, bool>{
  public:
    inline char&       Dagger()          {return std::get<0>(*this);}
    inline char const& Dagger() const    {return std::get<0>(*this);}
    
    inline char&       Component()       {return std::get<1>(*this);}
    inline char const& Component() const {return std::get<1>(*this);}
    
    inline bool&       IsFermion()        {return std::get<2>(*this);}
    inline bool const& IsFermion() const  {return std::get<2>(*this);}
    
    inline bool&       NotNambu()         {return std::get<3>(*this);}
    inline bool const& NotNambu() const   {return std::get<3>(*this);}

    QuantumField() {std::get<3>(*this) = true;}

  };


}//namespace ffd::quantum_field
