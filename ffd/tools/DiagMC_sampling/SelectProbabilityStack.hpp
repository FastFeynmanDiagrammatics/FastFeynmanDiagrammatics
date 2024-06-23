namespace ffd::diagmc_sampling{

  template<typename stack_t>
  int
  SelectProbabilityStack(stack_t stacks){

    
    Real p = ffd::user_space::Proba();

    
    auto stack_iterator = std::lower_bound( begin(stacks),
					    end(stacks),
					    p );

    
    return std::distance( begin(stacks),
			  stack_iterator );
  }


}//namespace
