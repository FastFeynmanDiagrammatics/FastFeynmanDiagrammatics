#pragma once

namespace ffd::diagmc_sampling::unit_test{


  void select_check(){

    std::vector<Real> stacks(2);
    stacks[0] = .5;
    stacks[1] = 1.;

    [[maybe_unused]] auto x = SelectProbabilityStack(stacks);
    
    

  }

}//namespace
