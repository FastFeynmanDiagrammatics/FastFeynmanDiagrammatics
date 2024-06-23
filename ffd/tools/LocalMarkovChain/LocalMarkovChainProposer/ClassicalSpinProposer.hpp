

namespace ffd::local_markov_chain{

  struct ClassicalSpinProposer{

    bool
    operator()(bool spin_){
      bool spin = !spin_;
      return spin;
    }

  };


}//namespace
