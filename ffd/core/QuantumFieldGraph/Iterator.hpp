

namespace ffd::qf_graph{
  
  class QFGraphIterator{
    QuantumFieldGraph const& Graph;
    
  public:
    char Vertex=0, Dot=0, Particle=0;

    QFGraphIterator(QuantumFieldGraph const& G_): Graph(G_) {}

    QFGraphIterator(QFGraphIterator const&) = default;

    QFGraphIterator& operator=(QFGraphIterator const& it_){
      Vertex = it_.Vertex;
      Dot = it_.Dot;
      Particle = it_.Particle;
      return *this;
    }

    QFGraphIterator& operator++(){
      if(Particle == char(size(Graph[Vertex][Dot]))-1){
	if(Dot == char(size(Graph[Vertex])) - 1){
	  if(Vertex == char(size(Graph)) -1){
	    Particle++;
	  }else{
	    Vertex++;
	    Dot = 0;
	    Particle = 0;
	  }
	}else{
	  Dot++;
	  Particle = 0;
	}
      }else{
	Particle++;
      }
      return *this;
    }

    
    QFGraphIterator operator+(char j_){
      assert(( j_>=0 ));
      QFGraphIterator ret = *this;
      for(char j=0; j<j_; ++j){
	++ret;
      }
      return ret;
    }

    
    bool NotEnd(){
      return Particle < char(size(Graph[Vertex][Dot])) &&
			Dot < char(size(Graph[Vertex])) &&
			      Vertex < char(size(Graph));
    }

    // std::pair<QuantumField, std::array<char, 2>>  operator()(){
    //   return {std::ref(Graph[Vertex][Dot][Particle]), {Vertex, Dot}};
    // }
    std::pair<QuantumField, std::array<char, 2>>  operator()(){
      return {Graph[Vertex][Dot][Particle], {Vertex, Dot}};
    }
    
  };


}//namespace ffd::qf_graph
