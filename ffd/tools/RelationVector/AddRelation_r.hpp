namespace ffd::relation_vector{

  template<class T>
  void
  AddRelation_r(std::array<T, 2> const& r,
		Container<T>& C){
    std::array<int, 2> sets;
    for(int s: {0, 1}){
      auto [set, discard] = WhichSet(r[s], C);
      sets[s] = set;
    }
    if(r[0] == r[1]){
      if(sets[0] == -1){
	C.push_back(std::vector<T>{r[0]});
      }
    }else if(sets[0] == -1 && sets[1] == -1){
      C.push_back(std::vector<T>(begin(r), end(r)));
    }else if(sets[0] == -1){
      std::size_t const pos =
	std::distance(begin(C.container),
		      std::find(begin(C.container),
				end(C.container),
				r[1]));
      C.container.resize(size(C.container)+1);
      for(uint j=size(C.container)-1; j>pos; --j){
	C.container[j] = C.container[j-1];
      }
      C.container[pos] = r[0];
      for(uint j=sets[1]+1; j<size(C.split); ++j){
	C.split[j] += 1;
      }
    }else if(sets[1] == -1){
      AddRelation_r({r[1], r[0]}, C);
    }else if(sets[0] < sets[1]){
      auto v0 = C[sets[0]];
      auto v1 = C[sets[1]];
      v0.insert(end(v0),
		begin(v1),
		end(v1));
      RemoveVector_r(sets[1], C);
      RemoveVector_r(sets[0], C);
      C.push_back(v0);
    }else if(sets[0] > sets[1]){
      AddRelation_r({r[1], r[0]}, C);
    }
  }

}//namespace
