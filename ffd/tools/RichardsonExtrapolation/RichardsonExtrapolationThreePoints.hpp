

namespace ffd::richardson_extrapolation{

  template<typename Field>
  std::pair<Field, Real>
  RichardsonExtrapolationThreePoints(std::vector<Field> sequence){

    
    const int size_seq = size(sequence);
    
    assert( size_seq%2 == 1 && size_seq >= 3 );

    
    for(int k=1; k<=size_seq/2; ++k){
      for(int j=0; j < size_seq - 2*k; ++j){
	Field seq_jp2 = sequence[j+2];
	Field seq_jp1 = sequence[j+1];
	Field delta_upp = seq_jp2 - seq_jp1;
	Field delta_low = seq_jp1 - sequence[j];
	
	sequence[j] = seq_jp2 -
	  delta_upp*delta_upp/(delta_upp-delta_low);
	
      }
    }

    
    return std::pair{sequence[0], std::abs(sequence[2]-sequence[0])};
  }
  
  
}//namespace
