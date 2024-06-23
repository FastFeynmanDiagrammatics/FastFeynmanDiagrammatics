namespace ffd::eigenvalues_vectors{

  template<typename Field>
  std::pair<std::vector<Real>, std::vector<std::vector<Field>>>
  EigenvaluesVectors(std::vector<std::vector<Field>> const& M_,
		     Real precision){
    Real const Precision = precision > 0 ?
      precision :
      20*std::numeric_limits<Real>::epsilon();
    

    int const size = std::size(M_);
    auto U_product = UnitMatrix<Complex>(size);



    std::vector<std::vector<Complex>> M_complex( size, std::vector<Complex>(size, 0.) );
    for(int j = 0; j < size; ++j){
      for(int k = 0; k < size; ++k){
	M_complex[j][k] = M_[j][k];
      }
    }

    
    int NumIter = 100;
    while(  ( MatrixNorm<2>(M_complex) - TraceNorm<2>(M_complex) )
    	    / TraceNorm<2>(M_complex)     >
    	   Precision && NumIter >= 0){
      NumIter--;
      //      std::cerr<<MatrixNorm<2>(M_complex)<<" "<< TraceNorm<2>(M_complex)<<std::endl;
      for(int j = 0; j < size; ++j){
	for(int k = j+1; k < size; ++k){
	  // std::cerr<<"("<<j<<", "<<k<<")"<<std::endl;
	  // std::cerr<<M_complex[j][j]<<" "<<M_complex[j][k]<<std::endl;
	  // std::cerr<<M_complex[k][j]<<" "<<M_complex[k][k]<<std::endl;
	  auto [theta, phi] =
	    Diagonalize_2x2({{
			      {{M_complex[j][j], M_complex[j][k]}},
			      {{M_complex[k][j], M_complex[k][k]}}
		}});
	  // std::cerr<<j<<" "<<k<<" "<<theta/ffd::core_math::Pi<<" "<<phi/ffd::core_math::Pi<<std::endl;
	  // std::cerr<<MatrixNorm<2>(M_complex)<<" "<< TraceNorm<2>(M_complex)<<std::endl;
	  //  std::cerr<<phi<<std::endl;

	  
	  auto U = MultiplyMatrices(RotationMatrix<Complex>(theta, {{j, k}}, size),
				    PhaseMatrix(phi, {{j, k}}, size));

	  // auto Identity = MultiplyMatrices(U, Dagger(U));
	  // for(int j=0; j < size; ++j){
	  //   for(int k=0; k < size; ++k){
	  //     std::cerr<<Identity[j][k]<<" ";
	  //   }
	  //   std::cerr<<std::endl;
	  // }

	  // for(int j=0; j < size; ++j){
	  //   for(int k=0; k < size; ++k){
	  //     std::cerr<<M_complex[j][k]<<" ";
	  //   }
	  //   std::cerr<<std::endl;
	  // }


	  U_product = MultiplyMatrices( U, U_product );

 
	  
	  
	  auto M_complex1 = MultiplyMatrices( U, M_complex );
	  auto M_complex2 = MultiplyMatrices( M_complex1, Dagger(U) );

	  // for(int j=0; j < size; ++j){
	  //   for(int k=0; k < size; ++k){
	  //     std::cerr<<M_complex2[j][k]<<" ";
	  //   }
	  //   std::cerr<<std::endl;
	  // }


	  M_complex = M_complex2;
	  
	  Hermitize(M_complex);
	}
      }
    }

    
    std::vector<Real> Eigenvalues(size, 0.);
    for(int j = 0; j < size; ++j){
      Eigenvalues[j] = std::real(M_complex[j][j]);
    }

    
    
    std::vector<std::vector<Field>> EigenvectorMatrixField( size, std::vector<Field>(size, 0.) );
    if constexpr(std::is_same<Field, Complex>::value){
	for(int j = 0; j < size; ++j){
	  for(int k = 0; k < size; ++k){
	    EigenvectorMatrixField[j][k] = std::conj( U_product[j][k] );
	  }
	}
      }else{
      for(int j = 0; j < size; ++j){
	Real phase = std::arg( U_product[j][0] );
	for(int k = 0; k < size; ++k){
	  EigenvectorMatrixField[j][k] = std::real( U_product[j][k]*ffd::user_space::expI(-phase) );
	}
      }
    }



    return {Eigenvalues, EigenvectorMatrixField};
  }


}//namespace
