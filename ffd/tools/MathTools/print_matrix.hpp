namespace ffd::math_tools{

  template<BinaryInt precision=6, typename matrix_t>
  void print_matrix(const matrix_t & matrix, BinaryInt rows = 0, BinaryInt cols = 0, BinaryInt offset = 0){

    std::ios oldState(nullptr);
    oldState.copyfmt(std::cout);

    if (cols==0) cols = rows;

    if (rows==0){
      rows = ffd::core_math::sqrt_int(matrix.size());
      cols = rows;
    }

    assert (offset+rows*cols <= (int)size(matrix));

    using std::cerr;
    using std::endl;
    cerr << std::scientific
    << std::setw(precision)
    << std::showpos
    << std::setprecision(precision);

    for (BinaryInt r=0; r<rows; ++r){
      for (BinaryInt c=0; c<cols; ++c){
        auto el = matrix[r*cols+c+offset];
        cerr << el << " ";
      }
      cerr << endl;
    }
    cerr << endl;
    std::cerr.copyfmt(oldState);

  }

  template<typename matrix_t>
  void print_matrix_signs(const matrix_t & matrix, BinaryInt rows = 0, BinaryInt cols = 0, BinaryInt offset = 0){

    if (cols==0) cols = rows;

    if (rows==0){
      rows = ffd::core_math::sqrt_int(matrix.size());
      cols = rows;
    }

    assert (offset+rows*cols <= (int)size(matrix));

    using std::cerr;
    using std::endl;

    for (BinaryInt r=0; r<rows; ++r){
      for (BinaryInt c=0; c<cols; ++c){
        auto el = matrix[r*cols+c+offset];
        if (el==0){
          cerr << "0 ";
        }
        else if (el>0){
          cerr << "+ ";
        }
        else if (el<0){
          cerr << "- ";
        }
        else {
          cerr << "N ";
        }
      }
      cerr << endl;
    }
    cerr << endl;

  }

}//namespace
