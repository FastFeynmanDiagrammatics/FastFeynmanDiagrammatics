namespace ffd {

  
  template < class F, std::size_t nx >
  struct array1d {
    std::array<F, nx> array;
    inline F operator[](std::size_t j) const { return array[j];}
    inline F& operator[](std::size_t j) { return array[j];}
    
    template < class int_t >
    inline F operator()(int_t x) const {
      assert(( x >= 0 && x < nx ));
      return array[x];}
    template < class int_t >
    inline F& operator()(int_t x) {
      assert(( x >= 0 && x < nx ));
      return array[x];
    }

    template < class int_t >
    inline F operator()(int_t x, char mode) const {
      if ( mode == 'c' ) {
	if ( x < 0 || x >= nx ) {
	  std::cerr << "!!!!@@@@@@@@@###$$$  ERROR: ffd::array1d<" << nx
		    << ">: (" << x << ") is out of bounds!!!\n\n\n";
	}
      } else if ( mode == 'p' ) {
	while(x<0) x += nx;
	while(x>=nx) x -= nx;
      } else if ( mode == 'z' ) {
	if ( x < 0 || x >= nx ) return 0;
      }
      return array[x];
    }

    template < class int_t >
    inline F& operator()(int_t x, char mode) {
      if ( mode == 'c' ) {
	if ( x < 0 || x >= nx ) {
	  std::cerr << "!!!!@@@@@@@@@###$$$  ERROR: ffd::array1d<" << nx
		    << ">: (" << x << ") is out of bounds!!!\n\n\n";
	}
      } else if ( mode == 'p' ) {
	while(x<0) x += nx;
	while(x>=nx) x -= nx;
      }
      return array[x];
    }

    auto
    indexes(std::size_t j) const {
      std::array<std::size_t, 1> ret;
      ret[0] = j;
      return ret;
    }
    
    template<class T>
    void
    fill(T x) { array.fill(x); }

    auto
    size() const { return nx; }
  };
  
  template < class F, std::size_t nx >
  auto constexpr
  size(array1d<F, nx> const&) { return nx;}

  template < class F, std::size_t nx >
  auto constexpr
  size_x(array1d<F, nx> const&) { return nx;}
  
  template < class F, std::size_t nx >
  auto
  sizes(array1d<F, nx> const&) { return std::array<std::size_t, 1>{nx}; }


  template < class F, std::size_t nx, std::size_t ny >
  struct array2d {
    std::array<F, nx*ny> array;
    
    inline F operator[](std::size_t j) const { return array[j];}
    inline F& operator[](std::size_t j) { return array[j];}
    
    template < class int0_t, class int1_t >
    inline F operator()(int0_t x, int1_t y) const {
      assert(( x>=0 && x < nx && y >= 0 && y < ny ));
      return array[x+y*nx];}
    template < class int0_t, class int1_t >
    inline F& operator()(int0_t x, int1_t y) {
      assert(( x>=0 && x < nx && y >= 0 && y < ny ));
      return array[x+y*nx];}

    template < class int_t >
    inline F operator()(int_t x, int_t y, char mode) const {
      if ( mode == 'z' ) {
	if (x < 0 || x >= nx) return 0;
	if (y < 0 || y >= ny) return 0;
      } else if ( mode == 'p' ) {
	while(x<0) x+= nx;
	while(x>=nx) x-= nx;
	while(y<0) y+= ny;
	while(y>=ny) y-= ny;
      } else if ( mode == 'c' ) {
	if (x < 0 || x >= nx || y < 0 || y >= ny) {
	  std::cerr << "!!!!@@@@@@@@@###$$$  ERROR: ffd::array2d<" << nx << ", " << ny
		    << ">: (" << x << ", " << y << ") is out of bounds!!!\n\n\n";
	}
      } else {
	std::cerr << "\n\n!!!!@@@@@@@@@###$$$  ERROR: ffd::array2d<" << nx << " " << ny
		  << ">: mode " << mode << " not available!!\n\n";
      }
      return array[x+y*nx];
    }

    template < class int_t >
    inline F& operator()(int_t x, int_t y, char mode) {
      if ( mode == 'p' ) {
	while(x<0) x+= nx;
	while(x>=nx) x-= nx;
	while(y<0) y+= ny;
	while(y>=ny) y-= ny;
      } else if ( mode == 'c' ) {
	if (x < 0 || x >= nx || y < 0 || y >= ny) {
	  std::cerr << "\n\n!!!!@@@@@@@@@###$$$  ERROR: ffd::array2d<" << nx << " " << ny
		    << ">: (" << x << ", " << y << ") is out of bounds!!!\n\n";
	}
      }
      else {
	std::cerr << "\n\n!!!!@@@@@@@@@###$$$  ERROR: ffd::array2d<" << nx << " " << ny
		  << ">: mode not available!!\n\n";
      }
      return array[x+y*nx];
    }


    auto
    indexes(std::size_t j) const {
      std::array<std::size_t, 2> ret;
      ret[0] = j%nx;
      ret[1] = j/nx;
      return ret;
    }

    
    template<class T>
    void
    fill(T x) { array.fill(x); }

    auto constexpr
    size() const { return nx*ny;}
  };

  template < class F, std::size_t nx, std::size_t ny >
  auto constexpr
  size(array2d<F, nx, ny> const&) { return nx*ny;}

  template < class F, std::size_t nx, std::size_t ny >
  auto
  sizes(array2d<F, nx, ny> const&) { return std::array<std::size_t, 2>{nx, ny}; }

  template < class F, std::size_t nx, std::size_t ny >
  auto constexpr
  size_x(array2d<F, nx, ny> const&) { return nx; }

  template < class F, std::size_t nx, std::size_t ny >
  auto constexpr
  size_y(array2d<F, nx, ny> const&) { return ny; }

  template < class F, std::size_t nx, std::size_t ny >
  auto
  transpose(array2d<F, nx, ny> const& a) {
    array2d<F, nx, ny> ret;
    for ( ulong x = 0; x < nx; ++x ) {
      for ( ulong y = 0; y < ny; ++y ) {
	ret(y, x) = a(x, y);
      } // for y in range(0, ny)
    } // for x in range(0, nx)
    return ret;
  }


  template < class F, std::size_t nx, std::size_t ny, std::size_t nz >
  struct array3d {
    std::array<F, nx*ny*nz> array;
    inline F operator[](std::size_t j) const { return array[j];}
    inline F& operator[](std::size_t j) { return array[j];}

    template < class int0_t, class int1_t, class int2_t >
    inline F operator()(int0_t x, int1_t y, int2_t z) const {
      assert(( x>=0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz));
      return array[x+nx*(y+ny*z)];}
    template < class int0_t, class int1_t, class int2_t >
    inline F& operator()(int0_t x, int1_t y, int2_t z) {
      assert(( x>=0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz));
      return array[x+nx*(y+ny*z)];}

    template < class int0_t, class int1_t, class int2_t >
    inline F operator()(int0_t x, int1_t y, int2_t z, char mode) const {
      if ( mode == 'z' ) {
	if (x < 0 || x >= nx) return 0;
	if (y < 0 || y >= ny) return 0;
	if (z < 0 || z >= nx) return 0;
	return array[x+y*nx];
      } else if ( mode == 'p' ) {
	while(x<0) x+= nx;
	while(x>=nx) x-= nx;
	while(y<0) y+= ny;
	while(y>=ny) y-= ny;
	while(z<0) z+= nz;
	while(z>=nz) z-= nz;
	return array[x+nx*(y+ny*z)];
      } else if ( mode == 'c' ) {
	if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz) {
	  std::cerr << "!!!!@@@@@@@@@###$$$  ERROR: ffd::array3d<" << nx << ", " << ny << ", " << nz
		    << ">: (" << x << ", " << y << ", " << z << ") is out of bounds!!!\n\n\n";
	  return 0;
	}
	return array[x+nx*(y+ny*z)];
      } else {
	return 0;
      }
    }
    
    template < class int0_t, class int1_t, class int2_t >
    inline F& operator()(int0_t x, int1_t y, int2_t z, char mode) {
      if ( mode == 'z' ) {
	return array[x+y*nx];
      } else if ( mode == 'p' ) {
	while(x<0) x+= nx;
	while(x>=nx) x-= nx;
	while(y<0) y+= ny;
	while(y>=ny) y-= ny;
	while(z<0) z+= nz;
	while(z>=nz) z-= nz;
	return array[x+nx*(y+ny*z)];
      } else if ( mode == 'c' ) {
	if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz) {
	  std::cerr << "!!!!@@@@@@@@@###$$$  ERROR: ffd::array3d<" << nx << ", " << ny << ", " << nz
		    << ">: (" << x << ", " << y << ", " << z << ") is out of bounds!!!\n\n\n";
	  return array[x+nx*(y+ny*z)];
	}
	return array[x+nx*(y+ny*z)];
      } else {
	return array[x+nx*(y+ny*z)];
      }
    }

    auto
    indexes(std::size_t j) const {
      std::array<std::size_t, 3> ret;
      ret[0] = j%nx;
      ret[1] = (j/nx)%ny;
      ret[2] = j/(nx*ny);
      return ret;
    }

    template<class T>
    void
    fill(T x) { array.fill(x); }

    auto
    size() { return size(array);}
  };
  
  template < class F, std::size_t nx, std::size_t ny, std::size_t nz >
  constexpr auto
  size(array3d<F, nx, ny, nz> const&) { return nx*ny*nz;}

  template < class F, std::size_t nx, std::size_t ny, std::size_t nz >
  auto
  sizes(array3d<F, nx, ny, nz> const&) { return std::array<std::size_t, 3>{nx, ny, nz}; }


}//namespace
