namespace ffd::bubbles::unit_test{

  
  Real epsilon_k(Real kx, Real ky, Real mu0){
    return -2*(cos(kx)+cos(ky)) - mu0;
  }

  
  auto
  G0_square_lattice(int x, int y, Real tau, Real beta, Real mu0, int L){
    Real ret = 0;

    
    using ffd::core_math::Pi;
    using namespace std;

    
    for(int kx=0; kx<L; ++kx){
      for(int ky=0; ky<L; ++ky){
	Real k_x = 2*kx*Pi/L;
	Real k_y = 2*ky*Pi/L;
	Real Ek = epsilon_k(k_x, k_y, mu0);
	ret -= cos(k_x*x)*cos(k_y*y)*exp(-tau*Ek)/(1+exp(-beta*Ek));
      }
    }
    return ret/L/L;
  }


}//namespace 
