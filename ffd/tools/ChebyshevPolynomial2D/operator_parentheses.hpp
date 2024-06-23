namespace ffd::chebyshev_polynomial_2d{
  
  template<typename Field>
  Field ChebyshevPolynomial2D<Field>::operator()(Real x, Real y) const{
    std::array<Real, 2> X = {{x, y}};
    std::array<Real, 2> t;
    for(int j: {0, 1}){
      t[j] = (X[j] -.5*(Domain[j][1]+Domain[j][0]))*2./(Domain[j][1]-Domain[j][0]);
    }

    std::vector<Field> Coef_x_evaluated(Order);
    for(int ky=0; ky < Order; ++ky){
      std::vector<Field> d_Cheby(Order+2);
      d_Cheby[Order+1] = 0;
      d_Cheby[Order] = 0;
      for(int kx=Order-1; kx >= 1; --kx){
	d_Cheby[kx] = 2*t[0]*d_Cheby[kx+1] - d_Cheby[kx+2] + Coef[kx+ky*Order];
      }
      Coef_x_evaluated[ky] = t[0]*d_Cheby[1] - d_Cheby[2] + 0.5*Coef[0+ky*Order];
    }

    std::vector<Field> d_Cheby(Order+2);
    d_Cheby[Order+1] = 0;
    d_Cheby[Order] = 0;
    for(int ky=Order-1; ky >= 1; --ky){
      d_Cheby[ky] = 2*t[1]*d_Cheby[ky+1] - d_Cheby[ky+2] + Coef_x_evaluated[ky];
    }
    return t[1]*d_Cheby[1] - d_Cheby[2] + 0.5*Coef_x_evaluated[0];
    
  }

}//namespace
