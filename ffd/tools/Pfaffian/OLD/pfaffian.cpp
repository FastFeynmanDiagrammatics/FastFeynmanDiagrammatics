#ifndef FFD_PFAFFIAN_HEADER_NOT_BEEN_HERE
#define FFD_PFAFFIAN_HEADER_NOT_BEEN_HERE



#include<iostream>
#include<cmath>
#include<iomanip>
#include<random>



class Pfaffian{
public:
  Pfaffian(Real* A_init, int n_init){
    n = n_init;
    A = new Real[2*n*n-n];
    for(int j=0;j<2*n*n-n;j++){
      A[j] = A_init[j];
    }
  }

  void update(Real* A_init){
    for(int j=0;j<2*n*n-n;j++){
      A[j] = A_init[j];
    }
  }
  
  Real el_out(int j, int k){
    if(j==k)
      return 0;
    return A[mat_conv(j,k)]*(2*(j<k)-1);
  }
  
  void el_in(int j, int k, Real value){
    if(j==k)
      return;
    A[mat_conv(j,k)] = value*(2*(j<k)-1);
  }

  Real pfaffian(){
    Real Pf_ret = 1;
    Real* row_col_temp = new Real[2*n];

    for(int u=0;u<n;u++){
      //determination of the pivot
      Real piv = el_out(2*u, 2*u+1);
      int m_piv = 2*u;
      for(int m=2*u+2;m<2*n-1;m++){
	if(fabs(el_out(m, 2*u+1))>fabs(piv)){
	  piv = el_out(m, 2*u+1);
	  m_piv = m;
	}
      }
      if( fabs(piv) == 0 )
	return 0;
    
      //row-column exchange
      if(m_piv != 2*u){
	//we save column M_{l,2u}
	for(int l=0;l<2*n;l++){
	  row_col_temp[l] = el_out(l, 2*u);
	}
	//we write column M_{l,m_piv} into M_{l,2u}
	for(int l=0;l<2*n;l++){
	  Real value_in = el_out(l, m_piv);
	  el_in(l, 2*u, value_in);
	}
	//we write column M_{l,2u} into M_{l,m}
	for(int l=0;l<2*n;l++){
	  el_in(l, m_piv, row_col_temp[l]);
	}
	//we write the special element M_{2u,m}
	el_in(m_piv, 2*u, -row_col_temp[m_piv]);
	//we multiply the Pfaffian by -1
	Pf_ret *= -1;
      }
      Pf_ret *= piv;
      // print();
      //Gaussian elimination
      for(int k=2*u+2;k<2*n;k++){
	Real lambda_k = - el_out(k, 2*u)/el_out(2*u+1, 2*u);
	//this cycle is useless for computations. It is useful if one wants to see the matrix in the canonical form
	/*	for(int l=0;l<2*u;l++){
	  el_in(k, l, el_out(k, l) + lambda_k*el_out(2*u+1, l));
	  }*/
	for(int l=2*u+1;l<2*n;l++){
	  el_in(k, l, el_out(k, l) + lambda_k*el_out(2*u+1, l));
	}
      }
    }
    //  print();
    return Pf_ret;
  }

  
  Real pfaffian_MC(){
    std::uniform_int_distribution<int> rand_int(0, 2*n-2);
    //we need to sample pairings between 2n points (0, 1, ..., 2n-1)
    //the buffer is composed by 2n points
    //A point at position 2v of value j (i.e. buffer[2v] = j) 
    //is connected to the point at position 2v+1 of value l (i.e. buffer[2v+1]=l)
    //This pairing gives a weight M_{jl}
    //There is an additional sign coming from the sign of the permutation: in order to 
    //determine it we start from the trivial configuration buffer[j] = j, which has sign = 1
    //then, we randomly exhange nearest-neighbors points, each exchange gives a minus sign
    //The Pfaffian will be given by the mean value of such configurations,
    //multiplied by the number of pairings (2n-1)!!
    int* buffer = new int[2*n];
    for(int j=0;j<2*n;j++){
      buffer[j] = j;
    }
    long int n_MC = 3e8;
    long int n_MC_it = n_MC;
    Real Pfaff_ret = 0;
    Real sign = 1;
    while(n_MC_it--){
      //compute the weight of the configuration
      Real Pfaff_conf = 1;
      for(int v=0;v<n;v++)
	Pfaff_conf *= el_out(buffer[2*v], buffer[2*v+1]);
      Pfaff_conf *= sign;
      Pfaff_ret += Pfaff_conf;
      
      //randomly exchange two points
      int point_ex = rand_int(gen);
      sign *= -1;
      int buff_temp = buffer[point_ex];
      buffer[point_ex] = buffer[point_ex+1];
      buffer[point_ex+1] = buff_temp;
    }
    Pfaff_ret /= n_MC;
    for(int j=1;j<2*n;j+=2){
      Pfaff_ret *= j;
    }
    return Pfaff_ret;
  }
  
  void print(){
    for(int j=0;j<2*n;j++)
      cout<<"--------";
    cout<<endl;
    for(int j=0;j<2*n;j++){
      for(int k=0;k<2*n;k++){
	cout<<el_out(j,k) <<"\t";
      }
      cout<<endl;
    }
    for(int j=0;j<2*n;j++)
      cout<<"--------";
    cout<<endl;
  }

private:
  Real* A;
  int n;
  int mat_conv(int j, int k){
    if(j>k){
      int l=j;
      j=k;
      k=l;
    }
    return 2*n*j - (j*(j+1))/2 + (k-j-1);
  }
};





int main(){
  Real a=rand_real(gen), b=rand_real(gen), c=rand_real(gen), d=rand_real(gen), e=rand_real(gen), f=rand_real(gen);
  int n = 2;
  Real* A = new Real[2*n*n-n];
  A[0] = a; A[1] = b; A[2] = c; A[3]=d; A[4] = e; A[5] = f;
  Pfaffian pf1(A, n);
  cout<<"For this matrix (4x4):"<<setprecision(2)<<endl;
  pf1.print();
  Real pfaff_MC = pf1.pfaffian_MC();
  Real pfaff_comp = pf1.pfaffian();
  Real pfaff_exact = a*f-b*e+c*d;
  cout<<"Exact      ="<<setprecision(50)<<pfaff_exact<<endl;
  cout<<"Gauss      ="<<setprecision(50)<<pfaff_comp<<endl;
  cout<<"Monte Carlo="<<pfaff_MC<<endl;
  n=4;
  Real* B = new Real[2*n*n-n];
  for(int j=0;j<2*n*n-n;j++)
    B[j] = rand_real(gen);
  Pfaffian pf2(B,n);
  cout<<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<<endl;
  cout<<"For this matrix (8x8):"<<setprecision(2)<<endl;
  pf2.print();
  pfaff_MC = pf2.pfaffian_MC();
  pfaff_comp = pf2.pfaffian();
  cout<<"Gauss      ="<<setprecision(50)<<pfaff_comp<<endl;
  cout<<"Monte Carlo="<<pfaff_MC<<endl;
  n=7;
  Real* C = new Real[2*n*n-n];
  for(int j=0;j<2*n*n-n;j++)
    C[j] = rand_real(gen);
  Pfaffian pf3(C,n);
  cout<<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<<endl;
  cout<<"For this matrix (14x14):"<<setprecision(2)<<endl;
  pf3.print();
  pfaff_MC = pf3.pfaffian_MC();
  pfaff_comp = pf3.pfaffian();
  cout<<"Gauss      ="<<setprecision(50)<<pfaff_comp<<endl;
  cout<<"Monte Carlo="<<pfaff_MC<<endl;

  return 0;
}
