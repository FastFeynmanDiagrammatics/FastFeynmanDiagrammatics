#include <iostream>
#include <sstream>
#include <array>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <time.h>
#include <chrono>
#include <random>
#include <string>
#include <stdio.h>
#include"../../../../std.hpp"
#include"../../../CoreMath.hpp"
#include"../../../../tools/RandomDistributions.hpp"
#include"../../../NilpotentPolynomial.hpp"
#include"../../../../tools/Timer.hpp"


#define N 10
#define set_N (1 << N)
#define set_N_1 (set_N-1)

using namespace std;
using namespace std::chrono;
using namespace ffd::user_space;

std::random_device rd;
std::mt19937 rng(rd());
std::uniform_real_distribution<double> dist(0.5,1);
int cnt_total = 1000;


using NilPoly = ffd::nilpotent_polynomial::NilpotentPolynomialF<N>;


//---------------------------------- POLY ------------------------------------//

inline int popcnt(unsigned int i)
{
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

class poly_t
{
public:
  double elmt[set_N] = {0};
  int mask = set_N_1;
  bool emptyset;
  void clear();
  void empty();
};

void poly_t::clear()
{
  for (unsigned int set=0; set<set_N; ++set)
  {
    elmt[set] = 0;
  }
  mask = 0;
  emptyset = false;
};

void poly_t::empty()
{
  elmt[0] = 0;
  mask = 0;
  emptyset = true;
};

inline void poly_scalar(poly_t &A, double V)
{
  unsigned int set_K = A.mask;
  for (unsigned int set = set_K; set != 0; set = ((set-1) & set_K))
  {
    A.elmt[set]*=V;
  }
  A.elmt[0]*=V;
}

inline void poly_shift(poly_t &A, poly_t &B, unsigned int set_sh)
{
  if (A.emptyset)
  {
    B=A;
    return;
  }
  B.clear();
  unsigned int set_K = A.mask - (A.mask & set_sh);
  for (unsigned int set = set_K; set != 0; set = ((set-1) & set_K))
  {
    B.elmt[set+set_sh] = A.elmt[set];
  }
  B.elmt[set_sh] = A.elmt[0];
  B.mask = A.mask | set_sh;
}

inline void poly_neg_(double *A, unsigned int set_K)
{
  for (unsigned int set = set_K; set != 0; set = ((set-1) & set_K))
  {
    A[set] *= -1.0;
  }
  A[0] *= -1.0;
}

inline void poly_add_(double *A, double *B, double *C, unsigned int set_K)
{
  for (unsigned int set = set_K; set != 0; set = ((set-1) & set_K))
  {
    C[set] = A[set] + B[set];
  }
  C[0] = A[0] + B[0];
}

inline void poly_sub_(double *A, double *B, double *C, unsigned int set_K)
{
  for (unsigned int set = set_K; set != 0; set = ((set-1) & set_K))
  {
    C[set] = A[set] - B[set];
  }
  C[0] = A[0] - B[0];
}

inline void poly_mul_(double *A, double *B, double *C, int set_K)
{
  unsigned int alt_set;
  double res;
  for (unsigned int set = set_K; set != 0; set = ((set-1) & set_K))
  {
    alt_set = set_K - set;
    res = 0;
    for (unsigned int sub_set = alt_set; sub_set != 0; sub_set = ((sub_set-1) & alt_set))
    {
      res += A[sub_set] * B[alt_set-sub_set];
    }
    C[alt_set] = res + A[0] * B[alt_set];
  }
  res = 0;
  for (unsigned int sub_set = set_K; sub_set != 0; sub_set = ((sub_set-1) & set_K))
  {
    res += A[sub_set] * B[set_K-sub_set];
  }
  C[set_K] = res + A[0] * B[set_K];
}

inline void poly_div_(double *A, double *B, double *C, unsigned int set_K)
{
  unsigned int alt_set;
  double res;
  for (unsigned int set = set_K; set != 0; set = ((set-1) & set_K))
  {
    alt_set = set_K - set;
    res = 0;
    for (unsigned int sub_set = alt_set; sub_set != 0; sub_set = ((sub_set-1) & alt_set))
    {
      res += B[sub_set] * C[alt_set-sub_set];
    }
    C[alt_set] = (A[alt_set] - res) / B[0];
  }
  res = 0;
  for (unsigned int sub_set = set_K; sub_set != 0; sub_set = ((sub_set-1) & set_K))
  {
    res += B[sub_set] * C[set_K-sub_set];
  }
  C[set_K] = (A[set_K] - res) / B[0];
}

inline void poly_neg(poly_t &A)
{
  poly_neg_(A.elmt, A.mask);
}

inline void poly_add(poly_t &A, poly_t &B, poly_t &C)
{
  if (A.emptyset)
  {
    C = B;
    return;
  }
  if (B.emptyset)
  {
    C = A;
    return;
  }
  C.mask = A.mask | B.mask;
  C.emptyset=false;
  poly_add_(A.elmt, B.elmt, C.elmt, C.mask);
}

inline void poly_sub(poly_t &A, poly_t &B, poly_t &C)
{
  if (A.emptyset)
  {
    C = B;
    if (!C.emptyset)
    {
      poly_neg(C);
    }
    return;
  }
  if (B.emptyset)
  {
    C = A;
    return;
  }
  C.mask = A.mask | B.mask;
  C.emptyset=false;
  poly_sub_(A.elmt, B.elmt, C.elmt, C.mask);
}

inline void poly_mul(poly_t &A, poly_t &B, poly_t &C)
{
  if (A.emptyset || B.emptyset)
  {
    C.empty();
    return;
  }
  C.mask = A.mask | B.mask;
  C.emptyset=false;
  poly_mul_(A.elmt, B.elmt, C.elmt, C.mask);
}

inline void poly_div(poly_t &A, poly_t &B, poly_t &C)
{
  if (B.emptyset)
  {
    cout << "Error: division by empty polynomial" << endl;
  }
  else if (A.emptyset)
  {
    C.empty();
    return;
  }
  C.mask = A.mask | B.mask;
  C.emptyset=false;
  poly_div_(A.elmt, B.elmt, C.elmt, C.mask);
}


// inline void poly_mul_even_(double *A, double *B, double *C, unsigned int set_K)
// {
//   unsigned int alt_set, alt_set_2, alt_sub_set;
//   double res;
//   for (unsigned int set = set_K; set != 0; set = ((set-1) & set_K))
//   {
//     if(!__builtin_parity(set))
//     {
//       alt_set = set_K - set;
//       alt_set_2 = alt_set/2;
//       C[alt_set] = 0;
//       for (unsigned int sub_set = alt_set; sub_set > alt_set_2; sub_set = ((sub_set-1) & alt_set))
//       {
//         {
//           alt_sub_set = alt_set-sub_set;
//           C[alt_set] += A[sub_set] * B[alt_sub_set];
//           C[alt_set] += A[alt_sub_set] * B[sub_set];
//         }
//       }
//     }
//   }
//   res = 0;
//   for (unsigned int sub_set = set_K; sub_set != 0; sub_set = ((sub_set-1) & set_K))
//   {
//     res += A[sub_set] * B[set_K-sub_set];
//   }
//   C[set_K] = res + A[0] * B[set_K];
// }
//
//
// inline void poly_mul_even(poly_t &A, poly_t &B, poly_t &C)
// {
//   if (A.emptyset || B.emptyset)
//   {
//     C.empty();
//     return;
//   }
//   C.mask = A.mask | B.mask;
//   C.emptyset=false;
//   poly_mul_even_(A.elmt, B.elmt, C.elmt, C.mask);
// }

inline void poly_print(poly_t &A)
{
  int set_K = A.mask;
  // int set_K = set_N-1;
  for (int set = set_K; set > 0; set = ((set-1) & set_K))
  {
    if (A.elmt[set_K-set] != 0.0)
    {
      cout << setprecision(4) << setw(8) << A.elmt[set_K-set] << "*["<<set_K-set<<"]" << " " << endl;
    }
  }
  cout << setprecision(4) << setw(8) << A.elmt[set_K] << "*["<<set_K<<"]   mask = " << A.mask << "  empty = " << A.emptyset << endl;
}



class mob_t
{
public:
  double elmt[set_N][N+1] = {0};
  int mask = set_N_1;
};

inline void mob_print(mob_t &A)
{
  for (unsigned int set = 0; set < set_N; ++set)
  {
    cout << set << ": ";
    for (unsigned int k=0; k<=N; ++k)
    {
      cout << setw(3) << A.elmt[set][k] << "["<<k<<"] ";
    }
    cout << endl;
  }
}

inline void mob_transform(poly_t &P, mob_t &M)
{
  int alt_set, pop_set, j_1, j_set, set_K = P.mask;
  double Fhat[set_N][N+1][N+1] = {0};
  for (unsigned int set = set_K; set > 0; set = ((set-1) & set_K))
  {
    alt_set = set_K-set;
    pop_set = popcnt(alt_set);
    for (unsigned int k=0; k<=N; ++k)
    {
      if (pop_set == k)
      {
        Fhat[alt_set][k][0] = P.elmt[alt_set];
      }
      j_1 = 0;
      j_set = 1;
      for(unsigned int j=1; j<=N; ++j)
      {
        Fhat[alt_set][k][j] = Fhat[alt_set][k][j_1];
        if (j_set & alt_set)
        {
          Fhat[alt_set][k][j] += Fhat[alt_set-j_set][k][j_1];
        }
        j_1++;
        j_set = (j_set<<1);
      }
      M.elmt[alt_set][k] = Fhat[alt_set][k][N];
    }
  }
  {
    pop_set = popcnt(set_K);
    for (unsigned int k=0; k<=N; ++k)
    {
      if (pop_set == k)
      {
        Fhat[set_K][k][0] = P.elmt[set_K];
      }
      j_1 = 0;
      j_set = 1;
      for(unsigned int j=1; j<=N; ++j)
      {
        Fhat[set_K][k][j] = Fhat[set_K][k][j_1];
        if (j_set & set_K)
        {
          Fhat[set_K][k][j] += Fhat[set_K-j_set][k][j_1];
        }
        j_1++;
        j_set = (j_set<<1);
      }
      M.elmt[set_K][k] = Fhat[set_K][k][N];
    }
  }
  M.mask = set_K;
}

inline void mob_invert(mob_t &M, poly_t &P)
{
  int alt_set, pop_set, j_1, j_set, set_K = M.mask;
  double Fhat[set_N][N+1][N+1] = {0};
  for (unsigned int set = set_K; set > 0; set = ((set-1) & set_K))
  {
    alt_set = set_K-set;
    pop_set = popcnt(alt_set);
    for (unsigned int k=0; k<=N; ++k)
    {
      Fhat[alt_set][k][0] = M.elmt[alt_set][k];
      j_1 = 0;
      j_set=1;
      for(unsigned int j=1; j<=N; ++j)
      {
        Fhat[alt_set][k][j] = Fhat[alt_set][k][j_1];
        if (j_set & alt_set)
        {
          Fhat[alt_set][k][j] -= Fhat[alt_set-j_set][k][j_1];
        }
        j_1++;
        j_set = (j_set<<1);
      }
    }
    P.elmt[alt_set] = Fhat[alt_set][pop_set][N];
  }
  {
    pop_set = popcnt(set_K);
    for (unsigned int k=0; k<=N; ++k)
    {
      Fhat[set_K][k][0] = M.elmt[set_K][k];
      j_1 = 0;
      j_set=1;
      for(unsigned int j=1; j<=N; ++j)
      {
        Fhat[set_K][k][j] = Fhat[set_K][k][j_1];
        if (j_set & set_K)
        {
          Fhat[set_K][k][j] -= Fhat[set_K-j_set][k][j_1];
        }
        j_1++;
        j_set = (j_set<<1);
      }
    }
    P.elmt[set_K] = Fhat[set_K][pop_set][N];
  }
  P.mask = set_K;
}

inline void mob_add(mob_t &A, mob_t &B, mob_t &C)
{
  int set_K = A.mask | B.mask;
  for (unsigned int set = set_K; set > 0; set = ((set-1) & set_K))
  {
    for (unsigned int k=0; k<=N; ++k)
    {
      C.elmt[set][k] = A.elmt[set][k] + B.elmt[set][k];
    }
  }
  {
    for (unsigned int k=0; k<=N; ++k)
    {
      C.elmt[0][k] = A.elmt[0][k] + B.elmt[0][k];
    }
  }
  C.mask = set_K;
}

inline void mob_sub(mob_t &A, mob_t &B, mob_t &C)
{
  int set_K = A.mask | B.mask;
  for (unsigned int set = set_K; set > 0; set = ((set-1) & set_K))
  {
    for (unsigned int k=0; k<=N; ++k)
    {
      C.elmt[set][k] = A.elmt[set][k] - B.elmt[set][k];
    }
  }
  {
    for (unsigned int k=0; k<=N; ++k)
    {
      C.elmt[0][k] = A.elmt[0][k] - B.elmt[0][k];
    }
  }
  C.mask = set_K;
}

inline void mob_neg(mob_t &A, mob_t &B)
{
  int set_K = A.mask;
  for (unsigned int set = set_K; set > 0; set = ((set-1) & set_K))
  {
    for (unsigned int k=0; k<=N; ++k)
    {
      B.elmt[set][k] = - A.elmt[set][k];
    }
  }
  {
    for (unsigned int k=0; k<=N; ++k)
    {
      B.elmt[0][k] = - A.elmt[0][k];
    }
  }
  B.mask = set_K;
}

inline void mob_mul(mob_t &A, mob_t &B, mob_t &C)
{
  int set_K = A.mask | B.mask;
  for (unsigned int set = set_K; set > 0; set = ((set-1) & set_K))
  {
    for (unsigned int k=0; k<=N; ++k)
    {
      C.elmt[set][k]=0;
      for (unsigned int j=0; j<=k; ++j)
      {
        C.elmt[set][k] += A.elmt[set][j] * B.elmt[set][k-j];
      }
    }
  }
  {
    for (unsigned int k=0; k<=N; ++k)
    {
      C.elmt[0][k]=0;
      for (unsigned int j=0; j<=k; ++j)
      {
        C.elmt[0][k] += A.elmt[0][j] * B.elmt[0][k-j];
      }
    }
  }
  C.mask = set_K;
}

inline void mob_div(mob_t &A, mob_t &B, mob_t &C)
{
  int set_K = A.mask | B.mask;
  for (unsigned int set = set_K; set > 0; set = ((set-1) & set_K))
  {
    for (unsigned int k=0; k<=N; ++k)
    {
      C.elmt[set][k] = A.elmt[set][k];
      for (unsigned int j=1; j<=k; ++j)
      {
        C.elmt[set][k] -= B.elmt[set][j] * C.elmt[set][k-j];
      }
      C.elmt[set][k] /= B.elmt[0][0];
    }
  }
  {
    for (unsigned int k=0; k<=N; ++k)
    {
      C.elmt[0][k] = A.elmt[0][k];
      for (unsigned int j=1; j<=k; ++j)
      {
        C.elmt[0][k] -= B.elmt[0][j] * C.elmt[0][k-j];
      }
      C.elmt[0][k] /= B.elmt[0][0];
    }
  }
  C.mask = set_K;
}

int main(int argc, char *argv[])
{
  poly_t P1,P2;
  NilPoly P3, P4;
  
  //SPEED TEST
  for (int i=0; i<set_N; ++i)
  {
    P1.elmt[i]=dist(rng);
    P2.elmt[i]=dist(rng);
    P3[i] = dist(rng);
    P4[i] = dist(rng);
  }
  P1.emptyset=false;
  P2.emptyset=false;
  P1.mask=set_N_1;
  P2.mask=set_N_1;

  high_resolution_clock::time_point time_start, time_start1;
  high_resolution_clock::time_point time_finish, time_finish1;
  duration<double> test_time, test_time1;

  time_start = high_resolution_clock::now();
  for (int i=0; i<cnt_total; ++i)
  {
    poly_mul(P1,P2,P1);
  }
  time_finish = high_resolution_clock::now();

  cout<<P1.elmt[(1<<N)-1]<<endl;
  
  test_time = duration_cast<duration<double>>(time_finish-time_start);

  time_start1 = high_resolution_clock::now();
  for (int i=0; i<cnt_total; ++i)
  {
    P4 *= P3;
  }
  time_finish1 = high_resolution_clock::now();

  test_time1 = duration_cast<duration<double>>(time_finish1-time_start1);

  cout<<P4[(1<<N)-1]<<endl;
  
  cout << ">> testing time multiplication " << 1000000*test_time.count()/((double)cnt_total) << " \u03BCs for " << cnt_total << " configurations" <<  endl;
    cout << ">> testing time multiplication (R) " << 1000000*test_time1.count()/((double)cnt_total) << "\u03BCs for " << cnt_total << " configurations" <<  endl;
}

