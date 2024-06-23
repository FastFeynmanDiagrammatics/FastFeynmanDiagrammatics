namespace ffd::alpha_f {

template <uint n,
          bool full = false,
          typename field_t = Real,
          typename PM_t = int>
auto AlphaFunction(PM_t& PM, uint shift = 0) {
  if constexpr (full) {
    const auto card = (1u << n);
    std::array<ffd::truncated_polynomial::P<n, field_t>, card> P_array;
    for (uint S = 0; S < card; ++S) {
      P_array[S][0] = PM[S + shift];
    }
    for (uint j = 1; j < card; j *= 2) {
      for (uint S = card - 1; S >= j; --S) {
        if (j & S) {
          const auto alt_S = S - j;
          const auto card_S = __builtin_popcount(S);
          for (uint k = 1; k <= card_S; ++k) {
            P_array[S][k] += P_array[alt_S][k - 1];
          }
        }
      }
    }
    return P_array;
  } else {
    const auto card = (1u << n);
    const int half_n = (n + (n & 1)) / 2;
    std::array<ffd::truncated_polynomial::P<half_n, field_t>, card> P_array;
    for (uint S = 0; S < card; ++S) {
      P_array[S][0] = PM[S + shift];
    }
    for (uint j = 1; j < card; j *= 2) {
      for (uint S = card - 1; S >= j; --S) {
        if (j & S) {
          const auto alt_S = S - j;
          for (uint k = 1; k <= std::min(__builtin_popcount(S), half_n); ++k) {
            P_array[S][k] += P_array[alt_S][k - 1];
          }
        }
      }
    }
    return P_array;
  }
}

template <uint order, uint poly_n, typename field_t = Real, typename PM_t = int>
auto AlphaFunctionSet(PM_t& PM, uint shift = 0) {
  const auto card = (1 << order);
  const auto V = card - 1 - shift;
  std::array<ffd::truncated_polynomial::P<poly_n, field_t>, card> P_array;
  for (uint S = V; S != 0; S = ((S - 1) & V)) {
    P_array[S][0] = PM[S + shift];
  }
  P_array[0][0] = PM[shift];
  for (uint j = 1; j < card; j *= 2) {
    for (uint S = V; S >= j; S = ((S - 1) & V)) {
      if (j & S) {
        const auto alt_S = S - j;
        const uint card_S = __builtin_popcount(S);
        for (uint k = 1; k <= std::min(card_S, poly_n); ++k) {
          P_array[S][k] += P_array[alt_S][k - 1];
        }
      }
    }
  }
  return P_array;
}

template <uint n,
          bool full = false,
          typename field_d = Real,
          typename arr_t = int,
          typename PM_t = int>
auto SignedAlphaFunction(PM_t& PM, const arr_t& S_array, uint shift = 0) {
  if constexpr (full) {
    const auto card = (1u << n);
    std::array<ffd::truncated_polynomial::P<n, field_d>, card> P_array;
    for (uint S = 0; S < card; ++S) {
      P_array[S][0] = PM[S + shift];
    }
    for (uint j = 0; j < n; ++j) {
      auto card_j = (1 << j);
      for (uint S = card - 1; S >= card_j; --S) {
        if (card_j & S) {
          const auto alt_S = S - card_j;
          const auto card_S = __builtin_popcount(S);
          for (uint k = 1; k <= card_S; ++k) {
            P_array[S][k] += P_array[alt_S][k - 1] * S_array[j];
          }
        }
      }
    }
    return P_array;
  } else {
    const auto card = (1u << n);
    const int half_n = (n + (n & 1)) / 2;
    std::array<ffd::truncated_polynomial::P<half_n, field_d>, card> P_array;
    for (uint S = 0; S < card; ++S) {
      P_array[S][0] = PM[S + shift];
    }
    for (uint j = 0; j < n; ++j) {
      auto card_j = (1 << j);
      for (uint S = card - 1; S >= card_j; --S) {
        if (card_j & S) {
          const auto card_S = __builtin_popcount(S);
          const auto alt_S = S - card_j;
          const auto k_max = (card_S < half_n) ? card_S : half_n;
          // #pragma GCC unroll 4
          // #pragma GCC ivdep
          for (uint k = 1; k <= k_max; ++k) {
            P_array[S][k] += P_array[alt_S][k - 1] * S_array[j];
          }
        }
      }
    }
    return P_array;
  }
}

template <uint order, uint poly_n, typename field_t = Real, typename PM_t = int>
auto AlphaFunctionSimple(PM_t& PM, uint shift = 0) {
  const auto card = (1 << order);
  const auto V = card - 1 - shift;
  std::array<field_t, card> P_array;

  for (uint S = V; S != 0; S = ((S - 1) & V)) {
    P_array[S] = PM[S + shift];
  }
  P_array[0] = PM[shift];

  return P_array;
}

}  // namespace ffd::alpha_f

// if constexpr (full){
//   const auto card = (1u<<n);
//   std::array<ffd::truncated_polynomial::P<n, full, field_d>, card> P_array;
//
//   // for (uint S=0; S<card; ++S){
//   //   auto pop_S =  __builtin_popcount(S);
//   //   P_array[S].order = pop_S;
//   //   for (uint Sp=S; Sp!=0; Sp=((Sp-1) & S)){
//   //     P_array[S][pop_S-__builtin_popcount(Sp)] += PM[Sp+shift];
//   //   }
//   //   P_array[S][pop_S] = PM[shift];
//   // }
//   // return P_array;
//
// 	for (uint S=0; S<card; ++S){
// 	  P_array[S][0] = PM[S+shift];
// 	}
// 	for (uint j=1; j<card; j*=2){
//  	  for (uint S=card-1; S>=j; --S){
// 			if (j & S){
// 			  const auto alt_S = S-j;
// 			  for (uint k=1; k<=__builtin_popcount(S); ++k){
// 				  P_array[S][k] += P_array[alt_S][k-1];
// 			  }
// 			}
// 	  }
// 	}
//   return P_array;
//
// }
// else{

// for (uint S=0; S<card; ++S){
//   int pop_S = __builtin_popcount(S);
//   P_array[S].order = std::min(pop_S, half_n);
//   for (uint Sp=S; Sp!=0; Sp=((Sp-1) & S)){
//     uint pop_Sp = pop_S - __builtin_popcount(Sp);
//     if (pop_Sp<=half_n){
//       P_array[S][pop_Sp] += PM[Sp+shift];
//     }
//   }
//   if (pop_S<=half_n){
//     P_array[S][pop_S] = PM[shift];
//   }
// }
// return P_array;
// }

//   for (uint j = 0; j < n; ++j) {
//     const auto card_j = (1 << j);
//     for (uint S_hi = 0; S_hi < (1u << (n - j - 1)); ++S_hi) {
//       const auto S_hi_part = card_j + (S_hi << (j + 1));
//       for (uint S_lo = 0; S_lo < card_j; ++S_lo) {
//         const auto S = S_lo + S_hi_part;
//         const auto alt_S = S - card_j;
//         const auto k_max =
//             (__builtin_popcount(S) < half_n) ? __builtin_popcount(S) :
//             half_n;
//         std::array<Real, half_n> P_acc;
//         P_acc.fill(S_array[j]);
// #pragma GCC unroll 4
// #pragma GCC ivdep
//         for (uint k = 1; k <= k_max; ++k) {
//           P_acc[k] *= P_array[alt_S][k - 1];
//         }
// #pragma GCC unroll 4
// #pragma GCC ivdep
//         for (uint k = 1; k <= k_max; ++k) {
//           P_array[S][k] += P_acc[k];
//         }
//       }
//     }
//   }
