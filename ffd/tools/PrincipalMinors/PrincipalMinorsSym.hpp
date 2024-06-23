namespace ffd::principal_minors {

template <uint order = 0, typename matrix_t>
auto Determinant_2_pow_n(matrix_t&& matrix) {
  using element_t =
      typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;

  const auto two_order = 2 * order;
  assert(two_order * two_order == (int)size(matrix));

  // DESCRIPTION
  //
  // Matrix:
  // |g11 g12 g13 g14 g15|
  // |g21 g22 g23 g24 g25|
  // |g31 g32 g33 g34 g35|
  // |g41 g42 g43 g44 g45|
  // |g51 g52 g53 g54 g55|
  //
  // In what follows we only write entries that are nonzero,
  // without specifying their coefficients.
  //
  // converting to binary:
  // (1,2,3,4) -> (1,2,4,8)
  // In the begining we have:
  // g11 = [a1,b1]                 equiv to (0,1)
  // g22 = [a2,b2]                 equiv to (0,2)
  // g33 = [a3,b3]                 equiv to (0,4)
  // g12 = [a1a2,b1a2,a1b2,b1b2]   equiv to (0,1,2,3)
  // g21 = [a2a1,b2a1,a2b1,b2b1]   equiv to (0,1,2,3) ->
  // g21 = [a1a2,b1a2,a1b2,b1b2]   equiv to (0,1,2,3)
  // g34 = [a3a4,b3a4,a3b4,b3b4]   equiv to (0,4,8,12)
  // g43 = [a3a4,b3a4,a3b4,b3b4]   equiv to (0,4,8,12)
  //
  // DIAGONAL:
  //
  // STEP 1:
  // g22 - g21 g11^-1 g12 =>
  // [a2,b2] (*)
  //([a1a2,b1a2,a1b2,b1b2] (*)
  // [a1,b1] (*)
  // [a1a2,b1a2,a1b2,b1b2]) =
  // [a1a2,b1a2,a1b2,b1b2]
  // computations
  // [0x0 ,1x1 ,0x2 ,1x3 ]
  // and
  // [0x0 ,0x1 ,1x2 ,1x3 ]
  //
  // STEP 2:
  // g33 - g32 g22^-1 g23 =>
  // [a1a3,b1a3,a1b3,b1b3] (*)
  // ([a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,b1b2a3,b1b2b3] (*)
  // [a1a2,b1a2,a1b2,b1b2] (*)
  // [a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,a1b2b3,b1b2b3]) =
  // [a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,a1b2b3,b1b2b3]
  // computations
  // [0x0   ,1x1   ,2x2   ,3x3   ,0x4   ,1x5   ,2x6   ,3x7   ]
  // and
  // [0x0   ,1x1   ,0x2   ,1x3   ,2x4   ,3x5   ,2x6   ,3x7   ]
  //
  // STEP 3:
  // g44 - g43 g33^-1 g34 =>
  // [a1a2a4,b1a2a4,a1b2a4,b1b2a4,a1a2b4,b1a2b4,a1b2b4,b1b2b4] (*)
  // [a1a2a3a4,b1a2a3a4,a1b2a3a4,b1b2a3a4,a1a2b3a4,b1a2b3a4,a1b2b3a4,b1b2b3a4,
  //  a1a2a3b4,b1a2a3b4,a1b2a3b4,b1b2a3b4,a1a2b3b4,b1a2b3b4,a1b2b3b4,b1b2b3b4]
  // (*)
  // [a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,a1b2b3,b1b2b3] (*)
  // [a1a2a3a4,b1a2a3a4,a1b2a3a4,b1b2a3a4,a1a2b3a4,b1a2b3a4,a1b2b3a4,b1b2b3a4,
  //  a1a2a3b4,b1a2a3b4,a1b2a3b4,b1b2a3b4,a1a2b3b4,b1a2b3b4,a1b2b3b4,b1b2b3b4]
  // computations
  // [0x0     ,1x1     ,2x2     ,3x3     ,4x4     ,5x5     ,6x6     ,7x7
  //  0x8     ,1x9    ,2x10     ,3x11    ,4x12    ,5x13    ,6x14    ,7x15 ]
  // and
  // [0x0     ,1x1     ,2x2     ,3x3     ,0x4     ,1x5     ,2x6     ,3x7
  //  4x8     ,5x9    ,6x10     ,7x11    ,4x12    ,5x13    ,6x14    ,7x15 ]
  //
  //  Rules:
  //  1. for off-diagonal g's increment by one up to n
  //  2. for diagonal inverted g's increment by one up n/2 twice
  //  3. for diagonal subtraction split into four groups, in each group
  //     increment n/2 steps twice
  //
  // OFF-DIAGONAL:
  //
  // STEP 1:
  // g34 - g31 g11^-1 g14 =>
  // [a3a4,b3a4,a3b4,b3b4] (*)
  //([a1a3,b1a3,a1b3,b1b3] (*)
  // [a1,b1] (*)
  // [a1a4,b1a4,a1b4,b1b4]) =
  // [a1a3a4,b1a3a4,a1b3a4,b1b3a4,a1a3b4,b1a3b4,a1b3b4,b1b3b4]
  // computations
  // [0x0   ,1x1   ,0x2   ,1x3]
  // and
  // [0x0   ,1x1   ,2x0   ,3x1   ,0x2   ,1x3   ,2x2   ,3x3   ]
  // and
  // [0x0   ,0x1   ,1x2   ,1x3   ,2x4   ,2x5   ,3x6   ,3x7   ]
  //
  // STEP 2:
  // g34 - g32 g22^-1 g24 =>
  // [a1a3a4,b1a3a4,a1b3a4,b1b3a4,a1a3b4,b1a3b4,a1b3b4,b1b3b4] (*)
  //([a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,a1b2b3,b1b2b3] (*)
  // [a1a2,b1a2,a1b2,b1b2] (*)
  // [a1a2a4,b1a2a4,a1b2a4,b1b2a4,a1a2b4,b1a2b4,a1b2b4,b1b2b4]) =
  // [a1a2a3a4,b1a2a3a4,a1b2a3a4,b1b2a3a4,a1a2b3a4,b1a2b3a4,a1b2b3a4,b1b2b3a4,
  //  a1a2a3b4,b1a2a3b4,a1b2a3b4,b1b2a3b4,a1a2b3b4,b1a2b3b4,a1b2b3b4,b1b2b3b4]
  // computations
  // [0x0     ,1x1     ,2x2     ,3x3     ,0x4     ,1x5     ,2x6     ,3x7 ]
  // and
  // [0x0     ,1x1     ,2x2     ,3x3     ,4x0     ,5x1     ,6x2     ,7x3
  //  0x4     ,1x5     ,2x6     ,3x7     ,4x4     ,5x5     ,6x6     ,7x7 ]
  // and
  // [0x0     ,1x1     ,0x2     ,1x3     ,2x4     ,3x5     ,2x6     ,3x7
  //  4x8     ,5x9     ,4x10    ,5x11    ,6x12    ,7x13    ,6x14    ,7x15 ]
  //
  // STEP 3:
  // g45 - g43 g33^-1 g35 =>
  // [a1a2a4a5,b1a2a4a5,a1b2a4a5,b1b2a4a5,a1a2b4a5,b1a2b4a5,a1b2b4a5,b1b2b4a5,
  //  a1a2a4b5,b1a2a4b5,a1b2a4b5,b1b2a4b5,a1a2b4b5,b1a2b4b5,a1b2b4b5,b1b2b4b5]
  // (*)
  //([a1a2a3a4,b1a2a3a4,a1b2a3a4,b1b2a3a4,a1a2b3a4,b1a2b3a4,a1b2b3a4,b1b2b3a4,
  //  a1a2a3b4,b1a2a3b4,a1b2a3b4,b1b2a3b4,a1a2b3b4,b1a2b3b4,a1b2b3b4,b1b2b3b4]
  // (*)
  // [a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,a1b2b3,b1b2b3]
  // (*)
  // [a1a2a3a5,b1a2a3a5,a1b2a3a5,b1b2a3a5,a1a2b3a5,b1a2b3a5,a1b2b3a5,b1b2b3a5,
  //  a1a2a3b5,b1a2a3b5,a1b2a3b5,b1b2a3b5,a1a2b3b5,b1a2b3b5,a1b2b3b5,b1b2b3b5])
  //  =
  // [a1a2a3a4a5,b1a2a3a4a5,a1b2a3a4a5,b1b2a3a4a5,a1a2b3a4a5,b1a2b3a4a5,a1b2b3a4a5,b1b2b3a4a5,
  // [a1a2a3b4a5,b1a2a3b4a5,a1b2a3b4a5,b1b2a3b4a5,a1a2b3b4a5,b1a2b3b4a5,a1b2b3b4a5,b1b2b3b4a5,
  //  a1a2a3a4b5,b1a2a3a4b5,a1b2a3a4b5,b1b2a3a4b5,a1a2b3a4b5,b1a2b3a4b5,a1b2b3a4b5,b1b2b3a4b5,
  //  a1a2a3b4b5,b1a2a3b4b5,a1b2a3b4b5,b1b2a3b4b5,a1a2b3b4b5,b1a2b3b4b5,a1b2b3b4b5,b1b2b3b4b5]
  // computations
  // [0x0,1x1,2x2,3x3,4x4,5x5,6x6,7x7,0x8,1x9,2x10,3x11,4x12,5x13,6x14,7x15]
  // and
  // [0x0,1x1,2x2,3x3,4x4,5x5,6x6,7x7,8x0,9x1,10x2,11x3,12x4,13x5,14x6,15x7,
  //  0x8,1x9,2x10,3x11,4x12,5x13,6x14,7x15,8x8,9x9,10x10,11x11,12x12,13x13,14x14,15x15]
  // and
  // [0x0,1x1,2x2,3x3,0x4,1x5,2x6,3x7,4x8,5x9,6x10,7x11,4x12,5x13,6x14,7x15,
  //  8x16,9x17,10x18,11x19,8x20,9x21,10x22,11x23,12x24,13x25,14x26,15x27,12x28,13x29,14x30,15x31]
  //
  // Rules:
  // 1. for inverted g's increment up to n/2 twice
  // 2. for left off-diagonal increment up to n/2 twice
  // 3. for right off-diagonal split into four groups, in each group
  //     increment n/2 steps twice
  // 4. diagonal g split into four groups, in each group increment n/4 steps
  //
  //
  //
  //
  //

  // correct minors
  // generate nxn matrix from 2nx2n matrix:
  // std::array<element_t, 1 << (order + 2)> m_c;
  // std::array<element_t, 1 << (order + 3)> m_up;
  // std::array<element_t, 1 << (order + 3)> m_dn;
  std::vector<element_t> m_c(1 << (order + 2));
  std::vector<element_t> m_up(1 << (order + 3));
  std::vector<element_t> m_dn(1 << (order + 3));

  for (uint j = 0; j < order; ++j) {
    const auto j_ = 2 * ((1 << j) - 1);
    m_c[j_] = matrix[2 * j * (two_order + 1)];
    m_c[j_ + 1] = matrix[(2 * j + 1) * (two_order + 1)];
  }
  for (uint j = 1; j < order; ++j) {
    const auto j2 = 2 * j;
    const auto j2_1 = j2 + 1;
    const auto j_ = (1 << j) - j - 1;
    for (uint k = 0; k < j; ++k) {
      const auto k2 = 2 * k;
      const auto k2_1 = k2 + 1;
      const auto k_ = (1 << k) - 1;
      const auto jk_ = 4 * (j_ + k_);
      m_up[jk_ + 0] = matrix[j2 + k2 * two_order];
      m_up[jk_ + 1] = matrix[j2 + k2_1 * two_order];
      m_up[jk_ + 2] = matrix[j2_1 + k2 * two_order];
      m_up[jk_ + 3] = matrix[j2_1 + k2_1 * two_order];
      m_dn[jk_ + 0] = matrix[k2 + j2 * two_order];
      m_dn[jk_ + 1] = matrix[k2_1 + j2 * two_order];
      m_dn[jk_ + 2] = matrix[k2 + j2_1 * two_order];
      m_dn[jk_ + 3] = matrix[k2_1 + j2_1 * two_order];
    }
  }

  // for storing minor results
  // std::array<element_t, (1 << (order + 1))> minors;
  std::vector<element_t> minors(1 << (order + 1));

  auto cnt = 2;
  for (uint S = 0; S < 2; ++S) {
    minors[cnt] = m_c[S];
    cnt++;
  }

  // vector sizes:
  // 2  4  4  4  4 ...
  // 4  4  8  8  8
  // 4  8  8 16 16
  // 4  8 16 16 32
  // 4  8 16 32 32
  // ...
  //
  // RULES:
  // diagonal (i>=0):
  // M_ii[S] = 2 * (2^i - 1) + S
  // 0 2 6 14 30 ...
  // off-diags (i>j>=1):
  // M_ij[S] = 4 * (2^i - i - 1 + 2^j - 1) + S
  // -  0  4 16 44 ...
  // -  -  8 20 48
  // -  -  - 28 56
  // -  -  -  - 88
  //
  // together:
  //  0  0  4 16 44
  //  0  2  8 20 48
  //  4  8  6 28 56
  // 16 20 28 14 88
  // 44 48 56 88 30
  //
  // Schur Complement loop
  for (auto l = 0; l < order - 1; ++l) {  // depth of tree
    const auto len_1 = (1 << l);
    const auto len_2 = 2 * len_1;
    const auto len_3 = 3 * len_1;
    const auto len_4 = 4 * len_1;
    const auto len_8 = 8 * len_1;
    const auto l_ = (1 << l) - 1;
    const auto ll = 2 * ((1 << l) - 1);

    // compute inverse
    for (auto S = 0; S < len_2; ++S) {
      m_c[ll + S] = 1. / m_c[ll + S];
    }

    // compute inverse * row
    for (auto j = l + 1; j < order; ++j) {
      const auto j_ = (1 << j) - j - 1;
      const auto jl = 4 * (j_ + l_);
      const auto p2 = jl + len_2;
      for (auto S = 0; S < len_2; ++S) {
        m_up[jl + S] *= m_c[ll + S];
        m_up[p2 + S] *= m_c[ll + S];
      }
    }

    // compute diagonal elements
    // [0x0,0x1,1x2,1x3]
    // [0x0,1x1,0x2,1x3,2x4,3x5,2x6,3x7]
    // [0x0,1x1,2x2,3x3,0x4,1x5,2x6,3x7,4x8,5x9,6x10,7x11,4x12,5x13,6x14,7x15]
    for (auto j = l + 1; j < order; ++j) {
      const auto j_ = (1 << j) - j - 1;
      const auto jl = 4 * (j_ + l_);
      const auto jj = 2 * ((1 << j) - 1);
      const auto p1 = jj - len_1;
      const auto p2 = jj - len_2;
#pragma GCC unroll 4
#pragma GCC ivdep
      for (auto S = len_3; S < len_4; ++S) {
        m_c[jj + S] = m_c[p2 + S] - m_up[jl + S] * m_dn[jl + S];
      }
#pragma GCC unroll 4
#pragma GCC ivdep
      for (auto S = len_2; S < len_3; ++S) {
        m_c[jj + S] = m_c[p1 + S] - m_up[jl + S] * m_dn[jl + S];
      }
#pragma GCC unroll 4
#pragma GCC ivdep
      for (auto S = len_1; S < len_2; ++S) {
        m_c[jj + S] = m_c[p1 + S] - m_up[jl + S] * m_dn[jl + S];
      }
#pragma GCC unroll 4
#pragma GCC ivdep
      for (auto S = 0; S < len_1; ++S) {
        m_c[jj + S] -= m_up[jl + S] * m_dn[jl + S];
      }
    }

    // compute off-diagonal elements (BOTTLENECK)
    // [0x0,1x1,2x0,3x1,0x2,1x3,2x2,3x3]
    // [0x0,1x1,2x2,3x3,4x0,5x1,6x2,7x3,0x4,1x5,2x6,3x7,4x4,5x5,6x6,7x7]
    // [0x0,1x1,2x2,3x3,4x4,5x5,6x6,7x7,8x0,9x1,10x2,11x3,12x4,13x5,14x6,15x7,
    //  0x8,1x9,2x10,3x11,4x12,5x13,6x14,7x15,8x8,9x9,10x10,11x11,12x12,13x13,14x14,15x15]
    //
    // [0x0,0x1,1x2,1x3,2x4,2x5,3x6,3x7]
    // 0 1 1 2 2 3 3 4
    // [0x0,1x1,0x2,1x3,2x4,3x5,2x6,3x7,4x8,5x9,4x10,5x11,6x12,7x13,6x14,7x15]
    // [0x0,1x1,2x2,3x3,0x4,1x5,2x6,3x7,4x8,5x9,6x10,7x11,4x12,5x13,6x14,7x15,
    //  8x16,9x17,10x18,11x19,8x20,9x21,10x22,11x23,12x24,13x25,14x26,15x27,12x28,13x29,14x30,15x31]
    for (auto j = l + 1; j < order; ++j) {
      const auto j_ = (1 << j) - j - 1;
      const auto lj = 4 * (j_ + l_);
      for (auto k = l + 1; k < j; ++k) {
        const auto k_ = (1 << k) - 1;
        const auto k__ = (1 << k) - k - 1;
        const auto lk = 4 * (l_ + k__);
        const auto jk = 4 * (j_ + k_);
        auto p2 = jk - len_4;
        auto p3 = lj - len_4;
        auto p4 = lk - len_4;
        auto S_max = 8 * len_1;
        auto S_min = S_max - len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
        for (auto S = S_min; S < S_max; ++S) {
          m_up[jk + S] = m_up[p2 + S] - m_up[p3 + S] * m_dn[p4 + S];
          m_dn[jk + S] = m_dn[p2 + S] - m_dn[p3 + S] * m_up[p4 + S];
        }
        S_min -= len_1;
        S_max -= len_1;
        p2 += len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
        for (auto S = S_min; S < S_max; ++S) {
          m_up[jk + S] = m_up[p2 + S] - m_up[p3 + S] * m_dn[p4 + S];
          m_dn[jk + S] = m_dn[p2 + S] - m_dn[p3 + S] * m_up[p4 + S];
        }
        S_min -= len_1;
        S_max -= len_1;
        p3 += len_2;
#pragma GCC unroll 4
#pragma GCC ivdep
        for (auto S = S_min; S < S_max; ++S) {
          m_up[jk + S] = m_up[p2 + S] - m_up[p3 + S] * m_dn[p4 + S];
          m_dn[jk + S] = m_dn[p2 + S] - m_dn[p3 + S] * m_up[p4 + S];
        }
        S_min -= len_1;
        S_max -= len_1;
        p2 += len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
        for (auto S = S_min; S < S_max; ++S) {
          m_up[jk + S] = m_up[p2 + S] - m_up[p3 + S] * m_dn[p4 + S];
          m_dn[jk + S] = m_dn[p2 + S] - m_dn[p3 + S] * m_up[p4 + S];
        }
        S_min -= len_1;
        S_max -= len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
        for (auto S = S_min; S < S_max; ++S) {
          m_up[jk + S] = m_up[p2 + S] - m_up[p3 + S] * m_dn[lk + S];
          m_dn[jk + S] = m_dn[p2 + S] - m_dn[p3 + S] * m_up[lk + S];
        }
        S_min -= len_1;
        S_max -= len_1;
        p2 += len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
        for (auto S = S_min; S < S_max; ++S) {
          m_up[jk + S] = m_up[p2 + S] - m_up[p3 + S] * m_dn[lk + S];
          m_dn[jk + S] = m_dn[p2 + S] - m_dn[p3 + S] * m_up[lk + S];
        }
        S_min -= len_1;
        S_max -= len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
        for (auto S = S_min; S < S_max; ++S) {
          m_up[jk + S] = m_up[p2 + S] - m_up[lj + S] * m_dn[lk + S];
          m_dn[jk + S] = m_dn[p2 + S] - m_dn[lj + S] * m_up[lk + S];
        }
        S_min -= len_1;
        S_max -= len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
        for (auto S = S_min; S < S_max; ++S) {
          m_up[jk + S] -= m_up[lj + S] * m_dn[lk + S];
          m_dn[jk + S] -= m_dn[lj + S] * m_up[lk + S];
        }
      }
    }

    const auto ll1 = 2 * ((1 << (l + 1)) - 1);
    for (auto S = 0; S < len_4; ++S) {
      minors[cnt] = m_c[ll1 + S];
      cnt++;
    }
  }
  // end main loop

  // [a1,b1] = [2,3]
  // [a1a2,b1a2,a1b2,b1b2] = [4,5,6,7]
  // [a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,a1b2b3,b1b2b3] =
  // [8,9,10,11,12,13,14,15]

  // correct minors
  for (auto l = 1; l < order; ++l) {
    for (auto S = (1 << l); S < (1 << (l + 1)); ++S) {
      minors[S + (1 << (l + 1))] *= minors[S];
      minors[S + (1 << l)] *= minors[S];
    }
  }

  return minors;
}

template <uint order = 0, typename matrix_t>
auto PrincipalMinors_2_pow_n(matrix_t&& matrix) {
  using element_t =
      typename std::decay<decltype(std::declval<matrix_t>()[0])>::type;
  const auto two_order = 2 * order;
  assert(two_order * two_order == (int)size(matrix));

  // all values for 3^n
  static constexpr std::array<int, 21> pow_3_n = {
      1,       3,        9,        27,        81,        243,       729,
      2187,    6561,     19683,    59049,     177147,    531441,    1594323,
      4782969, 14348907, 43046721, 129140163, 387420489, 1162261467};

  // compute a reasonable pivot coefficient
  element_t pivot_ = 0.0;
  for (auto j = 0; j < size(matrix); ++j) {
    pivot_ = std::max(abs(matrix[j]), pivot_);
  }
  // pivot_ /= (double)size(matrix);
  const element_t pivot = pivot_;
  // const element_t pivot = 0.;
  // std::cout << "pivot = " << pivot << std::endl;

  // DESCRIPTION
  //
  // Matrix:
  // |g11 g12 g13 g14 g15|
  // |g21 g22 g23 g24 g25|
  // |g31 g32 g33 g34 g35|
  // |g41 g42 g43 g44 g45|
  // |g51 g52 g53 g54 g55|
  //
  // In what follows we only write entries that are nonzero,
  // without specifying their coefficients.
  //
  // converting to binary:
  // (1,2,3,4) -> (1,2,4,8)
  // In the begining we have:
  // g11 = [a1,b1]                 equiv to (0,1)
  // g22 = [a2,b2]                 equiv to (0,2)
  // g33 = [a3,b3]                 equiv to (0,4)
  // g12 = [a1a2,b1a2,a1b2,b1b2]   equiv to (0,1,2,3)
  // g21 = [a2a1,b2a1,a2b1,b2b1]   equiv to (0,1,2,3) ->
  // g21 = [a1a2,b1a2,a1b2,b1b2]   equiv to (0,1,2,3)
  // g34 = [a3a4,b3a4,a3b4,b3b4]   equiv to (0,4,8,12)
  // g43 = [a3a4,b3a4,a3b4,b3b4]   equiv to (0,4,8,12)
  //
  // DIAGONAL:
  //
  // STEP 1:
  // g22 - g21 g11^-1 g12 =>
  // [a2,b2] (*)
  //([a1a2,b1a2,a1b2,b1b2] (*)
  // [a1,b1] (*)
  // [a1a2,b1a2,a1b2,b1b2]) =
  // [a1a2,b1a2,a1b2,b1b2]
  // computations
  // [0x0 ,1x1 ,0x2 ,1x3 ]
  // and
  // [0x0 ,0x1 ,1x2 ,1x3 ]
  //
  // STEP 2:
  // g33 - g32 g22^-1 g23 =>
  // [a1a3,b1a3,a1b3,b1b3] (*)
  // ([a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,b1b2a3,b1b2b3] (*)
  // [a1a2,b1a2,a1b2,b1b2] (*)
  // [a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,a1b2b3,b1b2b3]) =
  // [a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,a1b2b3,b1b2b3]
  // computations
  // [0x0   ,1x1   ,2x2   ,3x3   ,0x4   ,1x5   ,2x6   ,3x7   ]
  // and
  // [0x0   ,1x1   ,0x2   ,1x3   ,2x4   ,3x5   ,2x6   ,3x7   ]
  //
  // STEP 3:
  // g44 - g43 g33^-1 g34 =>
  // [a1a2a4,b1a2a4,a1b2a4,b1b2a4,a1a2b4,b1a2b4,a1b2b4,b1b2b4] (*)
  //
  // [a1a2a3a4,b1a2a3a4,a1b2a3a4,b1b2a3a4,a1a2b3a4,b1a2b3a4,a1b2b3a4,b1b2b3a4,
  // a1a2a3b4,b1a2a3b4,a1b2a3b4,b1b2a3b4,a1a2b3b4,b1a2b3b4,a1b2b3b4,b1b2b3b4]
  // (*)
  // [a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,a1b2b3,b1b2b3] (*)
  //
  // [a1a2a3a4,b1a2a3a4,a1b2a3a4,b1b2a3a4,a1a2b3a4,b1a2b3a4,a1b2b3a4,b1b2b3a4,
  // a1a2a3b4,b1a2a3b4,a1b2a3b4,b1b2a3b4,a1a2b3b4,b1a2b3b4,a1b2b3b4,b1b2b3b4]
  // computations
  // [0x0     ,1x1     ,2x2     ,3x3     ,4x4     ,5x5     ,6x6     ,7x7
  //  0x8     ,1x9    ,2x10     ,3x11    ,4x12    ,5x13    ,6x14    ,7x15 ]
  // and
  // [0x0     ,1x1     ,2x2     ,3x3     ,0x4     ,1x5     ,2x6     ,3x7
  //  4x8     ,5x9    ,6x10     ,7x11    ,4x12    ,5x13    ,6x14    ,7x15 ]
  //
  //  Rules:
  //  1. for off-diagonal g's increment by one up to n
  //  2. for diagonal inverted g's increment by one up n/2 twice
  //  3. for diagonal subtraction split into four groups, in each group
  //     increment n/2 steps twice
  //
  // OFF-DIAGONAL:
  //
  // STEP 1:
  // g34 - g31 g11^-1 g14 =>
  // [a3a4,b3a4,a3b4,b3b4] (*)
  //([a1a3,b1a3,a1b3,b1b3] (*)
  // [a1,b1] (*)
  // [a1a4,b1a4,a1b4,b1b4]) =
  // [a1a3a4,b1a3a4,a1b3a4,b1b3a4,a1a3b4,b1a3b4,a1b3b4,b1b3b4]
  // computations
  // [0x0   ,1x1   ,0x2   ,1x3]
  // and
  // [0x0   ,1x1   ,2x0   ,3x1   ,0x2   ,1x3   ,2x2   ,3x3   ]
  // and
  // [0x0   ,0x1   ,1x2   ,1x3   ,2x4   ,2x5   ,3x6   ,3x7   ]
  //
  // STEP 2:
  // g34 - g32 g22^-1 g24 =>
  // [a1a3a4,b1a3a4,a1b3a4,b1b3a4,a1a3b4,b1a3b4,a1b3b4,b1b3b4] (*)
  //([a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,a1b2b3,b1b2b3] (*)
  // [a1a2,b1a2,a1b2,b1b2] (*)
  // [a1a2a4,b1a2a4,a1b2a4,b1b2a4,a1a2b4,b1a2b4,a1b2b4,b1b2b4]) =
  //
  // [a1a2a3a4,b1a2a3a4,a1b2a3a4,b1b2a3a4,a1a2b3a4,b1a2b3a4,a1b2b3a4,b1b2b3a4,
  // a1a2a3b4,b1a2a3b4,a1b2a3b4,b1b2a3b4,a1a2b3b4,b1a2b3b4,a1b2b3b4,b1b2b3b4]
  // computations
  // [0x0     ,1x1     ,2x2     ,3x3     ,0x4     ,1x5     ,2x6     ,3x7 ]
  // and
  // [0x0     ,1x1     ,2x2     ,3x3     ,4x0     ,5x1     ,6x2     ,7x3
  //  0x4     ,1x5     ,2x6     ,3x7     ,4x4     ,5x5     ,6x6     ,7x7 ]
  // and
  // [0x0     ,1x1     ,0x2     ,1x3     ,2x4     ,3x5     ,2x6     ,3x7
  //  4x8     ,5x9     ,4x10    ,5x11    ,6x12    ,7x13    ,6x14    ,7x15 ]
  //
  // STEP 3:
  // g45 - g43 g33^-1 g35 =>
  //
  // [a1a2a4a5,b1a2a4a5,a1b2a4a5,b1b2a4a5,a1a2b4a5,b1a2b4a5,a1b2b4a5,b1b2b4a5,
  // a1a2a4b5,b1a2a4b5,a1b2a4b5,b1b2a4b5,a1a2b4b5,b1a2b4b5,a1b2b4b5,b1b2b4b5]
  // (*)
  //([a1a2a3a4,b1a2a3a4,a1b2a3a4,b1b2a3a4,a1a2b3a4,b1a2b3a4,a1b2b3a4,b1b2b3a4,
  // a1a2a3b4,b1a2a3b4,a1b2a3b4,b1b2a3b4,a1a2b3b4,b1a2b3b4,a1b2b3b4,b1b2b3b4]
  // (*)
  // [a1a2a3,b1a2a3,a1b2a3,b1b2a3,a1a2b3,b1a2b3,a1b2b3,b1b2b3]
  // (*)
  //
  // [a1a2a3a5,b1a2a3a5,a1b2a3a5,b1b2a3a5,a1a2b3a5,b1a2b3a5,a1b2b3a5,b1b2b3a5,
  //
  // a1a2a3b5,b1a2a3b5,a1b2a3b5,b1b2a3b5,a1a2b3b5,b1a2b3b5,a1b2b3b5,b1b2b3b5])
  //  =
  //
  // [a1a2a3a4a5,b1a2a3a4a5,a1b2a3a4a5,b1b2a3a4a5,a1a2b3a4a5,b1a2b3a4a5,a1b2b3a4a5,b1b2b3a4a5,
  //
  // [a1a2a3b4a5,b1a2a3b4a5,a1b2a3b4a5,b1b2a3b4a5,a1a2b3b4a5,b1a2b3b4a5,a1b2b3b4a5,b1b2b3b4a5,
  //
  // a1a2a3a4b5,b1a2a3a4b5,a1b2a3a4b5,b1b2a3a4b5,a1a2b3a4b5,b1a2b3a4b5,a1b2b3a4b5,b1b2b3a4b5,
  //
  // a1a2a3b4b5,b1a2a3b4b5,a1b2a3b4b5,b1b2a3b4b5,a1a2b3b4b5,b1a2b3b4b5,a1b2b3b4b5,b1b2b3b4b5]
  // computations
  // [0x0,1x1,2x2,3x3,4x4,5x5,6x6,7x7,0x8,1x9,2x10,3x11,4x12,5x13,6x14,7x15]
  // and
  // [0x0,1x1,2x2,3x3,4x4,5x5,6x6,7x7,8x0,9x1,10x2,11x3,12x4,13x5,14x6,15x7,
  //
  // 0x8,1x9,2x10,3x11,4x12,5x13,6x14,7x15,8x8,9x9,10x10,11x11,12x12,13x13,14x14,15x15]
  // and
  // [0x0,1x1,2x2,3x3,0x4,1x5,2x6,3x7,4x8,5x9,6x10,7x11,4x12,5x13,6x14,7x15,
  //
  // 8x16,9x17,10x18,11x19,8x20,9x21,10x22,11x23,12x24,13x25,14x26,15x27,12x28,13x29,14x30,15x31]
  //
  // Rules:
  // 1. for inverted g's increment up to n/2 twice
  // 2. for left off-diagonal increment up to n/2 twice
  // 3. for right off-diagonal split into four groups, in each group
  //     increment n/2 steps twice
  // 4. diagonal g split into four groups, in each group increment n/4 steps
  //
  //
  //
  //
  //

  // split matrix into diagonal and upper/lower triangles
  // one per level
  // std::array<std::array<element_t, 1 << (order + 2)>, order> m_c;
  // std::array<std::array<element_t, 1 << (order + 3)>, order> m_up;
  // std::array<std::array<element_t, 1 << (order + 3)>, order> m_dn;

  std::vector<std::vector<element_t>> m_c(order);
  std::vector<std::vector<element_t>> m_up(order);
  std::vector<std::vector<element_t>> m_dn(order);
  // std::array<std::vector<element_t>, order> m_c;
  // std::array<std::vector<element_t>, order> m_up;
  // std::array<std::vector<element_t>, order> m_dn;
  for (auto j = 0; j < order; ++j) {
    m_c[j].resize(1 << (order + 2));
    m_up[j].resize(1 << (order + 3));
    m_dn[j].resize(1 << (order + 3));
  }

  // generate nxn polynomial matrix from 2nx2n matrix:
  for (auto j = 0; j < order; ++j) {
    const auto j_ = 2 * ((1 << j) - 1);
    m_c[0][j_] = matrix[2 * j * (two_order + 1)];
    m_c[0][j_ + 1] = matrix[(2 * j + 1) * (two_order + 1)];
  }
  for (auto j = 1; j < order; ++j) {
    const auto j2 = 2 * j;
    const auto j2_1 = j2 + 1;
    const auto j_ = (1 << j) - j - 1;
    for (auto k = 0; k < j; ++k) {
      const auto k2 = 2 * k;
      const auto k2_1 = k2 + 1;
      const auto k_ = (1 << k) - 1;
      const auto jk_ = 4 * (j_ + k_);
      m_up[0][jk_ + 0] = matrix[j2 + k2 * two_order];
      m_up[0][jk_ + 1] = matrix[j2 + k2_1 * two_order];
      m_up[0][jk_ + 2] = matrix[j2_1 + k2 * two_order];
      m_up[0][jk_ + 3] = matrix[j2_1 + k2_1 * two_order];
      m_dn[0][jk_ + 0] = matrix[k2 + j2 * two_order];
      m_dn[0][jk_ + 1] = matrix[k2_1 + j2 * two_order];
      m_dn[0][jk_ + 2] = matrix[k2 + j2_1 * two_order];
      m_dn[0][jk_ + 3] = matrix[k2_1 + j2_1 * two_order];
    }
  }

  // array for storing minor results of size 3^n+1
  // std::array<element_t, pow_3_n[order] + 1> minors;
  std::vector<element_t> minors(pow_3_n[order] + 1);
  minors[0] = 1.;
  // index 2^n array positions and subset popcounts
  std::array<uint, (1 << order)> pos, pop, pop_1;
  pos[0] = 0;
  for (auto S = 0; S < (1 << order) - 1; ++S) {
    pop[S] = __builtin_popcount(S);
    pop_1[S] = (1 << pop[S]);
    pos[S + 1] = pos[S] + pop_1[S];
  }
  pop[(1 << order) - 1] = __builtin_popcount((1 << order) - 1);
  pop_1[(1 << order) - 1] = (1 << pop[(1 << order) - 1]);
  // store minors of subset 1
  for (auto S = 0; S < 2; ++S) {
    minors[pos[1] + S] = m_c[0][S];
  }

  // vector sizes:
  // 2  4  4  4  4 ...
  // 4  4  8  8  8
  // 4  8  8 16 16
  // 4  8 16 16 32
  // 4  8 16 32 32
  // ...
  //
  // RULES:
  // diagonal (i>=0):
  // M_ii[S] = 2 * (2^i - 1) + S
  // 0 2 6 14 30 ...
  // off-diags (i>j>=1):
  // M_ij[S] = 4 * (2^i - i - 1 + 2^j - 1) + S
  // -  0  4 16 44 ...
  // -  -  8 20 48
  // -  -  - 28 56
  // -  -  -  - 88
  //
  // together:
  //  0  0  4 16 44
  //  0  2  8 20 48
  //  4  8  6 28 56
  // 16 20 28 14 88
  // 44 48 56 88 30
  //
  // Schur Complement loop

  auto lev = 0;                         // tree level
  auto set = 1;                         // current subset
  std::array<long int, order> visited;  // last visites set at each level
  std::array<long int, order> tree_l;   // polynomial order at each level
  visited.fill(0);
  visited[0] = 1;
  tree_l.fill(0);

  // BINARY TREE
  while (visited[order - 1] != (1 << order) - 1) {
    // precompute some values
    const auto l = lev;
    const auto len_1 = (1 << tree_l[lev]);
    const auto len_2 = 2 * len_1;
    const auto len_3 = 3 * len_1;
    const auto len_4 = 4 * len_1;
    const auto len_8 = 8 * len_1;
    const auto l_ = (1 << l) - 1;
    const auto ll = 2 * ((1 << l) - 1);
    const auto ll1 = 2 * ((1 << (l + 1)) - 1);

    // move
    if (lev == order - 1)  // bottom
    {
      // move up in tree
      lev--;
      set = visited[lev];
    } else if (visited[lev + 1] < set + (1 << (lev + 1)) - (1 << lev))  // left
    {
      // copy sub-matrix
      {
        // diagonal
        for (auto j = l + 1; j < order; ++j) {
          const auto jj = 2 * ((1 << j) - 1);
          for (auto S = 0; S < len_2; ++S) {
            m_c[l + 1][jj + S] = m_c[l][jj + S];
          }
        }
        // off-diagonal
        for (auto j = l + 1; j < order; ++j) {
          const auto j_ = (1 << j) - j - 1;
          for (auto k = l + 1; k < j; ++k) {
            const auto k_ = (1 << k) - 1;
            const auto jk = 4 * (j_ + k_);
            for (auto S = 0; S < len_4; ++S) {
              m_up[l + 1][jk + S] = m_up[l][jk + S];
              m_dn[l + 1][jk + S] = m_dn[l][jk + S];
            }
          }
        }
      }

      // move left in tree
      set += (1 << (lev + 1)) - (1 << lev);
      lev++;
      visited[lev] = set;
      tree_l[lev] = tree_l[lev - 1];

      // copy minors
      for (auto S = 0; S < len_2; ++S) {
        minors[pos[set] + S] = m_c[l + 1][ll1 + S];
      }

    } else if (visited[lev + 1] < set + (1 << (lev + 1)))  // right
    {
      // do schur complement
      {
        if (__builtin_popcount(set) == 1) {
          for (auto S = 0; S < len_2; ++S) {
            m_c[l][ll + S] += pivot;
          }
        }
        // compute inverse
        for (auto S = 0; S < len_2; ++S) {
          m_c[l + 1][ll + S] = 1. / m_c[l][ll + S];
          // if (__builtin_popcount(set) == 1)
          //   std::cout << S << " " << m_c[l + 1][ll + S] << std::endl;
        }

        // compute inverse * row
        for (auto j = l + 1; j < order; ++j) {
          const auto j_ = (1 << j) - j - 1;
          const auto jl = 4 * (j_ + l_);
          const auto p2 = jl + len_2;
          for (auto S = 0; S < len_2; ++S) {
            m_up[l + 1][jl + S] = m_up[l][jl + S] * m_c[l + 1][ll + S];
            m_up[l + 1][p2 + S] = m_up[l][p2 + S] * m_c[l + 1][ll + S];
          }
        }

        // compute diagonal elements
        // [0x0,0x1,1x2,1x3]
        // [0x0,1x1,0x2,1x3,2x4,3x5,2x6,3x7]
        //
        // [0x0,1x1,2x2,3x3,0x4,1x5,2x6,3x7,4x8,5x9,6x10,7x11,4x12,5x13,6x14,7x15]
        for (auto j = l + 1; j < order; ++j) {
          const auto j_ = (1 << j) - j - 1;
          const auto jl = 4 * (j_ + l_);
          const auto jj = 2 * ((1 << j) - 1);
          const auto p1 = jj - len_1;
          const auto p2 = jj - len_2;
#pragma GCC unroll 4
#pragma GCC ivdep
          for (auto S = len_3; S < len_4; ++S) {
            m_c[l + 1][jj + S] =
                m_c[l][p2 + S] - m_up[l + 1][jl + S] * m_dn[l][jl + S];
          }
#pragma GCC unroll 4
#pragma GCC ivdep
          for (auto S = len_2; S < len_3; ++S) {
            m_c[l + 1][jj + S] =
                m_c[l][p1 + S] - m_up[l + 1][jl + S] * m_dn[l][jl + S];
          }
#pragma GCC unroll 4
#pragma GCC ivdep
          for (auto S = len_1; S < len_2; ++S) {
            m_c[l + 1][jj + S] =
                m_c[l][p1 + S] - m_up[l + 1][jl + S] * m_dn[l][jl + S];
          }
#pragma GCC unroll 4
#pragma GCC ivdep
          for (auto S = 0; S < len_1; ++S) {
            m_c[l + 1][jj + S] =
                m_c[l][jj + S] - m_up[l + 1][jl + S] * m_dn[l][jl + S];
          }
        }

        // compute off-diagonal elements (BOTTLENECK)
        // [0x0,1x1,2x0,3x1,0x2,1x3,2x2,3x3]
        // [0x0,1x1,2x2,3x3,4x0,5x1,6x2,7x3,0x4,1x5,2x6,3x7,4x4,5x5,6x6,7x7]
        //
        // [0x0,1x1,2x2,3x3,4x4,5x5,6x6,7x7,8x0,9x1,10x2,11x3,12x4,13x5,14x6,15x7,
        //
        // 0x8,1x9,2x10,3x11,4x12,5x13,6x14,7x15,8x8,9x9,10x10,11x11,12x12,13x13,14x14,15x15]
        //
        // [0x0,0x1,1x2,1x3,2x4,2x5,3x6,3x7]
        // 0 1 1 2 2 3 3 4
        //
        // [0x0,1x1,0x2,1x3,2x4,3x5,2x6,3x7,4x8,5x9,4x10,5x11,6x12,7x13,6x14,7x15]
        //
        // [0x0,1x1,2x2,3x3,0x4,1x5,2x6,3x7,4x8,5x9,6x10,7x11,4x12,5x13,6x14,7x15,
        //
        // 8x16,9x17,10x18,11x19,8x20,9x21,10x22,11x23,12x24,13x25,14x26,15x27,12x28,13x29,14x30,15x31]
        for (auto j = l + 1; j < order; ++j) {
          const auto j_ = (1 << j) - j - 1;
          const auto lj = 4 * (j_ + l_);
          for (auto k = l + 1; k < j; ++k) {
            const auto k_ = (1 << k) - 1;
            const auto k__ = (1 << k) - k - 1;
            const auto lk = 4 * (l_ + k__);
            const auto jk = 4 * (j_ + k_);
            auto p2 = jk - len_4;
            auto p3 = lj - len_4;
            auto p4 = lk - len_4;
            auto S_max = 8 * len_1;
            auto S_min = S_max - len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
            for (auto S = S_min; S < S_max; ++S) {
              m_up[l + 1][jk + S] =
                  m_up[l][p2 + S] - m_up[l + 1][p3 + S] * m_dn[l][p4 + S];
              m_dn[l + 1][jk + S] =
                  m_dn[l][p2 + S] - m_dn[l][p3 + S] * m_up[l + 1][p4 + S];
            }
            S_min -= len_1;
            S_max -= len_1;
            p2 += len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
            for (auto S = S_min; S < S_max; ++S) {
              m_up[l + 1][jk + S] =
                  m_up[l][p2 + S] - m_up[l + 1][p3 + S] * m_dn[l][p4 + S];
              m_dn[l + 1][jk + S] =
                  m_dn[l][p2 + S] - m_dn[l][p3 + S] * m_up[l + 1][p4 + S];
            }
            S_min -= len_1;
            S_max -= len_1;
            p3 += len_2;
#pragma GCC unroll 4
#pragma GCC ivdep
            for (auto S = S_min; S < S_max; ++S) {
              m_up[l + 1][jk + S] =
                  m_up[l][p2 + S] - m_up[l + 1][p3 + S] * m_dn[l][p4 + S];
              m_dn[l + 1][jk + S] =
                  m_dn[l][p2 + S] - m_dn[l][p3 + S] * m_up[l + 1][p4 + S];
            }
            S_min -= len_1;
            S_max -= len_1;
            p2 += len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
            for (auto S = S_min; S < S_max; ++S) {
              m_up[l + 1][jk + S] =
                  m_up[l][p2 + S] - m_up[l + 1][p3 + S] * m_dn[l][p4 + S];
              m_dn[l + 1][jk + S] =
                  m_dn[l][p2 + S] - m_dn[l][p3 + S] * m_up[l + 1][p4 + S];
            }
            S_min -= len_1;
            S_max -= len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
            for (auto S = S_min; S < S_max; ++S) {
              m_up[l + 1][jk + S] =
                  m_up[l][p2 + S] - m_up[l + 1][p3 + S] * m_dn[l][lk + S];
              m_dn[l + 1][jk + S] =
                  m_dn[l][p2 + S] - m_dn[l][p3 + S] * m_up[l + 1][lk + S];
            }
            S_min -= len_1;
            S_max -= len_1;
            p2 += len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
            for (auto S = S_min; S < S_max; ++S) {
              m_up[l + 1][jk + S] =
                  m_up[l][p2 + S] - m_up[l + 1][p3 + S] * m_dn[l][lk + S];
              m_dn[l + 1][jk + S] =
                  m_dn[l][p2 + S] - m_dn[l][p3 + S] * m_up[l + 1][lk + S];
            }
            S_min -= len_1;
            S_max -= len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
            for (auto S = S_min; S < S_max; ++S) {
              m_up[l + 1][jk + S] =
                  m_up[l][p2 + S] - m_up[l + 1][lj + S] * m_dn[l][lk + S];
              m_dn[l + 1][jk + S] =
                  m_dn[l][p2 + S] - m_dn[l][lj + S] * m_up[l + 1][lk + S];
            }
            S_min -= len_1;
            S_max -= len_1;
#pragma GCC unroll 4
#pragma GCC ivdep
            for (auto S = S_min; S < S_max; ++S) {
              m_up[l + 1][jk + S] =
                  m_up[l][jk + S] - m_up[l + 1][lj + S] * m_dn[l][lk + S];
              m_dn[l + 1][jk + S] =
                  m_dn[l][jk + S] - m_dn[l][lj + S] * m_up[l + 1][lk + S];
            }
          }
        }
      }

      // move right in tree
      set += (1 << (lev + 1));
      lev++;
      visited[lev] = set;
      tree_l[lev] = tree_l[lev - 1] + 1;

      // copy minors
      for (auto S = 0; S < len_4; ++S) {
        minors[pos[set] + S] = m_c[l + 1][ll1 + S];
      }

    } else  // up
    {
      // move up in tree
      visited[lev + 1] = 0;
      lev--;
      set = visited[lev];
    }
  }

  // for (auto S = 0; S < (1 << order); ++S) {
  //   for (auto j = 0; j < pop[S]; ++j) {
  //     std::cout << S << " " << j << " " << minors[pos[S] + j] << std::endl;
  //   }
  // }

  // add pivots to minors with popcount == 1
  for (auto set = 1; set < (1 << order); set *= 2) {
    const auto pos_1 = pos[set];
    for (auto S = 0; S < 2; ++S) {
      minors[pos_1 + S] += pivot;
    }
  }

  // multiply by elements from previous levels
  auto set_K = 1;
  for (auto l = 0; l < order; ++l) {
    const auto set_2K = 2 * set_K;

    for (auto set = set_K + 1; set < set_2K; ++set) {
      const auto pos_1 = pos[set];
      const auto pop = pop_1[set - set_K];
      const auto pos_1_ = pos_1 + pop;
      const auto pos_2 = pos[set - set_K];
      for (auto S = 0; S < pop; ++S) {
        minors[pos_1 + S] *= minors[pos_2 + S];
        minors[pos_1_ + S] *= minors[pos_2 + S];
      }
    }
    set_K = set_2K;
  }

  // pivot corrections
  for (auto set = (1 << (order - 1)); set >= 1; set /= 2) {
    for (auto set_ = set; set_ < (1 << order); set_ += 2 * set) {
      const auto pos_1 = pos[set_];
      const auto pos_2 = pos[set_ - set];
      const auto pop = pop_1[set_ - set];
      const auto pos_1_ = pos_1 + pop;
      for (auto S = 0; S < pop; ++S) {
        minors[pos_1 + 2 * S] -= pivot * minors[pos_2 + S];
        minors[pos_1 + 2 * S + 1] -= pivot * minors[pos_2 + S];
      }
    }
  }

  // Group minors
  // std::array<element_t, (1 << order)> res_minors;
  // res_minors.fill(0.);
  // for (auto S = 0; S < (1 << order); ++S) {
  //   for (auto j = 0; j < pop_1[S]; ++j) {
  //     res_minors[S] += minors[pos[S] + j];
  //   }
  //   res_minors[S] /= (double)pop_1[S];
  // }

  // // pivot corrections full?
  // for (auto set = (1 << order) - 1; set > 0; --set) {
  //   if (pivots[set]) {
  //     const auto set_K = (1 << (31 - __builtin_clz(set)));
  //     const auto set_K_2 = 2 * set_K;
  //     for (uint set_ = set; set_ < (1 << order); set_ += set_K_2) {
  //       const auto pos_1 = pos[set_];
  //       const auto pop = pop_1[set_ - set_K];
  //       const auto pos_1_ = pos_1 + pop;
  //       const auto pos_2 = pos[set_ - set_K];
  //       for (auto S = 0; S < pop; ++S) {
  //         minors[set_ + S] -= pivot * minors[set_ - set_K + S];
  //       }
  //     }
  //   }
  // }

  // return res_minors;
  return minors;
}

}  // namespace ffd::principal_minors