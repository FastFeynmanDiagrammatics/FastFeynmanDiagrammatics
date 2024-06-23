namespace ffd::truncated_polynomial{

  template<uint n, bool not_full=false, class field_d = Real>
  class P{
  public:
    std::array<field_d, n+1> coef;
    uint order = 0;

    P(){coef.fill(0.); order = 0;}

    P(int k, field_d x){
      coef.fill(0.);
      order=k;
      for(uint j=0;j<=order;++j){
        coef[j]=x;
      }
    }

    P(field_d x){coef.fill(0); coef[0] = x; order = 0;}

    field_d& operator[](std::size_t j){return coef[j];}

    field_d operator[](std::size_t j) const{return coef[j];}

    field_d operator()(field_d x) const{
	    field_d ret = 0.;
	     field_d pow_x = 1.;
	    for(uint j=0; j<=order; ++j){
		    ret += coef[j]*pow_x;
		   pow_x*=x;
	    }
	    return ret;
    }

  };

  template<uint n, bool not_full, class field_d>
  auto operator+(P<n, not_full, field_d> const& P1, P<n, not_full, field_d> const& P2){
    auto const order_ret = std::max(P1.order, P2.order);
    auto ret = P1;
    ret.order = order_ret;
    for(uint j=0; j<=P2.order; ++j){
      ret[j] += P2[j];
    }
    return ret;
  }

  template<uint n, bool not_full, class field_d>
  auto& operator+=(P<n, not_full, field_d>& P1, P<n, not_full, field_d> const& P2){
    P1 = P1 + P2;
    return P1;
  }

  template<uint n, bool not_full, class field_d>
  auto operator-(P<n, not_full, field_d> const& P1, P<n, not_full, field_d> const& P2){
    auto const order_ret = std::max(P1.order, P2.order);
    auto ret = P1;
    ret.order = order_ret;
    for(uint j=0; j<=P2.order; ++j){
      ret[j] -= P2[j];
    }
    return ret;
  }

  template<uint n, bool not_full, class field_d>
  auto& operator-=(P<n, not_full, field_d>& P1, P<n, not_full, field_d> const& P2){
    P1 = P1 - P2;
    return P1;
  }

  template<uint n, bool not_full, class field_d>
  auto operator*(P<n, not_full, field_d> const& P1, P<n, not_full, field_d> const& P2){
    // if constexpr (not_full){
    //   auto const order_ret = std::min(P1.order+P2.order, n);
    //   auto ret = P<n, not_full, field_d>(order_ret, 0.);
    //   for(uint j1=0; j1<=P1.order; ++j1){
    //     for(uint j2=0; j2<=std::min(order_ret-j1, P2.order); ++j2){
    //       ret[j1+j2] += P1[j1]*P2[j2];
    //     }
    //   }
    //   return ret;
    // }
    // else{
      auto ret = P1;
      ret.order = std::min(P1.order+P2.order, n);

      if constexpr (n>=0){
        ret[0] = P1[0]*P2[0];
      }
      if constexpr (n>=1){
        ret[1] = P1[1]*P2[0]+P1[0]*P2[1];
      }
      if constexpr (n>=2){
        ret[2] = P1[2]*P2[0]+P1[1]*P2[1]+P1[0]*P2[2];
      }
      if constexpr (n>=3){
        ret[3] = P1[3]*P2[0]+P1[2]*P2[1]+P1[1]*P2[2]+P1[0]*P2[3];
      }
      if constexpr (n>=4){
        ret[4] = P1[4]*P2[0]+P1[3]*P2[1]+P1[2]*P2[2]+P1[1]*P2[3]
                +P1[0]*P2[4];
      }
      if constexpr (n>=5){
        ret[5] = P1[5]*P2[0]+P1[4]*P2[1]+P1[3]*P2[2]+P1[2]*P2[3]
                +P1[1]*P2[4]+P1[0]*P2[5];
      }
      if constexpr (n>=6){
        ret[6] = P1[6]*P2[0]+P1[5]*P2[1]+P1[4]*P2[2]+P1[3]*P2[3]
                +P1[2]*P2[4]+P1[1]*P2[5]+P1[0]*P2[6];
      }
      if constexpr (n>=7){
        ret[7] = P1[7]*P2[0]+P1[6]*P2[1]+P1[5]*P2[2]+P1[4]*P2[3]
                +P1[3]*P2[4]+P1[2]*P2[5]+P1[1]*P2[6]+P1[0]*P2[7];
      }
      if constexpr (n>=8){
        ret[8] = P1[8]*P2[0]+P1[7]*P2[1]+P1[6]*P2[2]+P1[5]*P2[3]
                +P1[4]*P2[4]+P1[3]*P2[5]+P1[2]*P2[6]+P1[1]*P2[7]
                +P1[0]*P2[8];
      }
      if constexpr (n>=9){
        ret[9] = P1[9]*P2[0]+P1[8]*P2[1]+P1[7]*P2[2]+P1[6]*P2[3]
                +P1[5]*P2[4]+P1[4]*P2[5]+P1[3]*P2[6]+P1[2]*P2[7]
                +P1[1]*P2[8]+P1[0]*P2[9];
      }
      if constexpr (n>=10){
        ret[10]= P1[10]*P2[0]+P1[9]*P2[1]+P1[8]*P2[2]+P1[7]*P2[3]
                +P1[6]*P2[4]+P1[5]*P2[5]+P1[4]*P2[6]+P1[3]*P2[7]
                +P1[2]*P2[8]+P1[1]*P2[9]+P1[0]*P2[10];
      }

      return ret;
    // }
  }

  template<uint n, bool not_full, class field_d>
  auto& operator*=(P<n, not_full, field_d>& P1, P<n, not_full, field_d> const& P2){
    P1 = P1*P2;
    return P1;
  }

  template<uint n, bool not_full, class field_d>
  auto operator/(P<n, not_full, field_d> const& P1, field_d x){
    auto ret = P1;
    for(uint j=0; j<=P1.order; ++j){
      ret[j] /= x;
    }
    return ret;
  }

  template<uint n, bool not_full, class field_d>
  auto& operator/=(P<n, not_full, field_d> & P1, field_d x){
    P1 = P1/x;
    return P1;
  }

  template<uint n, bool not_full, class field_d>
  auto operator/(field_d x, P<n, not_full, field_d> const& P1){
    P<n, not_full, field_d> P2 = x;
    auto ret = P2/P1;
    return ret;
  }

  template<uint n, bool not_full, class field_d>
  auto operator/(P<n, not_full, field_d> const& P1, P<n, not_full, field_d> const& P3){
    auto ret = P1/P3[0];
    auto const P2 = P3/P3[0];
    ret.order = n;

    if constexpr (n>=1){
      ret[1] -= P2[1]*ret[0];
    }
    if constexpr (n>=2){
      ret[2] -= P2[2]*ret[0] + P2[1]*ret[1];
    }
    if constexpr (n>=3){
      ret[3] -= P2[3]*ret[0] + P2[2]*ret[1] + P2[1]*ret[2];
    }
    if constexpr (n>=4){
      ret[4] -= P2[4]*ret[0] + P2[3]*ret[1] + P2[2]*ret[2] + P2[1]*ret[3];
    }
    if constexpr (n>=5){
      ret[5] -= P2[5]*ret[0] + P2[4]*ret[1] + P2[3]*ret[2] + P2[2]*ret[3] + P2[1]*ret[4];
    }
    if constexpr (n>=6){
      ret[6] -= P2[6]*ret[0] + P2[5]*ret[1] + P2[4]*ret[2] + P2[3]*ret[3] + P2[2]*ret[4] + P2[1]*ret[5];
    }
    if constexpr (n>=7){
      ret[7] -= P2[7]*ret[0] + P2[6]*ret[1] + P2[5]*ret[2] + P2[4]*ret[3] + P2[3]*ret[4] + P2[2]*ret[5] + P2[1]*ret[6];
    }
    if constexpr (n>=8){
      ret[8] -= P2[8]*ret[0] + P2[7]*ret[1] + P2[6]*ret[2] + P2[5]*ret[3] + P2[4]*ret[4] + P2[3]*ret[5] + P2[2]*ret[6] + P2[1]*ret[7];
    }
    if constexpr (n>=9){
      ret[9] -= P2[9]*ret[0] + P2[8]*ret[1] + P2[7]*ret[2] + P2[6]*ret[3] + P2[5]*ret[4] + P2[4]*ret[5] + P2[3]*ret[6] + P2[2]*ret[7] + P2[1]*ret[8];
    }
    if constexpr (n>=10){
      ret[10] -= P2[10]*ret[0] + P2[9]*ret[1] + P2[8]*ret[2] + P2[7]*ret[3] + P2[6]*ret[4] + P2[5]*ret[5] + P2[4]*ret[6] + P2[3]*ret[7] + P2[2]*ret[8] + P2[1]*ret[9];
    }

    return ret;
  }

  template<uint n, bool not_full, class field_d>
  auto& operator/=(P<n, not_full, field_d>& P1, P<n, not_full, field_d> const& P2){
    P1 = P1/P2;
    return P1;
  }

  template<uint n, bool not_full, class field_d>
  auto operator*(P<n, not_full, field_d> const& P1, field_d x){
    auto ret = P1;
    for(uint j=0; j<=P1.order; ++j){
      ret[j] *= x;
    }
    return ret;
  }

  template<uint n, bool not_full, class field_d>
  auto operator*(field_d x,P<n, not_full, field_d> const& P1){
    return P1*x;
  }

  template<uint n, bool not_full, class field_d>
  auto operator*=(P<n, not_full, field_d> const& P1, field_d x){
    auto ret = P1*x;
    return ret;
  }

  template<uint n, bool not_full, class field_d>
  auto operator+(P<n, not_full, field_d> const& P1, field_d x){
    auto ret = P1;
    ret[0] += x;
    return ret;
  }

  template<uint n, bool not_full, class field_d>
  auto operator+=(P<n, not_full, field_d> const& P1, field_d x){
    auto ret = P1+x;
    return ret;
  }

  template<uint n, bool not_full, class field_d>
  auto operator-(P<n, not_full, field_d> const& P1, field_d x){
    auto ret = P1;
    ret[0] -= x;
    return ret;
  }

  template<uint n, bool not_full, class field_d>
  auto operator-=(P<n, not_full, field_d> const& P1, field_d x){
    auto ret = P1-x;
    return ret;
  }

  template<uint n, bool not_full, class field_d>
  auto operator-(P<n, not_full, field_d> const& P1){
    return P1 * -1.;
  }

  template<uint n, bool not_full, class field_d>
  auto operator+(P<n, not_full, field_d> const& P1){
    return P1;
  }

  template<uint n, bool not_full, class field_d>
  auto operator<(P<n, not_full, field_d> const& P1, P<n, not_full, field_d> const& P2){
    return (P1[0]<P2[0]);
  }

  template<uint n, bool not_full, class field_d>
  auto operator>(P<n, not_full, field_d> const& P1, P<n, not_full, field_d> const& P2){
    return (P1[0]>P2[0]);
  }

  template<uint n, bool not_full, class field_d>
  auto operator>(P<n, not_full, field_d> const& P1, field_d x){
    return (P1[0]>x);
  }

  template<uint n, bool not_full, class field_d>
  auto operator<(P<n, not_full, field_d> const& P1, field_d x){
    return (P1[0]<x);
  }

  template<uint n, bool not_full, class field_d>
  auto abs(P<n, not_full, field_d> const& P1){
    auto ret = P1;
    for (uint j=0; j<=P1.order; ++j){
	    ret[j] = std::abs(ret[j]);
    }
    return ret;
  }

}//namespace
