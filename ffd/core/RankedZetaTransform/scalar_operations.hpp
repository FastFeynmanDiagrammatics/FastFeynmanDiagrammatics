namespace ffd::ranked_zeta{

	// template<typename Field, SetT SetType>
    // void
    // RankedZeta<Field, SetType>::operator=(const Field number_) noexcept{
	//   cardinality = 1;
	//   rank_size=1;
    //   (*this)[0] = number_;
    //   return *this;
    // }

	template<typename Field, SetT SetType>
    RankedZeta<Field, SetType>&
    RankedZeta<Field, SetType>::operator+=(const Field number_) noexcept{
      (*this)[0] += number_;
      return *this;
    }


    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType>& RankedZeta<Field, SetType>::operator-=(const Field number_) noexcept{
      return *this += -number_;
    }


    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType> RankedZeta<Field, SetType>::operator+(const Field number_) const noexcept{
      RankedZeta<Field, SetType> ret = *this;
      return ret += number_;
    }


    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType> RankedZeta<Field, SetType>::operator-(const Field number_) const noexcept{
      RankedZeta<Field, SetType> ret = *this;
      return ret -= number_;
    }


    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType>& RankedZeta<Field, SetType>::operator*=(const Field number_) noexcept{
      if(this->NotZero){
        for(BinaryInt S=0; S < this->size(); ++S){
  	this->coef[S] *= number_;
        }
      }
      return *this;
    }


    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType> RankedZeta<Field, SetType>::operator*(const Field number_) const noexcept{
      RankedZeta<Field, SetType> ret = *this;
      return ret *= number_;
    }


    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType>& RankedZeta<Field, SetType>::operator/=(const Field number_) noexcept{
      *this *= 1./number_;
      return *this;
    }


    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType> RankedZeta<Field, SetType>::operator/(const Field number_) const noexcept{
      RankedZeta<Field, SetType> ret = *this;
      return ret /= number_;
    }


    //with this we can write -P[[\xi]]
    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType>
    RankedZeta<Field, SetType>::operator-() const noexcept{
      auto ret = *this;
      if( ret.NotZero ){
        for( std::size_t j = 0ul; j < ret.size(); ++j){
  	ret.coef[j] = -ret.coef[j];
        }
      }
      return ret;
    }


    //with this we can write +P[[\xi]]
    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType>
    RankedZeta<Field, SetType>::operator+() const noexcept{
      return *this;
    }


    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType> operator+(Field number_, RankedZeta<Field, SetType> const& zeta_) noexcept{
      return zeta_+number_;
    }


    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType> operator-(Field number_, RankedZeta<Field, SetType> const& zeta_) noexcept{
      return -zeta_ + number_;
    }


    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType> operator*(Field number_, RankedZeta<Field, SetType> const& zeta_) noexcept{
      return zeta_ * number_;
    }

    template<typename Field, SetT SetType>
    RankedZeta<Field, SetType> operator/(Field number_, RankedZeta<Field, SetType> const& zeta_) noexcept{
      RankedZeta<Field, SetType> ret(__builtin_popcount(zeta_.cardinality-1)); //this is necessary because we have no clue
	  ret.rank_size = zeta_.rank_size;
  	  ret.cardinality = zeta_.cardinality;

	  if (SetType==Even){
	    ret.resize(zeta_.cardinality*(zeta_.rank_size+1));
	  }

	  for (BinaryInt set=0; set<zeta_.cardinality; ++set){
	    ret[set] = number_;
  	  }
      return ret /= zeta_;
    }

}//namespace
