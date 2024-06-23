namespace ffd::user_space{

  template<typename T, SetT SetType>
    using RankedZeta = ffd::ranked_zeta::RankedZeta<T, SetType>;

    using RZeta = ffd::ranked_zeta::RankedZeta<Real, Full>;

    using FullZeta = ffd::ranked_zeta::RankedZeta<Real, Full>;

    using EvenZeta = ffd::ranked_zeta::RankedZeta<Real, Even>;

    using OddZeta = ffd::ranked_zeta::RankedZeta<Real, Odd>;

}//namespace ffd::user_space
