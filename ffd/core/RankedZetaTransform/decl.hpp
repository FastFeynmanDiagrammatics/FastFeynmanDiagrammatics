namespace ffd::ranked_zeta{

  template<typename Field = Real, SetT SetType = Full>
  class RankedZeta{
  private:
    std::vector<Field> coef;

  public:
    inline std::size_t size() const{
      std::size_t size = coef.size();
      assert(size == cardinality * (rank_size+1));
      return size;
  	}

	inline void resize(BinaryInt new_size) {
	  std::vector<Field> new_coef (new_size,0.);
	  coef = new_coef;
    }

  BinaryInt rank_size;
  BinaryInt cardinality;

  inline BinaryInt get_rank_size() const{return rank_size;}
  inline BinaryInt get_cardinality() const{return cardinality;}
  inline BinaryInt order() const{return ffd::set_theory::CardinalitySet(cardinality-1);}

  //"This is not the zero polynomial" (used to speed up arithmetic operations)
  bool NotZero = true;

  RankedZeta(): coef(1, 0.), rank_size(0), cardinality((1<<0)), NotZero(false) {}

  explicit RankedZeta(int order): coef((1<<order)*(order+1), 0.),rank_size(order), cardinality((1<<order)), NotZero(false) {}

  // explicit RankedZeta(int order, int rank): coef((1<<order)*order, 0.), NotZero(false), cardinality((1<<order)), rank_size(rank) {}

  // conversion.hpp
  RankedZeta(ffd::nilpotent_polynomial::NilpotentPolynomial<Field, SetType> const & P);
  //void operator= (ffd::nilpotent_polynomial::NilpotentPolynomial<Field, SetType> const & P) noexcept;

  RankedZeta(Field x): coef(1, 0.), rank_size(0), cardinality((1<<0)) {
      coef[0] = x;
  }

  //the [] operator uses the internal representation
  //whenever we touch the non-const one, the polynomial loses the eventual
  //"zero" status
  inline Field&       operator[](BinaryInt S){ NotZero = true; return coef[S];}
  inline Field const& operator[](BinaryInt S) const{ return coef[S];}

  //scalar_operations.hpp
                // RankedZeta  operator=  (const Field)       noexcept;
  [[nodiscard]] RankedZeta  operator+  ()            const noexcept;
  [[nodiscard]] RankedZeta  operator-  ()            const noexcept;
  [[nodiscard]] RankedZeta  operator+  (const Field) const noexcept;
  [[nodiscard]] RankedZeta  operator-  (const Field) const noexcept;
  [[nodiscard]] RankedZeta  operator*  (const Field) const noexcept;
  [[nodiscard]] RankedZeta  operator/  (const Field) const noexcept;
                RankedZeta& operator+= (const Field)       noexcept;
                RankedZeta& operator-= (const Field)       noexcept;
                RankedZeta& operator*= (const Field)       noexcept;
                RankedZeta& operator/= (const Field)       noexcept;

  //operator_sum.hpp
  [[nodiscard]] RankedZeta  operator+  (const RankedZeta&) const noexcept;
  [[nodiscard]] RankedZeta  operator-  (const RankedZeta&) const noexcept;
                RankedZeta& operator+= (const RankedZeta&)       noexcept;
                RankedZeta& operator-= (const RankedZeta&)       noexcept;

  //operator_comparison.hpp
  bool operator>  (const RankedZeta&) const noexcept;
  bool operator<  (const RankedZeta&) const noexcept;
  bool operator>= (const RankedZeta&) const noexcept;
  bool operator<= (const RankedZeta&) const noexcept;
  bool operator== (const RankedZeta&) const noexcept;
  bool operator>  (const Field)       const noexcept;
  bool operator<  (const Field)       const noexcept;
  bool operator>= (const Field)       const noexcept;
  bool operator<= (const Field)       const noexcept;
  bool operator== (const Field)       const noexcept;

  //operator_multiply.hpp
  [[nodiscard]] RankedZeta  operator*  (const RankedZeta&) const noexcept;
                RankedZeta& operator*= (const RankedZeta&)       noexcept;

  //operator_divide.hpp
  [[nodiscard]] RankedZeta  operator/  (const RankedZeta&) const noexcept;
                RankedZeta& operator/= (const RankedZeta&)       noexcept;

  //MaxAbsCoef.hpp
  //returns the maximal value of coef
  [[nodiscard]] Real MaxAbsCoef() noexcept; //why Real and not Field?

  //restrict.hpp
  [[nodiscard]] RankedZeta restrict_to_set  (BinaryInt set)  const noexcept;

  //extend.hpp
  [[nodiscard]] RankedZeta extend_to_set  (BinaryInt set, int order)  const noexcept;

  //shift.hpp
  [[nodiscard]] RankedZeta shift(BinaryInt set, int order) const noexcept;

};

// template<typename Field=Real>
// RankedZeta<Field, Even>::RankedZeta(RankedZeta<Field, Full> const & Z);
//
// template<typename Field=Real>
// RankedZeta<Field, Full>::RankedZeta(RankedZeta<Field, Even> const & Z);


}//namespace
