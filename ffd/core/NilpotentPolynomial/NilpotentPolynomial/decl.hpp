namespace ffd::nilpotent_polynomial{

  template <typename Field = Real, SetT SetType = Full>
  class NilpotentPolynomial{
  public:
    //this is used just for printing
    static const char variable_names[26];

  private:
    //the array of coefficients
    std::vector<Field> coef;

  public:
    inline std::size_t size() const{ return coef.size(); }

    inline int linear_size() const{return  ffd::set_theory::CardinalitySet(coef.size()-1); }


    //"This is not the zero polynomial" (used to speed up arithmetic operations)
    bool NotZero = true;

    NilpotentPolynomial(): coef(1, 0.), NotZero(false) {}

    explicit NilpotentPolynomial(int degree_): coef(1<<degree_, 0.), NotZero(false) {}

    NilpotentPolynomial(Field x): coef(1, 0.) { coef[0] = x; }

    // conversion.hpp
    NilpotentPolynomial(ffd::ranked_zeta::RankedZeta<Field, SetType> const & Z);
	//void operator= (ffd::ranked_zeta::RankedZeta<Field, SetType> const & P) noexcept;

    //the [] operator uses the internal representation
    //whenever we touch the non-const one, the polynomial loses the eventual
    //"zero" status
    inline Field& operator[](BinaryInt S){ NotZero = true; return coef[S];}
    inline Field const& operator[](BinaryInt S) const{ return coef[S];}


    //scalar_operations.hpp
    NilpotentPolynomial& operator+=(const Field) noexcept;
    NilpotentPolynomial& operator-=(const Field) noexcept;
    [[nodiscard]] NilpotentPolynomial  operator+(const Field) const noexcept;
    [[nodiscard]] NilpotentPolynomial  operator-(const Field) const noexcept;
    NilpotentPolynomial& operator*=(const Field) noexcept;
    [[nodiscard]] NilpotentPolynomial  operator*(const Field) const noexcept;
    NilpotentPolynomial& operator/=(const Field) noexcept;
    [[nodiscard]] NilpotentPolynomial  operator/(const Field) const noexcept;
    [[nodiscard]] NilpotentPolynomial  operator-() const noexcept;
    [[nodiscard]] NilpotentPolynomial  operator+() const noexcept;

    //operator_sum.hpp
    NilpotentPolynomial& operator+=(const NilpotentPolynomial&) noexcept;
    NilpotentPolynomial& operator-=(const NilpotentPolynomial&) noexcept;
    [[nodiscard]] NilpotentPolynomial  operator+(const NilpotentPolynomial&) const noexcept;
    [[nodiscard]] NilpotentPolynomial  operator-(const NilpotentPolynomial&) const noexcept;

	//operator_comparison.hpp
	bool operator>(const NilpotentPolynomial&) const noexcept;
	bool operator<(const NilpotentPolynomial&) const noexcept;
	bool operator>=(const NilpotentPolynomial&) const noexcept;
	bool operator<=(const NilpotentPolynomial&) const noexcept;
	bool operator==(const NilpotentPolynomial&) const noexcept;
	bool operator>(const Field) const noexcept;
	bool operator<(const Field) const noexcept;
	bool operator>=(const Field) const noexcept;
	bool operator<=(const Field) const noexcept;

    //operator_multiply.hpp
    [[nodiscard]] NilpotentPolynomial  operator*(const NilpotentPolynomial&) const noexcept;
    NilpotentPolynomial& operator*=(const NilpotentPolynomial&) noexcept;


    //operator_divide.hpp
    [[nodiscard]] NilpotentPolynomial  operator/(const NilpotentPolynomial&) const noexcept;
    NilpotentPolynomial& operator/=(const NilpotentPolynomial&) noexcept;


    //MaxAbsCoef.hpp
    //returns the maximal value of coef
    [[nodiscard]] Real MaxAbsCoef() noexcept;

    //restrict_to_set.hpp
    [[nodiscard]] NilpotentPolynomial
    restrict_to_set(BinaryInt set) const noexcept;

    //extend_to_set.hpp
    [[nodiscard]] NilpotentPolynomial
    extend_to_set(BinaryInt set, int order) const noexcept;

    //shift.hpp
    [[nodiscard]] NilpotentPolynomial
    shift(BinaryInt set, int order) const noexcept;

  };

  template <typename Field, SetT SetType>
  const char NilpotentPolynomial<Field, SetType>::variable_names[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};



}//namespace ffd::nilpotent_polynomial
