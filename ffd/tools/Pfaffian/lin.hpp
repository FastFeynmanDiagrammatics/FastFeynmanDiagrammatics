namespace ffd::gauss_pfaffian{

  inline int lin(int j, int k, int n){
    assert(( j < k )); assert(( j >= 0 )); assert(( n > 0 ));

    
    int const mjp1 = -(j+1);
    return (   (  ( (n<<2) + mjp1 )  *  j  )   >>   1   )    +
      k + mjp1;
  }


}//namespace
