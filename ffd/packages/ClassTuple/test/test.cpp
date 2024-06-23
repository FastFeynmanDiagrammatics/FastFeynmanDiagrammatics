#include<ffd/std.hpp>
#include<ffd/core.hpp>
#include"../../PackagesMath.hpp"
#include"../../ClassTuple.hpp"


using namespace ffd::class_tuple;


template<typename T, int d>
struct test{
  T x;
  int l = d;
};


int main(){
  ClassTupleTypeInt<test, 2, double, int, int> X;
  ffd::Get<0>(X).x = 3.14;
  ClassTupleTypeInt<test, 2, double, int, int> Y;
  ffd::Get<0>(Y).x = -3.14;
  auto Z = ffd::Merge(X, Y);
  std::cout<<ffd::Get<0>(Z).x<<" "<<ffd::Get<ffd::Size(X)>(Z).x<<std::endl;
}
