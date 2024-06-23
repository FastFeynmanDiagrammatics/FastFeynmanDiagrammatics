namespace ffd::periodic_coordinate::unit_test{

	void operator_plus_test(){

      using ffd::sign, ffd::get;
      using namespace ffd::user_space;

      PeriodicCoordinates<2, double, int, int> X;
      auto Y1 = X + std::make_tuple(0.5,1,-1);

	  auto Y2 = X;
  	  component<0>(Y2) = component<0>(X) + 0.5;
	  component<1>(Y2) = component<1>(X) + 1;
      component<2>(Y2) = component<2>(X) - 1;

      assert(component<0>(Y1) == component<0>(Y2));
      assert(component<1>(Y1) == component<1>(Y2));
      assert(component<2>(Y1) == component<2>(Y2));

      auto X1 = X;
	  X += std::make_tuple(0.5,1,-1);
	  X = X + std::make_tuple(0.3,2,-1);
	  X -= std::make_tuple(0.4,3,-2);
	  X = X - std::make_tuple(0.4,0,0);

	  assert(std::abs(component<0>(X) - component<0>(X1))<1e-10);
	  assert(component<1>(X) == component<1>(X1));
	  assert(component<2>(X) == component<2>(X1));

    }

}//namespace
