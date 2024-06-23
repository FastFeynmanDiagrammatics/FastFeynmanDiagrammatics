namespace ffd::distributions::unit_test{

void simple_test(){
    std::vector<Real> v(2,0.);
    auto f1 = gaussian(2.);
    auto f2 = lorentz(2.);
    auto f3 = pade_gauss(2.);
    assert(std::abs(f1(0)-1.)<1e-10);
    assert(std::abs(f2(0)-1.)<1e-10);
    assert(std::abs(f3(0)-1.)<1e-10);
    assert(std::abs(f1(v)-1.)<1e-10);
    assert(std::abs(f2(v)-1.)<1e-10);
    assert(std::abs(f3(v)-1.)<1e-10);
}

}//namespace
