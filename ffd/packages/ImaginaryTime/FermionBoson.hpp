#pragma once

namespace ffd::imaginary_time{

auto Fermionic(decltype(ffd::user_space::ImaginaryTimes(1.)) const & ImagTimes_){
    return ImagTimes_.second;
}
auto Bosonic(decltype(ffd::user_space::ImaginaryTimes(1.)) const & ImagTimes_){
    return ImagTimes_.first;
}

}//namespace
