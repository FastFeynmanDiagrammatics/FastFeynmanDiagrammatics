#include <ffd/std.hpp>

#ifdef FFD_SIMPLE_FFT_FLAG
#include "ffd/external/Simple-FFT/include/simple_fft/fft.hpp"
#endif

#include <ffd/core.hpp>
#include <ffd/tools.hpp>
#include <ffd/packages.hpp>

#include "include.hpp"
#include "call.hpp"

int main() {
  ffd::fft_matsubara_3d::unit_test::UnitTest();
}
