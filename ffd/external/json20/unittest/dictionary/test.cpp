#include "../../json20/include.hpp"

int main() {
  json20::dictionary jd, jd_2;
  jd["some value"] = 1;
  jd["copy"] = jd.at<int>("some value");
  jd["another value"] = 1.;
  jd["a vector"] = std::vector<double>{0.1, 0.2, 0.3, 0.3};
  std::vector<std::complex<double>> vec(10, 0.);
  vec[0] = std::complex<double>{1., 1.};
  jd["A complex vector"] = vec;
  jd_2["jd"] = jd;
  jd_2["some other stuff"] = "aaa";
  std::cerr << json20::to_string(jd_2);
  json20::to_file("file.json", jd_2);
}
