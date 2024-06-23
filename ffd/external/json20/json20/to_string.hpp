std::string to_string(dictionary);

std::string to_string(std::complex<double> z) {
  std::stringstream ss;
  ss << std::setprecision(16) << "[" << real(z) << ", " << imag(z) << "]";
  return ss.str();
}

std::string to_string(std::string s) {
  std::stringstream ss;
  if (!s.starts_with("\""))
    ss << "\"";
  ss << s;
  if (!s.ends_with("\""))
    ss << "\"";
  return ss.str();
}

std::string to_string(std::any a) {
  const char* char_a;
  std::stringstream ss;
  ss << std::setprecision(16);
  if (a.type() == typeid(dictionary{})) {
    ss << to_string(std::any_cast<dictionary>(a));
  } else if (a.type() == typeid(char_a)) {
    std::string bb = std::any_cast<const char*>(a);
    ss << to_string(bb);
  } else if (a.type() == typeid(int(0))) {
    ss << std::any_cast<int>(a);
  } else if (a.type() == typeid(long(0))) {
    ss << std::any_cast<long>(a);
  } else if (a.type() == typeid(0ul)) {
    ss << std::any_cast<unsigned long>(a);
  } else if (a.type() == typeid(0.)) {
    ss << std::setprecision(16) << std::any_cast<double>(a);
  } else if (a.type() == typeid(std::complex<double>(0., 0.))) {
    auto val = std::any_cast<std::complex<double>>(a);
    ss << "[" << std::setprecision(16) << real(val) << ", " << imag(val) << "]";
  } else if (a.type() == typeid('a')) {
    ss << "\"" << std::any_cast<char>(a) << "\"";
  } else if (a.type() == typeid(std::string(""))) {
    ss << to_string(std::any_cast<std::string>(a));
  } else if (a.type() == typeid(std::array<double, 1>{})) {
    auto vec = std::any_cast<std::array<double, 1>>(a);
    ss << "[";
    for (long j = 0; j < size(vec) - 1l; ++j) {
      ss << vec[j] << ", ";
    }  // for j in range(0, size(vec))
    if (size(vec) != 0)
      ss << vec[size(vec) - 1] << "]";
  } else if (a.type() == typeid(std::array<double, 2>{})) {
    auto vec = std::any_cast<std::array<double, 2>>(a);
    ss << "[";
    for (long j = 0; j < size(vec) - 1; ++j) {
      ss << vec[j] << ", ";
    }  // for j in range(0, size(vec))
    if (size(vec) != 0)
      ss << vec[size(vec) - 1] << "]";
  } else if (a.type() == typeid(std::array<double, 3>{})) {
    auto vec = std::any_cast<std::array<double, 3>>(a);
    ss << "[";
    for (long j = 0; j < size(vec) - 1; ++j) {
      ss << vec[j] << ", ";
    }  // for j in range(0, size(vec))
    if (size(vec) != 0)
      ss << vec[size(vec) - 1] << "]";
  } else if (a.type() == typeid(std::array<double, 4>{})) {
    auto vec = std::any_cast<std::array<double, 4>>(a);
    ss << "[";
    for (long j = 0; j < size(vec) - 1; ++j) {
      ss << vec[j] << ", ";
    }  // for j in range(0, size(vec))
    if (size(vec) != 0)
      ss << vec[size(vec) - 1] << "]";
  } else if (a.type() == typeid(std::array<double, 5>{})) {
    auto vec = std::any_cast<std::array<double, 5>>(a);
    ss << "[";
    for (long j = 0; j < size(vec) - 1; ++j) {
      ss << vec[j] << ", ";
    }  // for j in range(0, size(vec))
    if (size(vec) != 0)
      ss << vec[size(vec) - 1] << "]";
  } else if (a.type() == typeid(std::vector<int>{})) {
    auto vec = std::any_cast<std::vector<int>>(a);
    ss << "[";
    for (long j = 0; j < size(vec) - 1; ++j) {
      ss << vec[j] << ", ";
    }  // for j in range(0, size(vec))
    if (size(vec) != 0)
      ss << vec[size(vec) - 1] << "]";
  } else if (a.type() == typeid(std::vector<long>{})) {
    auto vec = std::any_cast<std::vector<long>>(a);
    ss << "[";
    for (long j = 0; j < size(vec) - 1; ++j) {
      ss << vec[j] << ", ";
    }  // for j in range(0, size(vec))
    if (size(vec) != 0)
      ss << vec[size(vec) - 1] << "]";
  } else if (a.type() == typeid(std::vector<unsigned long>{})) {
    auto vec = std::any_cast<std::vector<unsigned long>>(a);
    ss << "[";
    for (long j = 0; j < size(vec) - 1; ++j) {
      ss << vec[j] << ", ";
    }  // for j in range(0, size(vec))
    if (size(vec) != 0)
      ss << vec[size(vec) - 1] << "]";
  } else if (a.type() == typeid(std::vector<double>{})) {
    auto vec = std::any_cast<std::vector<double>>(a);
    ss << std::setprecision(16) << "[";
    for (long j = 0; j < size(vec) - 1; ++j) {
      ss << vec[j] << ", ";
    }  // for j in range(0, size(vec))
    if (size(vec) != 0)
      ss << vec[size(vec) - 1] << "]";
  } else if (a.type() == typeid(std::vector<std::complex<double>>{})) {
    auto vec = std::any_cast<std::vector<std::complex<double>>>(a);
    ss << std::setprecision(16) << "[";
    for (long j = 0; j < size(vec) - 1; ++j) {
      ss << to_string(vec[j]) << ", ";
    }  // for j in range(0, size(vec))
    if (size(vec) != 0)
      ss << to_string(vec[size(vec) - 1]) << "]";
  }
  return ss.str();
}

std::string to_string(dictionary dic) {
  std::stringstream ss;
  ss << "{\n";
  long counter = 0, dic_size = dic.data.size();
  for (auto const& [key, value] : dic.data) {
    ss << "\t\"" << key << "\": ";
    ss << to_string(value);
    if (counter != dic_size - 1) {
      ss << ",";
    }
    ss << "\n";
    ++counter;
  }
  ss << "}";
  return ss.str();
}
