struct file {
  std::stringstream ss;
  bool is_empty = true;
  std::string filename = "";

  file() { ss << "{\n"; }

  file(std::string fname) : filename(fname) { ss << "{\n"; }

  void write(std::string name, std::floating_point auto x) {
    if (!is_empty)
      ss << ",\n";
    ss << "\t\"" << name << "\": ";
    ss << x;
    is_empty = false;
  }

  void write(std::string name, std::integral auto x) {
    if (!is_empty)
      ss << ",\n";
    ss << "\t\"" << name << "\": ";
    ss << std::setprecision(20) << x;
    is_empty = false;
  }

  void write(std::string name, std::string x) {
    if (!is_empty)
      ss << ",\n";
    ss << "\t\"" << name << "\": ";
    ss << "\"" << x << "\"";
    is_empty = false;
  }

  template <class F>
  void write_vector(std::string name, std::vector<F> x) {
    assert((size(x) > 0));
    if (!is_empty)
      ss << ",\n";
    ss << "\t\"" << name << "\": [\n";
    if constexpr (std::is_same_v<F, std::string>) {
      for (long j = 0; j < long(x.size()) - 1; ++j) {
        ss << "\t\t\"" << x[j] << "\",\n";
      }  // for j in range(1, x.size()-1)
      ss << "\t\t\"" << x << "\"]";
    } else {
      for (long j = 0; j < long(x.size()) - 1; ++j) {
        ss << std::setprecision(20) << "\t\t" << x[j] << ",\n";
      }  // for j in range(1, x.size()-1)
      ss << std::setprecision(20) << "\t\t" << x[size(x) - 1] << "]";
    }
    is_empty = false;
  }

  void to_file(std::string fname) {
    std::string filename_json = fname;
    if (!filename_json.ends_with(".json"))
      filename_json += ".json";
    std::ofstream file_out(filename_json, std::ios::trunc);
    file_out << ss.str();
    file_out << "\n}\n";
  }

  ~file() {
    if (filename != "")
      this->to_file(filename);
  }
};
