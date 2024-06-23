void to_file(std::string filename, auto a) {
  std::ofstream ff(filename, std::ios_base::trunc);
  ff << to_string(a);
}
