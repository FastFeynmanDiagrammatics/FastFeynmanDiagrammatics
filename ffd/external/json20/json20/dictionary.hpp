struct dictionary {
  std::map<std::string, std::any> data;

  template <class T>
  auto at(std::string s) const {
    return std::any_cast<T>(data.at(s));
  }

  auto& operator[](std::string s) { return data[s]; }
};
