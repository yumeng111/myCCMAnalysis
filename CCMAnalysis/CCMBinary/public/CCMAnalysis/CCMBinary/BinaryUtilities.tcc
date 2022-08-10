#ifndef CCMAnalysis_BinaryUtilities_TCC
#define CCMAnalysis_BinaryUtilities_TCC

#include <string>
#include <vector>
#include <iostream>

template<typename T>
std::ostream & write_binary(std::ostream & os, std::vector<T> const & v) {
    constexpr size_t size_size = sizeof(uint64_t);
    uint64_t vector_size = v.size();
    os.write((char*)&vector_size, size_size);
    for(size_t i=0; i<vector_size; ++i) {
        write_binary(os, v[i]);
    }
    return os;
}

template<typename T>
std::istream & read_binary(std::istream & is, std::vector<T> & v) {
    constexpr size_t size_size = sizeof(uint64_t);
    uint64_t vector_size;
    is.read((char*)&vector_size, size_size);
    v.resize(vector_size);
    for(size_t i=0; i<vector_size; ++i) {
        read_binary(is, v[i]);
    }
    return is;
}

template<typename T>
std::ostream & write_binary_contiguous_vector(std::ostream & os, std::vector<T> const & v) {
    constexpr size_t size_size = sizeof(uint64_t);
    constexpr size_t T_size = sizeof(T);
    uint64_t vector_size = v.size();
    os.write((char*)&vector_size, size_size);
    os.write((char*)v.data(), vector_size * T_size);
    return os;
}

template<typename T>
std::istream & read_binary_contiguous_vector(std::istream & is, std::vector<T> & v) {
    constexpr size_t size_size = sizeof(uint64_t);
    constexpr size_t T_size = sizeof(T);
    uint64_t vector_size;
    is.read((char*)&vector_size, size_size);
    v.resize(vector_size);
    is.read((char*)v.data(), vector_size * T_size);
    return is;
}

template<>
std::ostream & write_binary(std::ostream & os, std::vector<uint32_t> const & v) {
    return write_binary_contiguous_vector<uint32_t>(os, v);
}
template<>
std::istream & read_binary(std::istream & is, std::vector<uint32_t> & v) {
    return read_binary_contiguous_vector<uint32_t>(is, v);
}

template<>
std::ostream & write_binary(std::ostream & os, std::vector<uint16_t> const & v) {
    return write_binary_contiguous_vector<uint16_t>(os, v);
}
template<>
std::istream & read_binary(std::istream & is, std::vector<uint16_t> & v) {
    return read_binary_contiguous_vector<uint16_t>(is, v);
}

template<>
std::ostream & write_binary(std::ostream & os, std::vector<float> const & v) {
    return write_binary_contiguous_vector<float>(os, v);
}
template<>
std::istream & read_binary(std::istream & is, std::vector<float> & v) {
    return read_binary_contiguous_vector<float>(is, v);
}

#endif // CCMAnalysis_BinaryUtilities_TCC
