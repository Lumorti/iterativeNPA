#pragma once
#include <vector>
#include <complex>
#include <iostream>

// Pretty print complex numbers
std::ostream& operator<<(std::ostream& os, const std::complex<double>& c);

// When printing a matrix of doubles, print it as a string
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<double>>& m);

// When printing a vector of something
template <typename T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {

    // Iterate through the vector
    os << "{";
    for (size_t i=0; i<v.size(); i++) {
        os << v[i];
        if (i < v.size()-1) {
            os << ", ";
        }
    }
    os << "}";

    // Return the stream
    return os;

}

