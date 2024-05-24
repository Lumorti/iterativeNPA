#include "printing.h"

// Pretty print complex numbers
std::ostream& operator<<(std::ostream& os, const std::complex<double>& c) {
    if (c.imag() == 0 && c.real() == 0) {
        os << "0";
    } else if (c.imag() == 0) {
        if (c.real() >= 0) {
            os << "+" << c.real();
        } else {
            os << c.real();
        }
    } else if (c.real() == 0) {
        if (c.imag() == 1) {
            os << "+i";
        } else if (c.imag() == -1) {
            os << "-i";
        } else {
            if (c.imag() > 0) {
                os << "+" << c.imag() << "i";
            } else {
                os << c.imag() << "i";
            }
        }
    } else {
        if (c.imag() >= 0) {
            os << "+(" << c.real() << "+" << c.imag() << "i)";
        } else {
            os << "+(" << c.real() << c.imag() << "i)";
        }
    }
    return os;
}

// When printing a matrix of doubles, print it as a string
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<double>>& m) {

    // Determine the maximum width of each column
    std::vector<int> columnWidths(m[0].size(), 0);
    for (size_t i=0; i<m.size(); i++) {
        for (size_t j=0; j<m[i].size(); j++) {
            std::stringstream ss;
            ss << m[i][j];
            int sizeWhenWritten = ss.str().size();
            columnWidths[j] = std::max(columnWidths[j], sizeWhenWritten);
        }
    }

    // Iterate through the matrix
    for (size_t i=0; i<m.size(); i++) {
        for (size_t j=0; j<m[i].size(); j++) {
            std::stringstream ss;
            ss << m[i][j];
            int sizeWhenWritten = ss.str().size();
            for (int k=0; k<columnWidths[j]-sizeWhenWritten; k++) {
                ss << " ";
            }
            ss << " ";
            os << ss.str();
        }
        os << std::endl;
    }

    // Return the stream
    return os;

}


