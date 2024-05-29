#include <algorithm>
#include <sstream>
#include <iostream>
#include "mon.h"
#include "settings.h"

// Useful constants
const std::complex<double> imag(0, 1);
const double zeroTol = 1e-13;

// If initialized from nothing
Mon::Mon() {
    monomial = {};
}

// If initialized from a vector
Mon::Mon(std::vector<std::pair<char, int>> monomial) {
    monomial = monomial;
}

// If initialized from a char and a number
Mon::Mon(char variable, int index) {
    monomial = {{variable, index}};
}

// If initialized from a single part of monomial
Mon::Mon(std::pair<char, int> part) {
    monomial = {part};
}

// If initialized from a string
Mon::Mon(std::string asString) {

    // If the string is empty, return an empty monomial
    if (asString == "" || asString.find("<") == std::string::npos) {
        return;
    }

    // Iterate through the string
    std::vector<std::pair<char, int>> toReturn;
    char currentVariable = ' ';
    std::string currentIndex;
    for (size_t i=0; i<asString.size(); i++) {

        // Check if this char is an integer
        if ((asString[i] >= '0' && asString[i] <= '9') || asString[i] == '-') {
            currentIndex += asString[i];

        // Otherwise it's a variable
        } else if (asString[i] >= 'A' && asString[i] <= 'Z') {
            currentIndex += asString[i];
            if (currentVariable != ' ') {
                toReturn.push_back(std::make_pair(currentVariable, std::stoi(currentIndex)));
            }
            currentVariable = asString[i];
            currentIndex = "";
        }

    }

    // Add the last variable
    toReturn.push_back(std::make_pair(currentVariable, std::stoi(currentIndex)));

    // Set the monomial
    monomial = toReturn;

}

// If asked to reverse the monomial in place
void Mon::reverse() {
    std::reverse(monomial.begin(), monomial.end());
}

// If asked to reverse the monomial
Mon Mon::reversed() const {
    Mon toReturn = *this;
    std::reverse(toReturn.monomial.begin(), toReturn.monomial.end());
    return toReturn;
}

// Multiply two monomials
Mon Mon::operator*(const Mon& other) const {
    Mon toReturn;
    toReturn.monomial = monomial;
    toReturn.monomial.insert(toReturn.monomial.end(), other.monomial.begin(), other.monomial.end());
    return toReturn;
}

// Bracket access
std::pair<char, int>& Mon::operator[](size_t index) {
    return monomial[index];
}
const std::pair<char, int>& Mon::operator[](size_t index) const {
    return monomial[index];
}

// Check equality
bool Mon::operator==(const Mon& other) const {
    return monomial == other.monomial;
}
bool Mon::operator!=(const Mon& other) const {
    return monomial != other.monomial;
}
bool Mon::operator==(const int other) const {
    return monomial.size() == 0;
}
bool Mon::operator!=(const int other) const {
    return monomial.size() != 0;
}

// Size of the monomial
size_t Mon::size() const {
    return monomial.size();
}

// Compare two sections of a monomial, but reversed
bool Mon::compareReversed(const std::pair<char, int>& a, const std::pair<char, int>& b) const {
    if (a.second == b.second) {
        return a.first < b.first;
    } else {
        return a.second < b.second;
    }
}

// Comparison between this and another monomial, but reversed
bool Mon::compareReversed(Mon& a, Mon& b) const {
    for (size_t i=0; i<std::min(a.size(), b.size()); i++) {
        if (a[i].second != b[i].second) {
            return a[i].second < b[i].second;
        }
    }
    return a.size() < b.size();
}

// Comparison between this and another monomial
bool Mon::operator<(const Mon& other) const {
    if (monomial.size() != other.monomial.size()) {
        return monomial.size() < other.monomial.size();
    } else {
        for (size_t i=0; i<monomial.size(); i++) {
            if (sortOrder == "numberFirst") {
                if (monomial[i].second != other.monomial[i].second) {
                    return monomial[i].second < other.monomial[i].second;
                } else if (monomial[i].first != other.monomial[i].first) {
                    return monomial[i].first < other.monomial[i].first;
                }
            } else if (sortOrder == "letterFirst") {
                if (monomial[i].first != other.monomial[i].first) {
                    return monomial[i].first < other.monomial[i].first;
                } else if (monomial[i].second != other.monomial[i].second) {
                    return monomial[i].second < other.monomial[i].second;
                }
            }
        }
    }
    return false;
}
bool Mon::operator>(const Mon& other) const {
    if (monomial.size() != other.monomial.size()) {
        return monomial.size() > other.monomial.size();
    } else {
        for (size_t i=0; i<monomial.size(); i++) {
            if (sortOrder == "numberFirst") {
                if (monomial[i].second != other.monomial[i].second) {
                    return monomial[i].second > other.monomial[i].second;
                } else if (monomial[i].first != other.monomial[i].first) {
                    return monomial[i].first > other.monomial[i].first;
                }
            } else if (sortOrder == "letterFirst") {
                if (monomial[i].first != other.monomial[i].first) {
                    return monomial[i].first > other.monomial[i].first;
                } else if (monomial[i].second != other.monomial[i].second) {
                    return monomial[i].second > other.monomial[i].second;
                }
            }
        }
    }
    return false;
}

// Reduce a monomial, returning the coefficient and the reduced monomial
std::pair<std::complex<double>, Mon> Mon::reduce() const {

    // Sort the monomial as much as we can
    Mon mon = *this;
    int monSize = int(mon.size());
    int monSize1 = monSize-1;
    if (commutating == "letters") {
        for (int i=0; i<monSize; i++) {
            for (int j=0; j<monSize1; j++) {
                if (mon[j].first != mon[j+1].first && mon[j] > mon[j+1]) {
                    std::swap(mon[j], mon[j+1]);
                }
            }
        }
    } else if (commutating == "numbers") {
        for (int i=0; i<monSize; i++) {
            for (int j=0; j<monSize1; j++) {
                if (mon[j].second != mon[j+1].second && mon[j].second > 0 && mon[j+1].second > 0 && !compareReversed(mon[j], mon[j+1])) {
                    std::swap(mon[j], mon[j+1]);
                }
            }
        }
    } else if (commutating == "both") {
        for (int i=0; i<monSize; i++) {
            for (int j=0; j<monSize1; j++) {
                if (mon[j] > mon[j+1]) {
                    std::swap(mon[j], mon[j+1]);
                }
            }
        }
    }
    
    // Simplify Pauli strings using the following rules (plus conjugates):
    // XY = iZ
    // ZX = iY
    // YZ = iX
    std::complex<double> coeff(1, 0);
    if (pauliReductions) {
        for (int i=mon.size()-1; i>0; i--) {
            if (mon[i-1].second == mon[i].second) {
                if (mon[i-1].first == 'X' && mon[i].first == 'Y') {
                    coeff *= imag;
                    mon[i-1] = std::make_pair('Z', mon[i-1].second);
                    mon.monomial.erase(mon.monomial.begin()+i);
                } else if (mon[i-1].first == 'X' && mon[i].first == 'Z') {
                    coeff *= -imag;
                    mon[i-1] = std::make_pair('Y', mon[i-1].second);
                    mon.monomial.erase(mon.monomial.begin()+i);
                } else if (mon[i-1].first == 'Y' && mon[i].first == 'Z') {
                    coeff *= imag;
                    mon[i-1] = std::make_pair('X', mon[i-1].second);
                    mon.monomial.erase(mon.monomial.begin()+i);
                } else if (mon[i-1].first == 'Y' && mon[i].first == 'X') {
                    coeff *= -imag;
                    mon[i-1] = std::make_pair('Z', mon[i-1].second);
                    mon.monomial.erase(mon.monomial.begin()+i);
                } else if (mon[i-1].first == 'Z' && mon[i].first == 'X') {
                    coeff *= imag;
                    mon[i-1] = std::make_pair('Y', mon[i-1].second);
                    mon.monomial.erase(mon.monomial.begin()+i);
                } else if (mon[i-1].first == 'Z' && mon[i].first == 'Y') {
                    coeff *= -imag;
                    mon[i-1] = std::make_pair('X', mon[i-1].second);
                    mon.monomial.erase(mon.monomial.begin()+i);
                } else if (mon[i-1].first == 'X' && mon[i].first == 'X') {
                    mon.monomial.erase(mon.monomial.begin()+i);
                    mon.monomial.erase(mon.monomial.begin()+i-1);
                    i--;
                } else if (mon[i-1].first == 'Y' && mon[i].first == 'Y') {
                    mon.monomial.erase(mon.monomial.begin()+i);
                    mon.monomial.erase(mon.monomial.begin()+i-1);
                    i--;
                } else if (mon[i-1].first == 'Z' && mon[i].first == 'Z') {
                    mon.monomial.erase(mon.monomial.begin()+i);
                    mon.monomial.erase(mon.monomial.begin()+i-1);
                    i--;
                }
            }
        }
    }

    // <A1A1> = <1>
    if (ASquaredEquals == "1") {
        int i = 0;
        while (i < int(mon.size())-1) {
            if (mon[i] == mon[i+1]) {
                mon.monomial.erase(mon.monomial.begin()+i+1);
                mon.monomial.erase(mon.monomial.begin()+i);
                i -= 2;
            }
            i++;
        }
    } else if (ASquaredEquals == "A") {
        int i = 0;
        while (i < int(mon.size())-1) {
            if (mon[i] == mon[i+1]) {
                mon.monomial.erase(mon.monomial.begin()+i+1);
                i--;
            }
            i++;
        }
    }

    // Flip it to see if it's smaller
    if (trySwap) {
        Mon monFlipped = mon.reversed();
        if (monFlipped < mon) {
            mon = monFlipped;
        }
    }

    return {coeff, mon};

}

// Pretty printing
std::ostream& operator<<(std::ostream& os, const Mon& m) {

    // Iterate through the monomial
    if (m.size() > 0) {
        os << "<";
        for (size_t i=0; i<m.size(); i++) {
            os << m.monomial[i].first << m.monomial[i].second;
        }
        os << ">";
    } else {
        os << "1";
    }

    // Return the stream
    return os;

}

// Convert a monomial to a string
std::string Mon::toString() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

// Pretty print part of monomial
std::ostream& operator<<(std::ostream& os, const std::pair<char, int>& p) {
    os << p.first << p.second;
    return os;
}

// Pretty print vector of monomials
std::ostream& operator<<(std::ostream& os, const std::vector<Mon>& m) {
    os << "{";
    for (size_t i=0; i<m.size(); i++) {
        os << m[i];
        if (i < m.size()-1) {
            os << ", ";
        }
    }
    os << "}";
    return os;
}

// Pretty print set of monomials
std::ostream& operator<<(std::ostream& os, const std::set<Mon>& m) {
    os << "{";
    for (auto it=m.begin(); it!=m.end(); it++) {
        os << *it;
        if (std::next(it) != m.end()) {
            os << ", ";
        }
    }
    os << "}";
    return os;
}

// Cycle the monomial to put a certain thing at the end
Mon Mon::cycleTo(char variable, int index) const {

    // Find this variable
    int location = -1;
    for (size_t i=0; i<monomial.size(); i++) {
        if (monomial[i].first == variable && monomial[i].second == index) {
            location = i;
            break;
        }
    }

    // Cycle by this amount
    Mon toReturn = *this;
    std::rotate(toReturn.monomial.begin(), toReturn.monomial.begin()+location+1, toReturn.monomial.end());
    return toReturn;

}

// Check if the monomial is empty (and thus represents 1)
bool Mon::isConstant() const {
    return size() == 0;
}


    
