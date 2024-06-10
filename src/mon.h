#pragma once
#include <vector>
#include <complex>
#include <string>
#include <set>

// Non-commuting monomial class
class Mon {
public:

    // Vars
    std::vector<std::pair<char, int>> monomial;

    // If initialized from nothing
    Mon();

    // If initialized from a vector
    Mon(std::vector<std::pair<char, int>> monomial);

    // If initialized from a char and a number
    Mon(char variable, int index);

    // If initialized from a single part of monomial
    Mon(std::pair<char, int> part);

    // If initialized from a string
    Mon(std::string asString);

    // If asked to reverse the monomial in place
    void reverse();
    
    // If asked to reverse the monomial
    Mon reversed() const;

    // If multiplied by another monomial
    Mon operator*(const Mon& other) const;

    // Bracket access
    std::pair<char, int>& operator[](size_t index);
    const std::pair<char, int>& operator[](size_t index) const;

    // Iterators
    std::vector<std::pair<char,int>>::iterator begin() {return monomial.begin();}
    std::vector<std::pair<char,int>>::iterator end()   {return monomial.end();}
    std::vector<std::pair<char,int>>::const_iterator begin() const {return monomial.begin();}
    std::vector<std::pair<char,int>>::const_iterator end() const {return monomial.end();}

    // Check equality
    bool operator==(const Mon& other) const;
    bool operator!=(const Mon& other) const;
    bool operator==(const int other) const;
    bool operator!=(const int other) const;

    // Size of the monomial
    size_t size() const;

    // Compare two sections of a monomial, but reversed
    bool compareReversed(const std::pair<char, int>& a, const std::pair<char, int>& b) const;

    // Comparison between this and another monomial, but reversed
    bool compareReversed(Mon& a, Mon& b) const;

    // Comparison between this and another monomial
    bool operator<(const Mon& other) const;
    bool operator>(const Mon& other) const;

    // Given a monomial, reduce it to its simplest form
    std::pair<std::complex<double>, Mon> reduce() const;

    // Pretty printing
    friend std::ostream& operator<<(std::ostream& os, const Mon& m);
    friend std::ostream& operator<<(std::ostream& os, const std::pair<char, int>& p);
    friend std::ostream& operator<<(std::ostream& os, const std::vector<Mon>& m);
    friend std::ostream& operator<<(std::ostream& os, const std::set<Mon>& m);

    // Convert a monomial to a string
    std::string toString() const;

    // Cycle the monomial to put a certain thing at the end
    Mon cycleTo(char variable, int index) const;

    // Check if the monomial is constant
    bool isConstant() const;

    // Check if contains a specific letter
    bool contains(char letter) const;

};

