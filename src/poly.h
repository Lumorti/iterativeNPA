#pragma once
#include <map>
#include "mon.h"

// Non-commuting polynomial class
class Poly {
public:

    // Vars
    std::map<Mon, std::complex<double>> polynomial;
    int verbosity = 1;

    // If initialized from nothing
    Poly();

    // If initialized from just a coeff
    Poly(std::complex<double> coeff);

    // If initialized from just a monomial
    Poly(Mon mon);

    // If initialized from a vector
    Poly (std::map<Mon, std::complex<double>> poly);

    // If initialized from a single coeff + monomial
    Poly(std::pair<std::complex<double>, Mon> pair);

    // If initialized from a coeff and a monomial
    Poly(std::complex<double> coeff, Mon mon);

    // If initialized from a single part of monomial
    Poly(std::pair<char, int> part);

    // If initialized from a coeff and a string
    Poly(std::complex<double> coeff, std::string asString);

    // If initialized from a string
    Poly(std::string asString);

    // When assigning manually with a vector
    Poly& operator=(const std::map<Mon, std::complex<double>>& poly);

    // Checking equality of two polynomials
    bool operator==(const Poly& other);
    bool operator!=(const Poly& other);

    // When summing two polynomials
    Poly operator+(const Poly& other);

    // When summing in-place
    Poly& operator+=(const Poly& other);

    // When subtracting two polynomials
    Poly operator-(const Poly& other);

    // When subtracting in-place
    Poly& operator-=(const Poly& other);

    // When multiplying two polynomials
    Poly operator*(const Poly& other) const;

    // When multiplying in-place
    Poly& operator*=(const Poly& other);

    // When multiplying by a constant
    Poly operator*(const std::complex<double>& other);
    friend Poly operator*(const std::complex<double>& other, const Poly& p);

    // When multiplying by a monomial
    Poly operator*(const Mon& other);
    friend Poly operator*(const Mon& other, const Poly& p);

    // When multiplying by a constant in-place
    Poly& operator*=(const std::complex<double>& other);

    // When dividing by a constant
    Poly operator/(const std::complex<double>& other);

    // When dividing by a constant in-place
    Poly& operator/=(const std::complex<double>& other);

    // Size of the polynomial
    const size_t size() const;

    // Allow bracket access
    const std::complex<double> operator[](Mon mon) const;
    std::complex<double>& operator[](Mon mon);

    // Allow getting a random term
    std::complex<double>& getValue();
    Mon getKey();

    // Iterators
    std::map<Mon,std::complex<double>>::iterator begin() {return polynomial.begin();}
    std::map<Mon,std::complex<double>>::iterator end()   {return polynomial.end();}
    std::map<Mon,std::complex<double>>::const_iterator begin() const {return polynomial.begin();}
    std::map<Mon,std::complex<double>>::const_iterator end() const {return polynomial.end();}

    // Self-negative
    Poly operator-();

    // Self-conjugate
    Poly conj() const;

    // Dagger operator
    Poly dagger() const;

    // Evaluate, given a list of variables and values
    std::complex<double> eval(std::vector<std::pair<Mon, std::complex<double>>>& vals);
    std::complex<double> eval(std::map<Mon, std::complex<double>>& vals);

    // Given a reduce each monomial and combine
    void reduce();
    Poly reduced();

    // Remove any zero terms from this polynomial
    void clean();

    // Pretty print
    friend std::ostream& operator<<(std::ostream& os, const Poly& p);
    friend std::ostream& operator<<(std::ostream& os, const std::vector<Poly>& v);
    friend std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<Poly>>& m);

    // Replace a monomial by a polynomial
    void replace(std::pair<char,int> mon, Poly replacement);
    Poly replaced(std::pair<char,int> mon, Poly replacement);
    void replace(Mon mon, Poly replacement);
    Poly replaced(Mon mon, Poly replacement);

    // Reduce a polynomial using plus/minus operators to Paulis
    void convertToPaulis();

    // Convert a polynomial to Paulis
    Poly convertedToPaulis();

    // Comuttator and anticommutator
    Poly commutator(Poly other);
    Poly anticommutator(Poly other);

    // Get the list of monomials
    std::vector<Mon> monomials();

    // Get the polynomial rearranged to equal a certain monomial
    Poly rearranged(Mon mon);

    // Check if a monomial is in the polynomial
    bool contains(const Mon mon) const;

    // Cycle each term such that a monomial is at the end
    void cycleTo(char variable, int index);
    void cycleToAndRemove(char variable, int index);

    // Check if the polynomial is constant
    bool isConstant() const;

    // Apply a monomial map
    Poly applyMap(std::map<Mon, Mon> map);

};

