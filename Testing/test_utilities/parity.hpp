//******************************************************************************
///  @file parity.hpp
///  @brief Contains classes to control sine and cosine parity of fourier
///  quantites.
//******************************************************************************

#ifndef parity_hpp
#define parity_hpp

#include <cmath>

//------------------------------------------------------------------------------
///  @brief Parity function interface.
//------------------------------------------------------------------------------
class parity {
public:
//------------------------------------------------------------------------------
///  @brief Evaluate parity function
///
///  @param[in] x Argument value.
///  @returns Value of the function.
//------------------------------------------------------------------------------
    virtual double f(const double x) const=0;
//------------------------------------------------------------------------------
///  @brief Evaluate parity derivative.
///
///  @param[in] x Argument value.
///  @returns Value of the function derivative.
//------------------------------------------------------------------------------
    virtual double df(const double x) const=0;
};

//------------------------------------------------------------------------------
///  @brief Sine Parity function interface.
//------------------------------------------------------------------------------
class sine : public parity {
public:
//------------------------------------------------------------------------------
///  @brief Evaluate sine function
///
///  @param[in] x Argument value.
///  @returns Value of the function.
//------------------------------------------------------------------------------
    inline double f(const double x) const final {
        return sin(x);
    };

//------------------------------------------------------------------------------
///  @brief Evaluate sine derivative.
///
///  @param[in] x Argument value.
///  @returns Value of the function derivative.
//------------------------------------------------------------------------------
    inline double df(const double x) const final {
        return cos(x);
    };
};

//------------------------------------------------------------------------------
///  @brief Cosine Parity function interface.
//------------------------------------------------------------------------------
class cosine : public parity {
public:
//------------------------------------------------------------------------------
///  @brief Evaluate cosine function
///
///  @param[in] x Argument value.
///  @returns Value of the function.
//------------------------------------------------------------------------------
    inline double f(const double x) const final {
        return cos(x);
    };

//------------------------------------------------------------------------------
///  @brief Evaluate cosine derivative.
///
///  @param[in] x Argument value.
///  @returns Value of the function derivative.
//------------------------------------------------------------------------------
    inline double df(const double x) const final {
        return -sin(x);
    };
};

#endif
