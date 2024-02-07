//******************************************************************************
///  @file grid_quantity.hpp
///  @brief Contains classes to interpolate full and half grid quanities.
//******************************************************************************

#ifndef grid_quantity_hpp
#define grid_quantity_hpp

#include <vector>
#include <cmath>
#include <cstddef>

#include "parity.hpp"

//------------------------------------------------------------------------------
///  @brief A radial quantity.
//------------------------------------------------------------------------------
class radial_quantity {
public:
///  Buffer containing radial quantity.
    const std::vector<double> buffer;
///  Grid spacing in s.
    const double ds;

//------------------------------------------------------------------------------
///  @brief Construct a radial quantity.
///
///  @param[in] buffer Buffer containing the radial quantity.
//------------------------------------------------------------------------------
    radial_quantity(const std::vector<double> &buffer);

//------------------------------------------------------------------------------
///  @brief vitual interface to get a value at a radial s position.
///
///  @param[in] s Radial s position.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
    virtual double get(const double s) const=0;
};

class half_grid;

//------------------------------------------------------------------------------
///  @brief A full grid quantity.
///
///  Full grid quantities are defined from the magnetic axis to the last closed
///  flux surface.
//------------------------------------------------------------------------------
class full_grid : public radial_quantity {
public:
//------------------------------------------------------------------------------
///  @brief Construct a full grid quantity.
///
///  @param[in] buffer Buffer containing the radial quantity.
//------------------------------------------------------------------------------
    full_grid(const std::vector<double> &buffer);

//------------------------------------------------------------------------------
///  @brief Get a value at a radial s position.
///
///  @param[in] s Radial s position.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
    double get(const double s) const final;

//------------------------------------------------------------------------------
///  @brief Get a radial derivative at a s position.
///
///  @returns The radial derivative on the half grid.
//------------------------------------------------------------------------------
    half_grid get_prime() const;
};

//------------------------------------------------------------------------------
///  @brief A half grid quantity.
///
///  Half grid quantities lay between the full grid points. For half grid
///  quantities the first index is invalid.
//------------------------------------------------------------------------------
class half_grid : public radial_quantity {
public:
//------------------------------------------------------------------------------
///  @brief Construct a half grid quantity.
///
///  @param[in] buffer Buffer containing the radial quantity.
//------------------------------------------------------------------------------
    half_grid(const std::vector<double> &buffer);

//------------------------------------------------------------------------------
///  @brief Get a value at a radial s position.
///
///  @param[in] s Radial s position.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
    double get(const double s) const final;

//------------------------------------------------------------------------------
///  @brief Get a radial derivative at a s position.
///
///  @param[in] m Poloidal mode number.
///  @returns The radial derivative on the half grid.
//------------------------------------------------------------------------------
    full_grid get_prime(const double m) const;
};

//------------------------------------------------------------------------------
///  @brief A radial vmec quantity.
///
///  Half grid quantities lay between the full grid points. For half grid
///  quantities the first index is invalid.
//------------------------------------------------------------------------------
template <class GRID_CLASS> class vmec_grid {
public:
///  Radial grid buffer.
    const GRID_CLASS grid;

//------------------------------------------------------------------------------
///  @brief Vmec radial quantity.
///
///  @param[in] buffer Buffer containing the radial quantity.
//------------------------------------------------------------------------------
    vmec_grid(const std::vector<double> &buffer) : grid(buffer) {}

//------------------------------------------------------------------------------
///  @brief Vmec radial quantity.
///
///  @param[in] grid A radial grid quantity.
//------------------------------------------------------------------------------
    vmec_grid(const GRID_CLASS &grid) : grid(grid.buffer) {}

//------------------------------------------------------------------------------
///  @brief Get a value at a radial s position.
///
///  @param[in] s Radial s position.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
    double get(const double s) const {
        return grid.get(s);
    }
};

//------------------------------------------------------------------------------
///  @brief A radial siesta quantity.
///
///  Half grid quantities lay between the full grid points. For half grid
///  quantities the first index is invalid.
//------------------------------------------------------------------------------
template <class GRID_CLASS> class siesta_grid {
public:
///  Radial grid buffer.
    const GRID_CLASS grid;

//------------------------------------------------------------------------------
///  @brief Siesta radial quantity.
///
///  @param[in] buffer Buffer containing the radial quantity.
//------------------------------------------------------------------------------
    siesta_grid(const std::vector<double> &buffer) : grid(buffer) {}

//------------------------------------------------------------------------------
///  @brief Siesta radial quantity.
///
///  @param[in] grid A radial grid quantity.
//------------------------------------------------------------------------------
    siesta_grid(const GRID_CLASS &grid) : grid(grid.buffer) {}

//------------------------------------------------------------------------------
///  @brief Get a value at a radial s position.
///
///  s_vmec = s_siesta^2
///
///  @param[in] s Radial s position in the vmec grid.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
    double get(const double s) const {
        return grid.get(sqrt(s));
    }
};

//------------------------------------------------------------------------------
///  @brief Type for vmec fourier quantities.
//------------------------------------------------------------------------------
template<class GRID_CLASS> using vmec_quantity = std::vector<vmec_grid<GRID_CLASS> >;

//------------------------------------------------------------------------------
///  @brief A cosine parity vmec quantity.
//------------------------------------------------------------------------------
template <class GRID_CLASS, class PARITY> class vmec_fourier {
public:
///  M modes.
    const std::vector<double> m;
///  N modes.
    const std::vector<double> n;
///  Mode amplitudes.
    const vmec_quantity<GRID_CLASS> quantity;
///  Parity function.
    const PARITY func;

//------------------------------------------------------------------------------
///  @brief Siesta radial quantity.
///
///  @param[in] buffer Buffer containing the radial quantity.
//------------------------------------------------------------------------------
    vmec_fourier(const vmec_quantity<GRID_CLASS> &buffer,
                 const std::vector<double> &m,
                 const std::vector<double> &n) :
    quantity(buffer), m(m), n(n) {}

//------------------------------------------------------------------------------
///  @brief Get a value at a radial s position.
///
///  s_vmec = s_siesta^2
///
///  @param[in] s Radial s position.
///  @param[in] u Poloidal u position.
///  @param[in] v Toroidal v position.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
    double get(const double s,
               const double u,
               const double v) const {
        double temp = 0.0;

        for (size_t i = 0, e = m.size(); i < e; i++) {
            temp += quantity[i].get(s)*func.f(m[i]*u - n[i]*v);
        }

        return temp;
    }

//------------------------------------------------------------------------------
///  @brief Get a poloidal derivative at a radial s position.
///
///  s_vmec = s_siesta^2
///
///  @param[in] s Radial s position.
///  @param[in] u Poloidal u position.
///  @param[in] v Toroidal v position.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
    double get_du(const double s,
                  const double u,
                  const double v) const {
        double temp = 0.0;

        for (size_t i = 0, e = m.size(); i < e; i++) {
            temp += m[i]*quantity[i].get(s)*func.df(m[i]*u - n[i]*v);
        }

        return temp;
    }

//------------------------------------------------------------------------------
///  @brief Get a toroidal derivative at a radial s position.
///
///  s_vmec = s_siesta^2
///
///  @param[in] s Radial s position.
///  @param[in] u Poloidal u position.
///  @param[in] v Toroidal v position.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
    double get_dv(const double s,
                  const double u,
                  const double v) const {
        double temp = 0.0;

        for (size_t i = 0, e = m.size(); i < e; i++) {
            temp -= n[i]*quantity[i].get(s)*func.df(m[i]*u - n[i]*v);
        }

        return temp;
    }
};

//------------------------------------------------------------------------------
///  @brief Type for siesta fourier quantities for a single dimension.
//------------------------------------------------------------------------------
template<class GRID_CLASS> using siesta_quantity_1d = std::vector<siesta_grid<GRID_CLASS> >;
//------------------------------------------------------------------------------
///  @brief Type for siesta fourier quantities.
//------------------------------------------------------------------------------
template<class GRID_CLASS> using siesta_quantity = std::vector<siesta_quantity_1d<GRID_CLASS> >;

//------------------------------------------------------------------------------
///  @brief A cosine parity vmec quantity.
//------------------------------------------------------------------------------
template <class GRID_CLASS, class PARITY> class siesta_fourier {
public:
///  M modes.
    const size_t mpol;
///  N modes.
    const size_t ntor;
///  Toroidal modes.
    const std::vector<int> tor_modes;
///  Number of field periods.
    const size_t nfp;
///  Mode amplitudes.
    const siesta_quantity<GRID_CLASS> quantity;
///  Parity function.
    const PARITY func;

//------------------------------------------------------------------------------
///  @brief Siesta radial quantity.
///
///  @param[in] buffer    Buffer containing the radial quantity.
///  @param[in] mpol      Number of poloidal modes.
///  @param[in] ntor      Number of toroidal modes.
///  @param[in] tor_modes Toroidal modes.
///  @param[in] nfp       Number of field periods.
//------------------------------------------------------------------------------
    siesta_fourier(const siesta_quantity<GRID_CLASS> &buffer,
                   const size_t mpol,
                   const size_t ntor,
                   const std::vector<int> &tor_modes,
                   const size_t nfp) :
    quantity(buffer), mpol(mpol), ntor(ntor), tor_modes(tor_modes), nfp(nfp) {}

//------------------------------------------------------------------------------
///  @brief Get a value at a radial s position.
///
///  s_vmec = s_siesta^2
///
///  @param[in] s Radial s position.
///  @param[in] u Poloidal u position.
///  @param[in] v Toroidal v position.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
    double get(const double s,
               const double u,
               const double v) const {
        double temp = 0.0;

        for (size_t m = 0, em = mpol; m <= em; m++) {
            for (long n = -ntor, en = ntor; n <= en; n++) {
                const size_t ni = n + ntor;
                temp += quantity[m][ni].get(s)*func.f(m*u + tor_modes[ni]*(nfp*v));
            }
        }

        return temp;
    }

//------------------------------------------------------------------------------
///  @brief Get a poloidal derivative at a radial s position.
///
///  s_vmec = s_siesta^2
///
///  @param[in] s Radial s position.
///  @param[in] u Poloidal u position.
///  @param[in] v Toroidal v position.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
    double get_du(const double s,
                  const double u,
                  const double v) const {
        double temp = 0.0;

        for (size_t m = 0, em = mpol; m <= em; m++) {
            for (long n = -ntor, en = ntor; n <= en; n++) {
                const size_t ni = n + ntor;
                temp += m*quantity[m][ni].get(s)*func.df(m*u + tor_modes[ni]*(nfp*v));
            }
        }

        return temp;
    }

//------------------------------------------------------------------------------
///  @brief Get a toroidal derivative at a radial s position.
///
///  s_vmec = s_siesta^2
///
///  @param[in] s Radial s position.
///  @param[in] u Poloidal u position.
///  @param[in] v Toroidal v position.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
    double get_dv(const double s,
                  const double u,
                  const double v) const {
        double temp = 0.0;

        for (size_t m = 0, em = mpol; m <= em; m++) {
            for (long n = -ntor, en = ntor; n <= en; n++) {
                const size_t ni = n + ntor;
                temp += tor_modes[ni]*(nfp*quantity[m][ni].get(s))*func.df(m*u + tor_modes[ni]*(nfp*v));
            }
        }

        return temp;
    }
};

#endif
