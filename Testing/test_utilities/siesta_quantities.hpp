//******************************************************************************
///  @file siesta_quantities.hpp
///  @brief Contains a context to store the siesta quantities.
//******************************************************************************

#ifndef siesta_quantities_hpp
#define siesta_quantities_hpp

#include <string>
#include <netcdf.h>

#include "grid_quantity.hpp"

//------------------------------------------------------------------------------
///  @brief Siesta quantities.
///
///  Class representing a quantities from siesta.
//------------------------------------------------------------------------------
class siesta_quantities {
public:
/// test function
    const siesta_grid<full_grid> test_full;
///  r
    const siesta_fourier<full_grid, cosine> r;
///  drds
    const siesta_fourier<half_grid, cosine> drds;
///  z
    const siesta_fourier<full_grid, sine> z;
///  dzds
    const siesta_fourier<half_grid, sine> dzds;
///  Radial derivative of toroidal flux.
    const siesta_grid<full_grid> phipf;
///  Radial derivative of poloidal flux.
    const siesta_grid<full_grid> chipf;
///  JB^s
    const siesta_fourier<half_grid, sine> jbsups;
///  dJB^sds
    const siesta_fourier<full_grid, sine> djbsupsds;
///  JB^u
    const siesta_fourier<half_grid, cosine> jbsupu;
///  JB^v
    const siesta_fourier<half_grid, cosine> jbsupv;
///  JB^s
    const siesta_fourier<half_grid, sine> bsups;
///  JB^u
    const siesta_fourier<half_grid, cosine> bsupu;
///  JB^v
    const siesta_fourier<half_grid, cosine> bsupv;
///  B_s
    const siesta_fourier<half_grid, sine> bsubs;
///  B_u
    const siesta_fourier<half_grid, cosine> bsubu;
///  B_v
    const siesta_fourier<half_grid, cosine> bsubv;
///  JK^s
    const siesta_fourier<full_grid, sine> jksups;
///  JK^u
    const siesta_fourier<full_grid, cosine> jksupu;
///  JK^v
    const siesta_fourier<full_grid, cosine> jksupv;
///  Total toroidal current.
    const double curtor;
///  Pressure
    const siesta_fourier<half_grid, cosine> p;

//------------------------------------------------------------------------------
///  @brief Siesta quantities.
///
///  For now this only supports symmetric cases only.
///
///  @param[in] restart_file File name of a restart file.
//------------------------------------------------------------------------------
    siesta_quantities(const std::string &restart_file) :
    test_full(make_full()),
    r(load_fourier<full_grid, cosine> (restart_file, "rmnc_m_n_r_")),
    drds(to_prime(r)),
    z(load_fourier<full_grid, sine> (restart_file, "zmns_m_n_r_")),
    dzds(to_prime(z)),
    phipf(load<full_grid> (restart_file, "phipf_r_")),
    chipf(load<full_grid> (restart_file, "chipf_r_")),
    jbsups(load_fourier_denorm<half_grid, sine> (restart_file, "JBsupssh_m_n_r_")),
    djbsupsds(to_prime(jbsups)),
    jbsupu(load_fourier_denorm<half_grid, cosine> (restart_file, "JBsupuch_m_n_r_")),
    jbsupv(load_fourier_denorm<half_grid, cosine> (restart_file, "JBsupvch_m_n_r_")),
    bsups(load_fourier<half_grid, sine> (restart_file, "bsupsmnsh_m_n_r_")),
    bsupu(load_fourier<half_grid, cosine> (restart_file, "bsupumnch_m_n_r_")),
    bsupv(load_fourier<half_grid, cosine> (restart_file, "bsupvmnch_m_n_r_")),
    bsubs(load_fourier<half_grid, sine> (restart_file, "bsubsmnsh_m_n_r_")),
    bsubu(load_fourier<half_grid, cosine> (restart_file, "bsubumnch_m_n_r_")),
    bsubv(load_fourier<half_grid, cosine> (restart_file, "bsubvmnch_m_n_r_")),
    jksups(load_fourier<full_grid, sine> (restart_file, "jksupsmnsf_m_n_r_")),
    jksupu(load_fourier<full_grid, cosine> (restart_file, "jksupumncf_m_n_r_")),
    jksupv(load_fourier<full_grid, cosine> (restart_file, "jksupvmncf_m_n_r_")),
    curtor(load_scalar(restart_file, "curtor")),
    p(load_fourier<half_grid, cosine> (restart_file, "pmnch_m_n_r_")) {}

//------------------------------------------------------------------------------
///  @brief Factory method to make a test quantity.
///
///  @returns A constructed test quantity.
//------------------------------------------------------------------------------
    static siesta_grid<full_grid> make_full() {
        std::vector<double> temp(100);

        for (size_t i = 0; i < 100; i++) {
            double s = static_cast<double> (i)/99.0;
            s *= s;

            temp[i] = cos(s);
        }

        return siesta_grid<full_grid> (temp);
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a siesta quantity.
///
///  @param[in] restart_file File name of a siesta restart file.
///  @param[in] name          Name of the quantity to load.
///  @returns A constructed siesta grid quantity.
//------------------------------------------------------------------------------
    template<class GIRD_CLASS>
    static siesta_grid<GIRD_CLASS> load(const std::string &restart_file,
                                        const std::string &name) {
        int ncid;
        int status = nc_open(restart_file.c_str(), NC_NOWRITE, &ncid);
        if (status) {
            std::cout << "Failed to open " << restart_file << std::endl;
            exit(status);
        }

        int varid;
        int nrad;
        nc_inq_varid(ncid, "nrad", &varid);
        nc_get_var(ncid, varid, &nrad);

        std::vector<double> temp(nrad);
        nc_inq_varid(ncid, name.c_str(), &varid);
        nc_get_var(ncid, varid, temp.data());

        nc_close(ncid);

        return siesta_grid<GIRD_CLASS> (temp);
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a scalar siesta quantity.
///
///  @param[in] restart_file File name of a siesta restart file.
///  @param[in] name          Name of the quantity to load.
///  @returns A constructed scalar quantity.
//------------------------------------------------------------------------------
    static double load_scalar(const std::string &restart_file,
                              const std::string &name) {
        int ncid;
        int status = nc_open(restart_file.c_str(), NC_NOWRITE, &ncid);
        if (status) {
            std::cout << "Failed to open " << restart_file << std::endl;
            exit(status);
        }

        int varid;
        double temp;
        nc_inq_varid(ncid, "curtor", &varid);
        nc_get_var(ncid, varid, &temp);

        nc_close(ncid);

        return temp;
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a fourier siesta quantity.
///
///  @param[in] restart_file File name of a siesta restart file.
///  @param[in] name         Name of the quantity to load.
///  @returns A constructed siesta fourier quantity.
//------------------------------------------------------------------------------
    template<class GIRD_CLASS, class PARITY>
    static siesta_fourier<GIRD_CLASS, PARITY> load_fourier(const std::string &restart_file,
                                                           const std::string &name) {
        int ncid;
        int status = nc_open(restart_file.c_str(), NC_NOWRITE, &ncid);
        if (status) {
            std::cout << "Failed to open " << restart_file << std::endl;
            exit(status);
        }

        int varid;
        int nrad;
        nc_inq_varid(ncid, "nrad", &varid);
        nc_get_var(ncid, varid, &nrad);

        int mpol;
        nc_inq_varid(ncid, "mpol", &varid);
        nc_get_var(ncid, varid, &mpol);

        int ntor;
        nc_inq_varid(ncid, "ntor", &varid);
        nc_get_var(ncid, varid, &ntor);

        int nfp;
        nc_inq_varid(ncid, "nfp", &varid);
        nc_get_var(ncid, varid, &nfp);

        siesta_quantity<GIRD_CLASS> quantity;
        for (size_t m = 0, em = mpol; m <= em; m++) {
            for (long n = -ntor, en = ntor; n <= en; n++) {
                const size_t ni = n + ntor;

                std::vector<double> temp(nrad);
                nc_inq_varid(ncid, name.c_str(), &varid);

                const std::array<size_t, 3> start = {0, ni, m};
                const std::array<size_t, 3> end = {static_cast<size_t> (nrad), 1, 1};
                nc_get_vara(ncid, varid,
                            start.data(), end.data(),
                            temp.data());

                if (n == -ntor) {
                    siesta_quantity_1d<GIRD_CLASS> new_grid;
                    new_grid.push_back(temp);
                    quantity.push_back(new_grid);
                } else {
                    quantity[m].push_back(siesta_grid<GIRD_CLASS> (temp));
                }
            }
        }

        nc_close(ncid);

        return siesta_fourier<GIRD_CLASS, PARITY> (quantity, mpol, ntor, nfp);
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a fourier denormalized siesta quantity.
///
///  @param[in] restart_file File name of a siesta restart file.
///  @param[in] name         Name of the quantity to load.
///  @returns A constructed denormalized siesta fourier quantity.
//------------------------------------------------------------------------------
    template<class GIRD_CLASS, class PARITY>
    static siesta_fourier<GIRD_CLASS, PARITY> load_fourier_denorm(const std::string &restart_file,
                                                                  const std::string &name) {
        int ncid;
        int status = nc_open(restart_file.c_str(), NC_NOWRITE, &ncid);
        if (status) {
            std::cout << "Failed to open " << restart_file << std::endl;
            exit(status);
        }

        int varid;
        int nrad;
        nc_inq_varid(ncid, "nrad", &varid);
        nc_get_var(ncid, varid, &nrad);

        int mpol;
        nc_inq_varid(ncid, "mpol", &varid);
        nc_get_var(ncid, varid, &mpol);

        int ntor;
        nc_inq_varid(ncid, "ntor", &varid);
        nc_get_var(ncid, varid, &ntor);

        int nfp;
        nc_inq_varid(ncid, "nfp", &varid);
        nc_get_var(ncid, varid, &nfp);

        double bfactor;
        nc_inq_varid(ncid, "b_factor", &varid);
        nc_get_var(ncid, varid, &bfactor);

        siesta_quantity<GIRD_CLASS> quantity;
        for (size_t m = 0, em = mpol; m <= em; m++) {
            for (long n = -ntor, en = ntor; n <= en; n++) {
                const size_t ni = n + ntor;

                std::vector<double> temp(nrad);
                nc_inq_varid(ncid, name.c_str(), &varid);

                const std::array<size_t, 3> start = {0, ni, m};
                const std::array<size_t, 3> end = {static_cast<size_t> (nrad), 1, 1};
                nc_get_vara(ncid, varid,
                            start.data(), end.data(),
                            temp.data());

                for (size_t i = 0; i < nrad; i++) {
                    if (m == 0 && n == 0) {
                        temp[i] = temp[i]/bfactor;
                    } else {
                        temp[i] = sqrt(2.0)*temp[i]/bfactor;
                    }
                }

                if (n == -ntor) {
                    siesta_quantity_1d<GIRD_CLASS> new_grid;
                    new_grid.push_back(temp);
                    quantity.push_back(new_grid);
                } else {
                    quantity[m].push_back(siesta_grid<GIRD_CLASS> (temp));
                }
            }
        }

        nc_close(ncid);

        return siesta_fourier<GIRD_CLASS, PARITY> (quantity, mpol, ntor, nfp);
    }

//------------------------------------------------------------------------------
///  @brief Convert full grid quantity to half grid primed.
///
///  @param[in] siesta A full grid siesta quantity to take the derivative of.
///  @returns A half grid derivative.
//------------------------------------------------------------------------------
    template<class PARITY>
    static siesta_fourier<half_grid, PARITY> to_prime(const siesta_fourier<full_grid, PARITY> siesta) {
        siesta_quantity<half_grid> half_quanity;

        for (const siesta_quantity_1d<full_grid> slice : siesta.quantity) {
            siesta_quantity_1d<half_grid> half_slice;

            for (const siesta_grid<full_grid> &q : slice) {
                half_slice.push_back(siesta_grid<half_grid> (q.grid.get_prime()));
            }

            half_quanity.push_back(half_slice);
        }

        return siesta_fourier<half_grid, PARITY> (half_quanity,
                                                  siesta.mpol,
                                                  siesta.ntor,
                                                  siesta.nfp);
    }

//------------------------------------------------------------------------------
///  @brief Convert full grid quantity to half grid primed.
///
///  @param[in] siesta A full grid siesta quantity to take the derivative of.
///  @returns A full grid derivative.
//------------------------------------------------------------------------------
    template<class PARITY>
    static siesta_fourier<full_grid, PARITY> to_prime(const siesta_fourier<half_grid, PARITY> siesta) {
        siesta_quantity<full_grid> full_quanity;

        for (size_t m = 0, em = siesta.mpol; m <= em; m++) {
            siesta_quantity_1d<full_grid> full_slice;

            for (long n = - siesta.ntor, en = siesta.ntor; n <= en; n++) {
                const size_t ni = n + siesta.ntor;
                full_slice.push_back(siesta_grid<full_grid> (siesta.quantity[m][ni].grid.get_prime(static_cast<double> (m))));
            }

            full_quanity.push_back(full_slice);
        }

        return siesta_fourier<full_grid, PARITY> (full_quanity,
                                                  siesta.mpol,
                                                  siesta.ntor,
                                                  siesta.nfp);
    }

//------------------------------------------------------------------------------
///  @brief Get R quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns R at the s,u,v position.
//------------------------------------------------------------------------------
    double get_r(const double s, const double u, const double v) const {
        return r.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get R' quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns R' at the s,u,v position.
//------------------------------------------------------------------------------
    double get_r_prime(const double s, const double u, const double v) const {
        return drds.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get dRdu quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns dRdu at the s,u,v position.
//------------------------------------------------------------------------------
    double get_dr_du(const double s, const double u, const double v) const {
        return r.get_du(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get dRdv quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns dRdv at the s,u,v position.
//------------------------------------------------------------------------------
    double get_dr_dv(const double s, const double u, const double v) const {
        return r.get_dv(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get Z quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns Z at the s,u,v position.
//------------------------------------------------------------------------------
    double get_z(const double s, const double u, const double v) const {
        return z.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get Z' quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns R' at the s,u,v position.
//------------------------------------------------------------------------------
    double get_z_prime(const double s, const double u, const double v) const {
        return dzds.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get dZdu quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns dZdu at the s,u,v position.
//------------------------------------------------------------------------------
    double get_dz_du(const double s, const double u, const double v) const {
        return z.get_du(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get dZdv quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns dZdv at the s,u,v position.
//------------------------------------------------------------------------------
    double get_dz_dv(const double s, const double u, const double v) const {
        return z.get_dv(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get phipf quantity.
///
///  @param[in] s Radial position.
///  @returns phipf at the s position.
//------------------------------------------------------------------------------
    double get_phipf(const double s) const {
        return phipf.get(s);
    }

//------------------------------------------------------------------------------
///  @brief Get chipf quantity.
///
///  @param[in] s Radial position.
///  @returns chipf at the s position.
//------------------------------------------------------------------------------
    double get_chipf(const double s) const {
        return chipf.get(s);
    }

//------------------------------------------------------------------------------
///  @brief Get bsups quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns bsups at the s,u,v position.
//------------------------------------------------------------------------------
    double get_bsups(const double s, const double u, const double v) const {
        return bsups.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get bsupu quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns bsupu at the s,u,v position.
//------------------------------------------------------------------------------
    double get_bsupu(const double s, const double u, const double v) const {
        return bsupu.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get bsupv quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns bsupv at the s,u,v position.
//------------------------------------------------------------------------------
    double get_bsupv(const double s, const double u, const double v) const {
        return bsupv.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get bsubs quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns bsubs at the s,u,v position.
//------------------------------------------------------------------------------
    double get_bsubs(const double s, const double u, const double v) const {
        return bsubs.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get bsubu quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns bsubu at the s,u,v position.
//------------------------------------------------------------------------------
    double get_bsubu(const double s, const double u, const double v) const {
        return bsubu.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get bsubv quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns bsubv at the s,u,v position.
//------------------------------------------------------------------------------
    double get_bsubv(const double s, const double u, const double v) const {
        return bsubv.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get jksups quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns jksups at the s,u,v position.
//------------------------------------------------------------------------------
    double get_jksups(const double s, const double u, const double v) const {
        return jksups.get(s, u, v)/(M_PI*4.0E-7);
    }

//------------------------------------------------------------------------------
///  @brief Get jksupu quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns jksupu at the s,u,v position.
//------------------------------------------------------------------------------
    double get_jksupu(const double s, const double u, const double v) const {
        return jksupu.get(s, u, v)/(M_PI*4.0E-7);
    }

//------------------------------------------------------------------------------
///  @brief Get jksupv quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns jksupv at the s,u,v position.
//------------------------------------------------------------------------------
    double get_jksupv(const double s, const double u, const double v) const {
        return jksupv.get(s, u, v)/(M_PI*4.0E-7);
    }

//------------------------------------------------------------------------------
///  @brief Get jbsups quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns jbsups at the s,u,v position.
//------------------------------------------------------------------------------
    double get_jbsups(const double s, const double u, const double v) const {
        return jbsups.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get jbsupu quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns jbsupu at the s,u,v position.
//------------------------------------------------------------------------------
    double get_jbsupu(const double s, const double u, const double v) const {
        return jbsupu.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get jbsupv quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns jbsupv at the s,u,v position.
//------------------------------------------------------------------------------
    double get_jbsupv(const double s, const double u, const double v) const {
        return jbsupv.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get jacobian.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns jacobian at the s,u,v position.
//------------------------------------------------------------------------------
    double get_jacobian(const double s, const double u, const double v) const {
        return get_r(s, u, v)*(get_dr_du(s, u, v)*get_z_prime(s, u, v) -
                               get_r_prime(s, u, v)*get_dz_du(s, u, v));
    }

//------------------------------------------------------------------------------
///  @brief Get total toroidal current.
///
///  @returns The total toroidal current.
//------------------------------------------------------------------------------
    double get_curtor() const {
        return curtor;
    }

//------------------------------------------------------------------------------
///  @brief Get djbsupsds quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns djbsupsds at the s,u,v position.
//------------------------------------------------------------------------------
    double get_djbsups_ds(const double s, const double u, const double v) const {
        return djbsupsds.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get djbsupudu quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns djbsupudu at the s,u,v position.
//------------------------------------------------------------------------------
    double get_djbsupu_du(const double s, const double u, const double v) const {
        return jbsupu.get_du(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get djbsupvdv quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns djbsupvdv at the s,u,v position.
//------------------------------------------------------------------------------
    double get_djbsupv_dv(const double s, const double u, const double v) const {
        return jbsupv.get_dv(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get divergence of B.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns Div B at the s,u,v position.
//------------------------------------------------------------------------------
    double get_divb(const double s, const double u, const double v) const {
        return get_djbsups_ds(s,u,v) + get_djbsupu_du(s,u,v) + get_djbsupv_dv(s,u,v);
    }

//------------------------------------------------------------------------------
///  @brief Get pressure.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns Pressure at the s,u,v position.
//------------------------------------------------------------------------------
    double get_pressure(const double s, const double u, const double v) const {
        return p.get(s,u,v);
    }

//------------------------------------------------------------------------------
///  @brief Get test quantity.
///
///  @param[in] s Radial position.
///  @returns Test value at the s position.
//------------------------------------------------------------------------------
    double get_test(const double s) const {
        return test_full.get(s);
    }
};

#endif
