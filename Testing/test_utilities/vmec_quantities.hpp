//******************************************************************************
///  @file vmec_quantities.hpp
///  @brief Contains a context to store the vmec quantities.
//******************************************************************************

#ifndef vmec_quantities_hpp
#define vmec_quantities_hpp

#include <string>
#include <netcdf.h>
#include <array>

#include "grid_quantity.hpp"

//------------------------------------------------------------------------------
///  @brief Vmec quantities.
///
///  Class representing a quantities from vmec.
//------------------------------------------------------------------------------
class vmec_quantities {
public:
/// test function
    const vmec_grid<full_grid> test_full;
///  r
    const vmec_fourier<full_grid, cosine> r;
///  drds
    const vmec_fourier<half_grid, cosine> drds;
///  z
    const vmec_fourier<full_grid, sine> z;
///  dzds
    const vmec_fourier<half_grid, sine> dzds;
///  Radial derivative of toroidal flux.
    const vmec_grid<full_grid> phipf;
///  Radial derivative of poloidal flux.
    const vmec_grid<full_grid> chipf;
///  Sign of vmec jacobian.
    const double signj;
///  Jacobian
    const vmec_fourier<half_grid, cosine> j;
///  B^u
    const vmec_fourier<half_grid, cosine> bsupu;
///  B^v
    const vmec_fourier<half_grid, cosine> bsupv;
///  B_s
    const vmec_fourier<half_grid, sine> bsubs;
///  B_u
    const vmec_fourier<half_grid, cosine> bsubu;
///  B_v
    const vmec_fourier<half_grid, cosine> bsubv;
///  JK^u
    const vmec_fourier<full_grid, cosine> jksupu;
///  JK^v
    const vmec_fourier<full_grid, cosine> jksupv;
///  Total toroidal current.
    const double curtor;
///  Pressure
    const vmec_grid<full_grid> p;

//------------------------------------------------------------------------------
///  @brief Vmec quantities.
///
///  For now this only supports symmetric cases only.
///
///  @param[in] wout_file File name of a vmec wout file.
//------------------------------------------------------------------------------
    vmec_quantities(const std::string &wout_file) :
    test_full(make_full()),
    r(load_fourier<full_grid, cosine> (wout_file, "rmnc")),
    drds(to_prime(r)),
    z(load_fourier<full_grid, sine> (wout_file, "zmns")),
    dzds(to_prime(z)),
    phipf(load<full_grid> (wout_file, "phipf")),
    chipf(load<full_grid> (wout_file, "chipf")),
    signj(load_scalar<int> (wout_file, "signgs")),
    j(load_fourier_nyq<half_grid, cosine> (wout_file, "gmnc")),
    bsupu(load_fourier_nyq<half_grid, cosine> (wout_file, "bsupumnc")),
    bsupv(load_fourier_nyq<half_grid, cosine> (wout_file, "bsupvmnc")),
    bsubs(load_fourier_nyq<half_grid, sine> (wout_file, "bsubsmns")),
    bsubu(load_fourier_nyq<half_grid, cosine> (wout_file, "bsubumnc")),
    bsubv(load_fourier_nyq<half_grid, cosine> (wout_file, "bsubvmnc")),
    jksupu(load_fourier_nyq<full_grid, cosine> (wout_file, "currumnc")),
    jksupv(load_fourier_nyq<full_grid, cosine> (wout_file, "currvmnc")),
    curtor(load_scalar<double> (wout_file, "ctor")),
    p(load<full_grid> (wout_file, "presf")) {}

//------------------------------------------------------------------------------
///  @brief Factory method to make a test quantity.
///
///  @returns A constructed test quantity.
//------------------------------------------------------------------------------
    static vmec_grid<full_grid> make_full() {
        std::vector<double> temp(100);

        for (size_t i = 0; i < 100; i++) {
            const double s = static_cast<double> (i)/99.0;

            temp[i] = cos(s);
        }

        return vmec_grid<full_grid> (temp);
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a scalar vmec quantity.
///
///  @param[in] wout_file File name of a vmec wout file.
///  @param[in] name      Name of the quantity to load.
///  @returns A constructed scalar quantity.
//------------------------------------------------------------------------------
    template<typename TYPE>
    static double load_scalar(const std::string &wout_file,
                              const std::string &name) {
        int ncid;
        nc_open(wout_file.c_str(), NC_NOWRITE, &ncid);

        int varid;
        TYPE temp;
        nc_inq_varid(ncid, name.c_str(), &varid);
        nc_get_var(ncid, varid, &temp);

        nc_close(ncid);

        return static_cast<double> (temp);
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a vmec quantity.
///
///  @param[in] wout_file File name of a vmec wout file.
///  @param[in] name      Name of the quantity to load.
///  @returns A constructed vmec grid quantity.
//------------------------------------------------------------------------------
    template<class GIRD_CLASS>
    static vmec_grid<GIRD_CLASS> load(const std::string &wout_file,
                                      const std::string &name) {
        int ncid;
        nc_open(wout_file.c_str(), NC_NOWRITE, &ncid);

        int varid;
        int ns;
        nc_inq_varid(ncid, "ns", &varid);
        nc_get_var(ncid, varid, &ns);

        std::vector<double> temp(ns);
        nc_inq_varid(ncid, name.c_str(), &varid);
        nc_get_var(ncid, varid, temp.data());

        nc_close(ncid);

        return vmec_grid<GIRD_CLASS> (temp);
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a fourier vmec quantity.
///
///  @param[in] wout_file File name of a vmec wout file.
///  @param[in] name      Name of the quantity to load.
///  @returns A constructed vmec fourier quantity.
//------------------------------------------------------------------------------
    template<class GIRD_CLASS, class PARITY>
    static vmec_fourier<GIRD_CLASS, PARITY> load_fourier(const std::string &wout_file,
                                                         const std::string &name) {
        int ncid;
        nc_open(wout_file.c_str(), NC_NOWRITE, &ncid);

        int varid;
        int ns;
        nc_inq_varid(ncid, "ns", &varid);
        nc_get_var(ncid, varid, &ns);

        int mnmax;
        nc_inq_varid(ncid, "mnmax", &varid);
        nc_get_var(ncid, varid, &mnmax);

        std::vector<double> xm(mnmax);
        nc_inq_varid(ncid, "xm", &varid);
        nc_get_var(ncid, varid, xm.data());

        std::vector<double> xn(mnmax);
        nc_inq_varid(ncid, "xn", &varid);
        nc_get_var(ncid, varid, xn.data());

        vmec_quantity<GIRD_CLASS> quantity;
        for (size_t i = 0; i < mnmax; i++) {
            std::vector<double> temp(ns);
            nc_inq_varid(ncid, name.c_str(), &varid);

            const std::array<size_t, 2> start = {0, i};
            const std::array<size_t, 2> end = {static_cast<size_t> (ns), 1};
            nc_get_vara(ncid, varid,
                        start.data(), end.data(),
                        temp.data());

            quantity.push_back(vmec_grid<GIRD_CLASS> (temp));
        }

        nc_close(ncid);

        return vmec_fourier<GIRD_CLASS, PARITY> (quantity, xm, xn);
    }

//------------------------------------------------------------------------------
///  @brief Factory method to load a nyquest fourier vmec quantity.
///
///  @param[in] wout_file File name of a vmec wout file.
///  @param[in] name      Name of the nyquest quantity to load.
///  @returns A constructed vmec nyquest fourier quantity.
//------------------------------------------------------------------------------
    template<class GIRD_CLASS, class PARITY>
    static vmec_fourier<GIRD_CLASS, PARITY> load_fourier_nyq(const std::string &wout_file,
                                                             const std::string &name) {
        int ncid;
        nc_open(wout_file.c_str(), NC_NOWRITE, &ncid);

        int varid;
        int ns;
        nc_inq_varid(ncid, "ns", &varid);
        nc_get_var(ncid, varid, &ns);

        int mnmax;
        nc_inq_varid(ncid, "mnmax_nyq", &varid);
        nc_get_var(ncid, varid, &mnmax);

        std::vector<double> xm(mnmax);
        nc_inq_varid(ncid, "xm_nyq", &varid);
        nc_get_var(ncid, varid, xm.data());

        std::vector<double> xn(mnmax);
        nc_inq_varid(ncid, "xn_nyq", &varid);
        nc_get_var(ncid, varid, xn.data());

        vmec_quantity<GIRD_CLASS> quantity;
        for (size_t i = 0; i < mnmax; i++) {
            std::vector<double> temp(ns);
            nc_inq_varid(ncid, name.c_str(), &varid);

            const std::array<size_t, 2> start = {0, i};
            const std::array<size_t, 2> end = {static_cast<size_t> (ns), 1};
            nc_get_vara(ncid, varid,
                        start.data(), end.data(),
                        temp.data());

            quantity.push_back(vmec_grid<GIRD_CLASS> (temp));
        }

        nc_close(ncid);

        return vmec_fourier<GIRD_CLASS, PARITY> (quantity, xm, xn);
    }

//------------------------------------------------------------------------------
///  @brief Convert full grid quantity to half grid primed.
///
///  @param[in] vmec A full grid vmec quantity to take the derivative of.
///  @returns A half grid derivative.
//------------------------------------------------------------------------------
    template<class PARITY>
    static vmec_fourier<half_grid, PARITY> to_prime(const vmec_fourier<full_grid, PARITY> vmec) {
        vmec_quantity<half_grid> half_quanity;

        for (const vmec_grid<full_grid> &q : vmec.quantity) {
            half_quanity.push_back(vmec_grid<half_grid> (q.grid.get_prime()));
        }

        return vmec_fourier<half_grid, PARITY> (half_quanity,  vmec.m,  vmec.n);
    }

//------------------------------------------------------------------------------
///  @brief Convert full grid quantity to half grid primed.
///
///  @param[in] vmec A full grid vmec quantity to take the derivative of.
///  @returns A full grid derivative.
//------------------------------------------------------------------------------
    template<class PARITY>
    static vmec_fourier<full_grid, PARITY> to_prime(const vmec_fourier<half_grid, PARITY> vmec) {
        vmec_quantity<full_grid> full_quanity;

        for (size_t i = 0, e = vmec.m.size(); i < e; i++) {
            full_quanity.push_back(vmec_grid<full_grid> (vmec.quantity[i].grid.get_prime(vmec.m[i])));
        }

        return vmec_fourier<full_grid, PARITY> (full_quanity,  vmec.m,  vmec.n);
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
        return 2.0*sqrt(s)*drds.get(s, u, v);
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
        return 2.0*sqrt(s)*dzds.get(s, u, v);
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
        return signj*sqrt(s)*phipf.get(s)/M_PI;
    }

//------------------------------------------------------------------------------
///  @brief Get chipf quantity.
///
///  @param[in] s Radial position.
///  @returns chipf at the s position.
//------------------------------------------------------------------------------
    double get_chipf(const double s) const {
        return signj*sqrt(s)*chipf.get(s)/M_PI;
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
        return 0.0;
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
        return 0.0;
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
        return 2.0*sqrt(s)*jksupu.get(s, u, v);
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
        return 2.0*sqrt(s)*jksupv.get(s, u, v);
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
        return 0.0;
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
        return 2.0*sqrt(s)*j.get(s, u, v)*get_bsupu(s, u, v);
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
        return 2.0*sqrt(s)*j.get(s, u, v)*get_bsupv(s, u, v);
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
        return 2.0*sqrt(s)*j.get(s, u, v);
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
///  @brief Get dbsupudu quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns dbsupudu at the s,u,v position.
//------------------------------------------------------------------------------
    double get_dbsupu_du(const double s, const double u, const double v) const {
        return 2.0*sqrt(s)*j.get(s, u, v)*bsupu.get_du(s, u, v) + 2.0*sqrt(s)*j.get_du(s, u, v)*bsupu.get(s, u, v);
    }

//------------------------------------------------------------------------------
///  @brief Get dbsupvdv quantity.
///
///  @param[in] s Radial position.
///  @param[in] u Radial position.
///  @param[in] v Radial position.
///  @returns dbsupvdv at the s,u,v position.
//------------------------------------------------------------------------------
    double get_dbsupv_dv(const double s, const double u, const double v) const {
        return 2.0*sqrt(s)*j.get(s, u, v)*bsupv.get_dv(s, u, v) + 2.0*sqrt(s)*j.get_dv(s, u, v)*bsupv.get(s, u, v);
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
        return get_dbsupu_du(s,u,v) + get_dbsupv_dv(s,u,v);
    }

//------------------------------------------------------------------------------
///  @brief Get pressure quantity.
///
///  @param[in] s Radial position.
///  @returns Pressure at the s position.
//------------------------------------------------------------------------------
    double get_pressure(const double s) const {
        return p.get(s);
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
