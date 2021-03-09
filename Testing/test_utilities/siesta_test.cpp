#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include <functional>

#include "vmec_quantities.hpp"
#include "siesta_quantities.hpp"
#include "commandline_parser.hpp"

#define RUN_TEST_3D(vmec, siesta, method, name, args)                 \
if (args.is_set("-dump")) {                                           \
    dump([&vmec]   (const double s, const double u, const double v) { \
             return vmec.method(s,u,v);                               \
         },                                                           \
         [&siesta] (const double s, const double u, const double v) { \
             return siesta.method(s,u,v);                             \
         },                                                           \
         args.get<double> ("-u"),                                     \
         args.get<double> ("-v"),                                     \
         args.is_set("-relative"));                                   \
} else {                                                              \
    test([&vmec]   (const double s, const double u, const double v) { \
             return vmec.method(s,u,v);                               \
         },                                                           \
         [&siesta] (const double s, const double u, const double v) { \
             return siesta.method(s,u,v);                             \
         },                                                           \
         name,                                                        \
         args.get<double> ("-tol"),                                   \
         args.get<double> ("-min"),                                   \
         args.get<double> ("-max"),                                   \
         args.is_set("-relative"));                                   \
}

#define RUN_TEST_13D(vmec, siesta, method, name, args) \
if (args.is_set("-dump")) {                            \
    dump([&vmec]   (const double s) {                  \
            return vmec.method(s);                     \
         },                                            \
         [&siesta] (const double s) {                  \
            return siesta.method(s, 0.0, 0.0);         \
         },                                            \
         args.is_set("-relative"));                    \
} else {                                               \
    test([&vmec]   (const double s) {                  \
            return vmec.method(s);                     \
         },                                            \
         [&siesta] (const double s) {                  \
            return siesta.method(s, 0.0, 0.0);         \
         },                                            \
         name,                                         \
         args.get<double> ("-tol"),                    \
         args.get<double> ("-min"),                    \
         args.get<double> ("-max"),                    \
         args.is_set("-relative"));                    \
}

#define RUN_TEST_1D(vmec, siesta, method, name, arsg)           \
if (args.is_set("-dump")) {                                     \
    dump([&vmec]   (const double s) {return vmec.method(s);},   \
         [&siesta] (const double s) {return siesta.method(s);}, \
         args.is_set("-relative"));                             \
} else {                                                        \
    test([&vmec]   (const double s) {return vmec.method(s);},   \
         [&siesta] (const double s) {return siesta.method(s);}, \
         name,                                                  \
         args.get<double> ("-tol"),                             \
         args.get<double> ("-min"),                             \
         args.get<double> ("-max"),                             \
         args.is_set("-relative"));                             \
}

#define RUN_TEST_0D(vmec, siesta, method, name, args) \
if (args.is_set("-dump")) {                           \
    dump([&vmec]   () {return vmec.method();},        \
         [&siesta] () {return siesta.method();},      \
         args.is_set("-relative"));                   \
} else {                                              \
    test([&vmec]   () {return vmec.method();},        \
         [&siesta] () {return siesta.method();},      \
         name, args.get<double> ("-tol"),             \
         args.is_set("-relative"));                   \
}

#define RUN_TEST_DIV(eq, method, name, args)                  \
test([&eq] (const double s, const double u, const double v) { \
        return eq.method(s,u,v);                              \
     },                                                       \
     name,                                                    \
     args.get<double> ("-tol"),                               \
     args.get<double> ("-min"),                               \
     args.get<double> ("-max"))

//------------------------------------------------------------------------------
///  @brief Dump vmec and siesta scalar quantities.
///
///  @param[in] vmec     Getter function for the vmec value.
///  @param[in] siesta   Getter function for the siesta value.
///  @param[in] relative Use relative error.
//------------------------------------------------------------------------------
void dump(std::function<double()> vmec,
          std::function<double()> siesta,
          const bool relative) {
    const double vmec_value = vmec();
    const double siesta_value = siesta();
    double diff = abs(vmec_value - siesta_value);

    if (relative) {
        diff /= std::max(abs(vmec_value), abs(siesta_value));
    }

    std::cout << vmec_value   << " ";
    std::cout << siesta_value << " ";
    std::cout << diff << std::endl;
}

//------------------------------------------------------------------------------
///  @brief Dump vmec and siesta 1D profile quantities.
///
///  @param[in] vmec     Getter function for the vmec value.
///  @param[in] siesta   Getter function for the siesta value.
///  @param[in] relative Use relative error.
//------------------------------------------------------------------------------
void dump(std::function<double(const double s)> vmec,
          std::function<double(const double s)> siesta,
          const bool relative) {
    for (double s = 0; s <= 1.0; s += 0.005) {
        const double vmec_value = vmec(s);
        const double siesta_value = siesta(s);
        double diff = abs(vmec_value - siesta_value);

        if (relative) {
            diff /= std::max(abs(vmec_value), abs(siesta_value));
        }

        std::cout << s << " ";
        std::cout << vmec_value   << " ";
        std::cout << siesta_value << " ";
        std::cout << diff << std::endl;
    }
}

//------------------------------------------------------------------------------
///  @brief Dump vmec and siesta 3D profile quantities.
///
///  @param[in] vmec     Getter function for the vmec value.
///  @param[in] siesta   Getter function for the siesta value.
///  @param[in] u        Poloidal angle.
///  @param[in] v        Toroidal angle.
///  @param[in] relative Use relative error.
//------------------------------------------------------------------------------
void dump(std::function<double(const double, const double, const double)> vmec,
          std::function<double(const double, const double, const double)> siesta,
          const double u, const double v,
          const bool relative) {
    for (double s = 0; s <= 1.0; s += 0.005) {
        const double vmec_value = vmec(s, u, v);
        const double siesta_value = siesta(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        if (relative) {
            diff /= std::max(abs(vmec_value), abs(siesta_value));
        }

        std::cout << s << " ";
        std::cout << vmec_value   << " ";
        std::cout << siesta_value << " ";
        std::cout << diff << std::endl;
    }
}

//------------------------------------------------------------------------------
///  @brief Test comparitor.
///
///  Randomly samples s,u,v positions to check that the value fall with in the
///  tolerance value.
///
///  @param[in] vmec      Getter function for the vmec value.
///  @param[in] siesta    Getter function for the siesta value.
///  @param[in] name      Name of the test.
///  @param[in] tolarance Tolarnce threshold that values must fall with in.
///  @param[in] smin      Minimum radial extent.
///  @param[in] smax      Maximum radial extent.
///  @param[in] relative  Use relative error.
//------------------------------------------------------------------------------
void test(std::function<double(const double,const double,const double)> vmec,
          std::function<double(const double,const double,const double)> siesta,
          const std::string &name,
          const double tolarance,
          const double smin,
          const double smax,
          const bool relative) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> min_to_max(smin, smax);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = min_to_max(engine);
        const double u = zero_to_twopi(engine);
        const double v = zero_to_twopi(engine);
        double vmec_value = vmec(s, u, v);
        double siesta_value = siesta(s, u, v);
        double diff = abs(vmec_value - siesta_value);

        if (relative) {
            diff /= std::max(abs(vmec_value), abs(siesta_value));
        }

        const bool pass = diff < tolarance;

        if (!pass) {
            std::cout << name << ": failed at s = " << s
                      <<                    " u = " << u
                      <<                    " v = " << v << std::endl;
            std::cout <<            "vmec value = " << vmec_value
                      <<         " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.001) {
                vmec_value = vmec(s, u, v);
                siesta_value = siesta(s, u, v);
                diff = abs(vmec_value - siesta_value);

                if (relative) {
                    diff /= std::max(abs(vmec_value), abs(siesta_value));
                }

                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << diff << std::endl;
            }

            exit(1);
        }
    }

    std::cout << name << " : passed" << std::endl;
}

//------------------------------------------------------------------------------
///  @brief Test comparitor.
///
///  Randomly samples s,u,v positions to check that the value fall with in the
///  tolerance value.
///
///  @param[in] div       Getter function for the diveregnce value.
///  @param[in] name      Name of the test.
///  @param[in] tolarance Tolarnce threshold that values must fall with in.
///  @param[in] smin      Minimum radial extent.
///  @param[in] smax      Maximum radial extent.
//------------------------------------------------------------------------------
void test(std::function<double(const double,const double,const double)> div,
          const std::string &name,
          const double tolarance,
          const double smin,
          const double smax) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> min_to_max(smin, smax);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = min_to_max(engine);
        const double u = zero_to_twopi(engine);
        const double v = zero_to_twopi(engine);
        double diff = abs(div(s, u, v));

        const bool pass = diff < tolarance;

        if (!pass) {
            std::cout << name << ": failed at s = " << s
                      <<                    " u = " << u
                      <<                    " v = " << v << std::endl;
            std::cout << " divergence = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.00001) {
                diff = abs(div(s, u, v));
                std::cout << s << " "
                          << diff << std::endl;
            }

            exit(1);
        }
    }

    std::cout << name << " : passed" << std::endl;
}

//------------------------------------------------------------------------------
///  @brief Test comparitor.
///
///  Randomly samples s positions to check that the value fall with in the
///  tolerance value.
///
///  @param[in] vmec      Getter function for the vmec value.
///  @param[in] siesta    Getter function for the siesta value.
///  @param[in] name      Name of the test.
///  @param[in] tolarance Tolarnce threshold that values must fall with in.
///  @param[in] smin      Minimum radial extent.
///  @param[in] smax      Maximum radial extent.
///  @param[in] relative  Use relative error.
//------------------------------------------------------------------------------
void test(std::function<double(const double)> vmec,
          std::function<double(const double)> siesta,
          const std::string &name,
          const double tolarance,
          const double smin,
          const double smax,
          const bool relative) {
    std::default_random_engine engine(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> min_to_max(smin, smax);
    std::uniform_real_distribution<double> zero_to_twopi(0.0, 2.0*M_PI);

    for (size_t i = 0; i < 1000000; i++) {
        double s = min_to_max(engine);
        double vmec_value = vmec(s);
        double siesta_value = siesta(s);
        double diff = abs(vmec_value - siesta_value);

        if (relative) {
            diff /= std::max(abs(vmec_value), abs(siesta_value));
        }

        const bool pass = diff < tolarance;

        if (!pass) {
            std::cout << name << ": failed at s = " << std::endl;
            std::cout <<            "vmec value = " << vmec_value
                      <<         " siesta value = " << siesta_value
                      << " difference between vmec and siesta = "
                      << diff << std::endl;

            for (s = 0; s <= 1.0; s += 0.00001) {
                vmec_value = vmec(s);
                siesta_value = siesta(s);
                diff = abs(vmec_value - siesta_value);

                if (relative) {
                    diff /= std::max(abs(vmec_value), abs(siesta_value));
                }

                std::cout << s << " "
                          << vmec_value << " " << siesta_value << " "
                          << diff << std::endl;
            }

            exit(1);
        }
    }

    std::cout << name << " : passed" << std::endl;
}

//------------------------------------------------------------------------------
///  @brief Test comparitor.
///
///  Randomly samples s positions to check that the value fall with in the
///  tolerance value.
///
///  @param[in] vmec      Getter function for the vmec value.
///  @param[in] siesta    Getter function for the siesta value.
///  @param[in] name      Name of the test.
///  @param[in] tolarance Tolarnce threshold that values must fall with in.
///  @param[in] smin      Minimum radial extent.
///  @param[in] smax      Maximum radial extent.
///  @param[in] relative  Use relative error.
//------------------------------------------------------------------------------
void test(std::function<double(void)> vmec,
          std::function<double(void)> siesta,
          const std::string &name,
          const double tolarance,
          const bool relative) {

    double vmec_value = vmec();
    double siesta_value = siesta();
    double diff = abs(vmec_value - siesta_value);

    if (relative) {
        diff /= std::max(abs(vmec_value), abs(siesta_value));
    }

    const bool pass = diff < tolarance;

    if (!pass) {
        std::cout << name << ": failed " << std::endl;
        std::cout <<     "vmec value = " << vmec_value
                  <<  " siesta value = " << siesta_value
                  << " difference between vmec and siesta = "
                  << diff << std::endl;

        exit(1);
    }

    std::cout << name << " : passed" << std::endl;
}

//------------------------------------------------------------------------------
///  @brief Main test program.
///
///  @param[in] argc Number of commandline arguments.
///  @param[in] argv Array of arguments strings.
//------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
    const commandline_parser args(argc, argv);

    const vmec_quantities vmec(args.get<std::string> ("-wout_file"));
    const siesta_quantities siesta(args.get<std::string> ("-restart_file"));

    const std::string test_name = args.get<std::string> ("-test_name");

    if (test_name == "r") {
        RUN_TEST_3D(vmec, siesta, get_r, test_name, args);
    } else if (test_name == "z") {
        RUN_TEST_3D(vmec, siesta, get_z, test_name, args);
    } else if (test_name == "drdu") {
        RUN_TEST_3D(vmec, siesta, get_dr_du, test_name, args);
    } else if (test_name == "drdv") {
        RUN_TEST_3D(vmec, siesta, get_dr_dv, test_name, args);
    } else if (test_name == "dzdu") {
        RUN_TEST_3D(vmec, siesta, get_dz_du, test_name, args);
    } else if (test_name == "dzdv") {
        RUN_TEST_3D(vmec, siesta, get_dz_dv, test_name, args);
    } else if (test_name == "phipf") {
        RUN_TEST_1D(vmec, siesta, get_phipf, test_name, args);
    } else if (test_name == "chipf") {
        RUN_TEST_1D(vmec, siesta, get_chipf, test_name, args);
    } else if (test_name == "bsups") {
        RUN_TEST_3D(vmec, siesta, get_bsups, test_name, args);
    } else if (test_name == "bsupu") {
        RUN_TEST_3D(vmec, siesta, get_bsupu, test_name, args);
    } else if (test_name == "bsupv") {
        RUN_TEST_3D(vmec, siesta, get_bsupv, test_name, args);
    } else if (test_name == "bsubs") {
        RUN_TEST_3D(vmec, siesta, get_bsubs, test_name, args);
    } else if (test_name == "bsubu") {
        RUN_TEST_3D(vmec, siesta, get_bsubu, test_name, args);
    } else if (test_name == "bsubv") {
        RUN_TEST_3D(vmec, siesta, get_bsubv, test_name, args);
    } else if (test_name == "jbsups") {
        RUN_TEST_3D(vmec, siesta, get_jbsups, test_name, args);
    } else if (test_name == "jbsupu") {
        RUN_TEST_3D(vmec, siesta, get_jbsupu, test_name, args);
    } else if (test_name == "jbsupv") {
        RUN_TEST_3D(vmec, siesta, get_jbsupv, test_name, args);
    } else if (test_name == "jksups") {
        RUN_TEST_3D(vmec, siesta, get_jksups, test_name, args);
    } else if (test_name == "jksupu") {
        RUN_TEST_3D(vmec, siesta, get_jksupu, test_name, args);
    } else if (test_name == "jksupv") {
        RUN_TEST_3D(vmec, siesta, get_jksupv, test_name, args);
    } else if (test_name == "jacobian") {
        RUN_TEST_3D(vmec, siesta, get_jacobian, test_name, args);
    } else if (test_name == "curtor") {
        RUN_TEST_0D(vmec, siesta, get_curtor, test_name, args);
    } else if (test_name == "vmec_divb") {
        RUN_TEST_DIV(vmec, get_divb, test_name, args);
    } else if (test_name == "siesta_divb") {
        RUN_TEST_DIV(siesta, get_divb, test_name, args);
    } else if (test_name == "pressure") {
        RUN_TEST_13D(vmec, siesta, get_pressure, test_name, args);
    } else if (test_name == "test") {
        RUN_TEST_1D(vmec, siesta, get_test, test_name, args);
    } else {
        std::cout << "Unknown test name " << test_name << "." << std::endl;
    }
}
