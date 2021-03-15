//******************************************************************************
///  @file grid_quantity.cpp
///  @brief Contains implementations to interpolate full and half grid quanities.
//******************************************************************************

#include "grid_quantity.hpp"

//******************************************************************************
///  radial_quantity
//******************************************************************************
//------------------------------------------------------------------------------
///  @brief Construct a radial quantity.
///
///  @param[in] buffer Buffer containing the radial quantity.
//------------------------------------------------------------------------------
radial_quantity::radial_quantity(const std::vector<double> &buffer) :
buffer(buffer), ds(1.0/(buffer.size() - 1.0)) {}

//******************************************************************************
///  full_grid
//******************************************************************************
//------------------------------------------------------------------------------
///  @brief Construct a full grid quantity.
///
///  @param[in] buffer Buffer containing the radial quantity.
//------------------------------------------------------------------------------
full_grid::full_grid(const std::vector<double> &buffer) :
radial_quantity(buffer) {}

//------------------------------------------------------------------------------
///  @brief Get a value at a radial s position.
///
///  @param[in] s Radial s position.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
double full_grid::get(const double s) const {
    if (s >= 1.0) {
        return buffer.back();
    } else if (s <= 0.0) {
        return buffer.front();
    } else {
        const double index = s*(buffer.size() - 1.0);

        const size_t i_low = static_cast<size_t> (floor(index));
        const size_t i_high = i_low + 1;

        const double weight_high = index - i_low;
        const double weight_low = 1.0 - weight_high;

        return weight_high*buffer[i_high] + weight_low*buffer[i_low];
    }
}

//------------------------------------------------------------------------------
///  @brief Get a radial derivative at a s position.
///
///  @returns The radial derivative on the half grid.
//------------------------------------------------------------------------------
half_grid full_grid::get_prime() const {
    std::vector<double> temp(buffer.size());

    for (size_t i = 1, e = buffer.size(); i < e; i++) {
        temp[i] = (buffer[i] - buffer[i - 1])/ds;
    }

    return half_grid(temp);
}

//******************************************************************************
///  half_grid
//******************************************************************************
//------------------------------------------------------------------------------
///  @brief Construct a half grid quantity.
///
///  @param[in] buffer Buffer containing the radial quantity.
//------------------------------------------------------------------------------
half_grid::half_grid(const std::vector<double> &buffer) :
radial_quantity(buffer) {}

//------------------------------------------------------------------------------
///  @brief Get a value at a radial s position.
///
///  @param[in] s Radial s position.
///  @returns The interpolated s position.
//------------------------------------------------------------------------------
double half_grid::get(const double s) const {
    if (s >= 1.0 - ds*0.5) {
        return buffer.back();
    } else if (s <= ds*0.5) {
        return buffer[1];
    } else {
        const double index = (2.0 - 4.0*s - ds*buffer.size() + 2.0*s*buffer.size())
                           / (2.0 - 2.0*ds);

        const size_t i_low = static_cast<size_t> (floor(index));
        const size_t i_high = i_low + 1;

        const double weight_high = index - i_low;
        const double weight_low = 1.0 - weight_high;

        return weight_high*buffer[i_high] + weight_low*buffer[i_low];
    }
}

//------------------------------------------------------------------------------
///  @brief Get a radial derivative at a s position.
///
///  @param[in] m Poloidal mode number.
///  @returns The radial derivative on the half grid.
//------------------------------------------------------------------------------
full_grid half_grid::get_prime(const double m) const {
    std::vector<double> temp(buffer.size());

    for (size_t i = 1, e = buffer.size(); i < e; i++) {
        temp[i] = (buffer[i] - buffer[i - 1])/ds;
    }
    temp[temp.size() - 1] = temp[temp.size() - 2];
    if (m == 1.0) {
        temp[0] = temp[1];
    } else {
        temp[0] = 0.0;
    }

    return full_grid(temp);
}
