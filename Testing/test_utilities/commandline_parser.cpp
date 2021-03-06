//******************************************************************************
///  @file commandline_parser.cpp
///  @brief Contains implementations to interpolate full and half grid quanities.
//******************************************************************************

#include "commandline_parser.hpp"

//------------------------------------------------------------------------------
///  @brief Get the string value of the agument.
///
///  @param[in] key Commandline key to check.
///  @returns String value of the argument.
//------------------------------------------------------------------------------
template<> std::string commandline_parser::get<std::string> (const std::string &key) const {
    if (!is_set(key)) {
        help();
    }

    return commands.at(key);
}
