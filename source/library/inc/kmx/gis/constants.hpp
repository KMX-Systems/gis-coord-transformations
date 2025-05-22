/// Copyright (c) 2025 - present KMX Systems. All rights reserved.
/// @file gis/constants.hpp
/// @brief Defines core constants for GIS operations.
///
/// @note This file assumes a C++20 compliant compiler.
/// @warning Does not use PCH if PCH is not defined.
#pragma once
#ifndef PCH
    #include <cmath> ///< For std::sqrt, std::sin, std::cos, std::atan, std::atan2, std::hypot, std::abs, std::copysign, std::pow, std::isfinite, std::fmod
    #include <limits>      ///< For std::numeric_limits
    #include <numbers>     ///< Requires C++20, for std::numbers::pi_v
    #include <type_traits> ///< For std::is_floating_point_v
#endif

/// @namespace kmx::gis
/// @brief Namespace for Geographic Information System (GIS) related functionalities.
///
/// This namespace contains structures, constants, and algorithms commonly used
/// in geodetic and cartographic computations.
namespace kmx::gis
{
    /// @brief Provides compile-time constants for GIS calculations, templated on floating-point type.
    ///
    /// This struct centralizes commonly used numerical values, mathematical constants,
    /// conversion factors, tolerances, and series coefficients to ensure consistency
    /// and precision across different GIS operations. It is organized into nested structs
    /// for clarity.
    ///
    /// @tparam T The floating-point type (e.g., `float`, `double`, `long double`) for the constants.
    ///           Must satisfy `std::is_floating_point_v<T>`.
    template <typename T>
    struct constants
    {
        /// @brief Static assertion to ensure T is a floating-point type.
        static_assert(std::is_floating_point_v<T>, "constants requires a floating-point type.");
        /// @brief Deleted default constructor to prevent instantiation.
        constants() = delete;

        /// @brief Holds common numerical values used in calculations.
        struct values
        {
            /// @brief Deleted default constructor to prevent instantiation.
            values() = delete;

            static constexpr T zero = T(0);            ///< Value 0.
            static constexpr T half = T(0.5);          ///< Value 0.5.
            static constexpr T one = T(1);             ///< Value 1.
            static constexpr T two = T(2);             ///< Value 2.
            static constexpr T four = T(4);            ///< Value 4.
            static constexpr T five = T(5);            ///< Value 5.
            static constexpr T six = T(6);             ///< Value 6.
            static constexpr T seven = T(7);           ///< Value 7.
            static constexpr T eight = T(8);           ///< Value 8.
            static constexpr T ten = T(10);            ///< Value 10.
            static constexpr T twelve = T(12.0);       ///< Value 12.
            static constexpr T thirteen = T(13.0);     ///< Value 13.
            static constexpr T fifteen = T(15.0);      ///< Value 15.
            static constexpr T twenty_four = T(24.0);  ///< Value 24.
            static constexpr T twenty_nine = T(29.0);  ///< Value 29.
            static constexpr T forty_eight = T(48.0);  ///< Value 48.
            static constexpr T eighty_one = T(81.0);   ///< Value 81.
            static constexpr T ninety = T(90.0);       ///< Value 90 (degrees).
            static constexpr T hundred = T(100.0);     ///< Value 100.
            static constexpr T one_twenty = T(120.0);  ///< Value 120.
            static constexpr T one_eighty = T(180.0);  ///< Value 180 (degrees).
            static constexpr T two_forty = T(240.0);   ///< Value 240.
            static constexpr T three_sixty = T(360.0); ///< Value 360 (degrees).
            static constexpr T million = T(1e6);       ///< Value 1,000,000.
        };

        /// @brief Holds core mathematical constants.
        struct math
        {
            /// @brief Deleted default constructor to prevent instantiation.
            math() = delete;

            /// @brief The mathematical constant Pi (π), using C++20 `std::numbers`.
            static constexpr T pi = std::numbers::pi_v<T>;
            /// @brief Two times Pi (2π).
            static constexpr T two_pi = values::two * pi;
            /// @brief Pi divided by two (π/2).
            static constexpr T half_pi = pi / values::two;
            /// @brief Pi divided by four (π/4).
            static constexpr T quarter_pi = pi / values::four;
            /// @brief Machine epsilon for the floating-point type T. Represents the smallest value such that 1.0 + epsilon != 1.0.
            static constexpr T epsilon = std::numeric_limits<T>::epsilon();
        };

        /// @brief Holds factors for converting between different units.
        struct conversions
        {
            /// @brief Deleted default constructor to prevent instantiation.
            conversions() = delete;

            /// @brief Factor to convert degrees to arcseconds (3600).
            static constexpr T sec_per_degree = T(3600);
            /// @brief Factor to convert degrees to radians (π / 180).
            static constexpr T deg_to_rad = math::pi / values::one_eighty;
            /// @brief Factor to convert radians to degrees (180 / π).
            static constexpr T rad_to_deg = values::one_eighty / math::pi;
            /// @brief Factor to convert arcseconds to radians (deg_to_rad / sec_per_degree).
            static constexpr T arcsec_to_rad = deg_to_rad / sec_per_degree;
            /// @brief Factor to convert parts per million (ppm) to a dimensionless scale factor (1 / 1,000,000).
            static constexpr T ppm_to_scale = values::one / values::million;
        };

        /// @brief Holds tolerance values for comparisons and iteration limits.
        struct tolerance
        {
            /// @brief Deleted default constructor to prevent instantiation.
            tolerance() = delete;

            /// @brief Default tolerance for iterative geodetic calculations (e.g., geocentric to geodetic).
            static constexpr T default_geodetic = T(1e-11);
            /// @brief A small value near zero, based on machine epsilon, used for singularity checks.
            static constexpr T near_zero = math::epsilon * values::hundred;
            /// @brief Maximum number of iterations for iterative algorithms (e.g., geocentric to geodetic).
            static constexpr int max_iterations = 10;
        };

        /// @brief Holds coefficients used in series expansions (e.g., for projection calculations).
        struct series
        {
            /// @brief Deleted default constructor to prevent instantiation.
            series() = delete;

            // C2 Coefficients (used in geodetic latitude series)
            /// @brief Numerator of the first term for C2 coefficient calculation.
            static constexpr T c2_term1_num = values::seven;
            /// @brief Denominator of the first term for C2 coefficient calculation.
            static constexpr T c2_term1_den = values::forty_eight;
            /// @brief Numerator of the second term for C2 coefficient calculation.
            static constexpr T c2_term2_num = values::twenty_nine;
            /// @brief Denominator of the second term for C2 coefficient calculation.
            static constexpr T c2_term2_den = values::two_forty;
            /// @brief Numerator of the third term for C2 coefficient calculation.
            static constexpr T c2_term3_num = T(811);
            /// @brief Denominator of the third term for C2 coefficient calculation.
            static constexpr T c2_term3_den = T(11520);
            // C3 Coefficients (used in geodetic latitude series)
            /// @brief Numerator of the first term for C3 coefficient calculation.
            static constexpr T c3_term1_num = values::seven;
            /// @brief Denominator of the first term for C3 coefficient calculation.
            static constexpr T c3_term1_den = values::one_twenty;
            /// @brief Numerator of the second term for C3 coefficient calculation.
            static constexpr T c3_term2_num = values::eighty_one;
            /// @brief Denominator of the second term for C3 coefficient calculation.
            static constexpr T c3_term2_den = T(1120);
            // C4 Coefficients (used in geodetic latitude series)
            /// @brief Numerator of the first term for C4 coefficient calculation.
            static constexpr T c4_term1_num = T(4279);
            /// @brief Denominator of the first term for C4 coefficient calculation.
            static constexpr T c4_term1_den = T(161280);
        };
    };
} // namespace gis
