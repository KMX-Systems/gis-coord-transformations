/// Copyright (c) 2025 - present KMX Systems. All rights reserved.
/// @file gis/coordiante/stereo70.hpp
/// @brief Defines Stereo70 data structures and functions for GIS operations.
///
/// @note This file assumes a C++20 compliant compiler.
/// @warning Does not use PCH if PCH is not defined.
#pragma once
#ifndef PCH
    #include <kmx/gis/coordinate/base.hpp>
    #include <string> ///< For std::to_string in error messages
#endif

/// @namespace kmx::gis::stereo70
/// @brief Contains types and parameters specific to the Romanian Stereo70 projection system (EPSG:31700).
///
/// Stereo70 is based on the Krasovsky 1940 ellipsoid and uses an Oblique Stereographic projection.
/// The coordinate system definition (EPSG:31700) uses X for Northing and Y for Easting.
namespace kmx::gis::stereo70
{
    /// @brief Alias for `projected_coord<T>` representing coordinates in the Stereo70 system.
    /// @note For Stereo70 (EPSG:31700), the `x` member holds Northing and the `y` member holds Easting.
    /// @tparam T The floating-point type for coordinate values.
    template <typename T>
    struct coordinate: gis::coordinate::projected_coord<T>
    {
        using base_t = gis::coordinate::projected_coord<T>;

        constexpr coordinate() = default;
        constexpr coordinate(const coordinate&) = default;
        constexpr coordinate(coordinate&&) = default;

        coordinate& operator=(const coordinate&) = default;
        coordinate& operator=(coordinate&&) = default;

        constexpr coordinate(const T x, const T y, const T z) noexcept: base_t {x, y, z} // Initialize base using aggregate init
        {
        }

        using base_t::x;
        using base_t::y;
        using base_t::z;

        constexpr void validate() const noexcept(false)
        {
            /// @brief Defines the approximate valid range for Stereo70 Easting coordinates within Romania.
            constexpr T min_easting = T(110000.0);
            constexpr T max_easting = T(890000.0);

            /// @brief Defines the approximate valid range for Stereo70 Northing coordinates within Romania.
            constexpr T min_northing = T(220000.0);
            constexpr T max_northing = T(780000.0);

            if ((y < min_easting) || (y > max_easting))
                throw std::runtime_error("invalid easting value " + std::to_string(y));

            if ((x < min_northing) && (x > max_northing))
                throw std::runtime_error("invalid northing value " + std::to_string(x));
        }

        std::ostream& nice_print(std::ostream& out) const
        {
            using std::setw;
            out << std::fixed << std::setprecision(3);
            out << "[X: " << setw(12) << x << " m, Y: " << setw(12) << y << " m, Z: " << setw(9) << z << " m]";
            return out;
        }

        std::ostream& nice_println(std::ostream& out) const { return nice_print(out) << std::endl; }

        std::ostream& print(std::ostream& out) const
        {
            using std::setw;
            out << std::fixed << std::setprecision(3);
            out << x << ", " << y << ", " << z;
            // out << y << ", " << x;
            return out;
        }

        std::ostream& println(std::ostream& out) const { return print(out) << std::endl; }
    };

    template <typename Coord>
    std::ostream& operator<<(std::ostream& out, const Coord& coord)
    {
        return coord.print(out);
    }

    /// @brief Defines the core projection parameters for the Stereo70 system (EPSG:31700).
    ///
    /// These are the fundamental constants defining the Stereo70 projection,
    /// such as the latitude/longitude of origin, scale factor, and false easting/northing.
    ///
    /// @tparam T The floating-point type (e.g., `float`, `double`, `long double`) for the parameters.
    ///           Must satisfy `std::is_floating_point_v<T>`.
    template <typename T>
    struct projection_params
    {
        /// @brief Static assertion to ensure T is a floating-point type.
        static_assert(std::is_floating_point_v<T>);
        /// @brief Latitude of the projection origin (degrees). EPSG:8801
        static constexpr T lat0_deg = T(46.0);
        /// @brief Longitude of the projection origin / Central Meridian (degrees). EPSG:8802
        static constexpr T lon0_deg = T(25.0);
        /// @brief Scale factor at the projection origin. EPSG:8805
        static constexpr T k0 = T(0.99975);
        /// @brief False Easting applied to the Y coordinate (meters). EPSG:8806
        static constexpr T fe = T(500000.0);
        /// @brief False Northing applied to the X coordinate (meters). EPSG:8807
        static constexpr T fn = T(500000.0);

        // Tolerances specific to Stereo70 calculations
        /// @brief Tolerance for checking if a point is at the projection origin (degrees).
        static constexpr T origin_tol_deg = T(1e-9);
        /// @brief Tolerance for checking proximity to the antipodal point of the projection origin.
        static constexpr T antipodal_tol = T(1e-10);
        /// @brief Tolerance for checking if a projected coordinate is at the origin (meters).
        static constexpr T origin_proj_tol = T(1e-6);
    };
}
