/// Copyright (c) 2025 - present KMX Systems. All rights reserved.
/// @file gis/coordinate/base.hpp
/// @brief Defines core constants, data structures, and conversion functions for GIS operations.
///
/// This header provides fundamental building blocks for Geographic Information Systems (GIS)
/// calculations, including template-based constants, ellipsoid definitions, coordinate
/// structures (geodetic, geocentric/XYZ), Helmert transformation parameters, and functions
/// for coordinate conversions between different systems (WGS84, Stereo70, Geocentric).
/// It utilizes C++20 features like `<numbers>` and concepts implicitly via `static_assert`.
///
/// @note This file assumes a C++20 compliant compiler.
/// @warning Does not use PCH if PCH is not defined.
#pragma once
#ifndef PCH
    #include <algorithm> ///< For std::clamp
    #include <cmath> ///< For std::sqrt, std::sin, std::cos, std::atan, std::atan2, std::hypot, std::abs, std::copysign, std::pow, std::isfinite, std::fmod
    #include <iomanip>     ///< For std::fixed, std::setprecision
    #include <limits>      ///< For std::numeric_limits
    #include <numbers>     ///< Requires C++20, for std::numbers::pi_v
    #include <stdexcept>   ///< For std::invalid_argument, std::runtime_error
    #include <string>      ///< For std::to_string in error messages
    #include <tuple>       ///< For std::tuple
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

    /// @brief Represents a geodetic coordinate (latitude, longitude, ellipsoidal altitude).
    ///
    /// Coordinates are stored in degrees for latitude and longitude, and meters for altitude.
    /// Latitude ranges from -90 to +90 degrees. Longitude is typically normalized to -180 to +180 degrees.
    /// Altitude is relative to the surface of the reference ellipsoid.
    /// Coordinates are mutable, allowing for normalization or modification after calculation.
    ///
    /// @tparam T The floating-point type (e.g., `float`, `double`, `long double`) for the coordinate values.
    ///           Must satisfy `std::is_floating_point_v<T>`.
    template <typename T>
    struct geodetic_coord
    {
        /// @brief Static assertion to ensure T is a floating-point type.
        static_assert(std::is_floating_point_v<T>);
        /// @brief Alias for constants of the corresponding floating-point type.
        using K = constants<T>;
        /// @brief Alias for floating-point value type.
        using value_type = T;

        // Coordinates are mutable (e.g., result of calculation, can be normalized)
        /// @brief Geodetic latitude in degrees. Expected range [-90, 90]. Positive north.
        value_type latitude {};
        /// @brief Geodetic longitude in degrees. Expected range [-180, 180]. Positive east.
        value_type longitude {};
        /// @brief Ellipsoidal altitude in meters. Altitude above the reference ellipsoid surface.
        value_type altitude {};

        /// @brief Normalizes the latitude and longitude coordinates to standard ranges.
        ///
        /// Clamps latitude to [-90, 90] degrees.
        /// Wraps longitude to the range (-180, 180] degrees.
        /// Modifies the object's state.
        /// @note Uses `std::clamp` for latitude and loops for longitude wrapping.
        constexpr void normalize() noexcept
        {
            // Clamp latitude to the valid range [-90, 90]
            latitude = std::clamp(latitude, -K::values::ninety, K::values::ninety);

            // Normalize longitude to the range (-180, 180]
            // Use fmod for potentially better performance on some platforms, handling edge cases.
            if (!std::isfinite(longitude))
                return; // Avoid fmod with non-finite values

            longitude = std::fmod(longitude, K::values::three_sixty); // -> (-360, 360)

            if (longitude > K::values::one_eighty)
                longitude -= K::values::three_sixty;      // -> (-180, 180]
            else if (longitude <= -K::values::one_eighty) // Catches <= -180
                longitude += K::values::three_sixty;      // -> (-180, 180] (maps -180 to +180)
            // Result is in (-180, 180]
        }

        /// @brief Validates the geodetic coordinates.
        ///
        /// Checks if latitude, longitude, and altitude are finite numbers.
        /// Checks if latitude is within the valid range [-90, 90] degrees.
        /// Does not check longitude range as it can be normalized.
        /// Does not modify the object's state.
        ///
        /// @throws std::invalid_argument if any coordinate is not finite.
        /// @throws std::invalid_argument if latitude is outside the [-90, 90] degree range.
        constexpr void validate() const noexcept(false)
        {
            using std::isfinite;
            if (!isfinite(latitude) || !isfinite(longitude) || !isfinite(altitude))
                throw std::invalid_argument("Geodetic coordinates (latitude, longitude, altitude) must be finite");

            // Use a small tolerance for latitude range check if needed due to floating point inaccuracies?
            // Sticking to strict comparison for now. Add epsilon if required.
            constexpr value_type lat_tolerance = K::tolerance::near_zero; // Use a small tolerance
            if (latitude < (-K::values::ninety - lat_tolerance) || latitude > (K::values::ninety + lat_tolerance))
                throw std::invalid_argument("Latitude must be between -90 and 90 degrees (got " + std::to_string(latitude) + ")");
            // Longitude validation might be needed depending on context, but normalize() handles wrapping.
            // Altitude validation (e.g., reasonableness) is context-dependent and not done here.
        }

        std::ostream& nice_print(std::ostream& out) const
        {
            using std::setw;
            out << std::fixed << std::setprecision(8);
            out << "[lat: " << setw(12) << latitude << " deg, lon: " << setw(12) << longitude << " deg, alt: " << std::setprecision(4)
                << setw(12) << altitude << " m]";
            return out;
        }

        std::ostream& nice_println(std::ostream& out) const { return nice_print(out) << std::endl; }

        std::ostream& print(std::ostream& out) const
        {
            using std::setw;
            out << std::fixed << std::setprecision(8);
            out << latitude << ", " << longitude << ", " << altitude;
            return out;
        }

        std::ostream& println(std::ostream& out) const { return print(out) << std::endl; }
    };

    /// @brief Represents a 3D Cartesian coordinate (X, Y, Z).
    ///
    /// This structure is used for various Cartesian systems, such as Earth-Centered, Earth-Fixed (ECEF)
    /// geocentric coordinates or projected coordinates (often using X for Easting, Y for Northing,
    /// and Z optionally for altitude or elevation, **but axis meaning can vary by CRS definition**).
    /// Coordinates are typically in meters. Coordinates are mutable.
    ///
    /// @tparam T The floating-point type (e.g., `float`, `double`, `long double`) for the coordinate values.
    ///           Must satisfy `std::is_floating_point_v<T>`.
    template <typename T>
    struct xyz_coord
    {
        /// @brief Static assertion to ensure T is a floating-point type.
        static_assert(std::is_floating_point_v<T>);

        using value_type = T;

        value_type x {}; ///< X coordinate value (e.g., meters in ECEF; Northing in EPSG:31700).
        value_type y {}; ///< Y coordinate value (e.g., meters in ECEF; Easting in EPSG:31700).
        value_type z {}; ///< Z coordinate value (e.g., meters in ECEF or Altitude/Elevation).

        /// @brief Validates that the X, Y, Z coordinates are finite numbers.
        ///
        /// Checks if x, y, and z are finite. Does not modify the object's state.
        ///
        /// @throws std::invalid_argument if any coordinate (x, y, z) is not finite (NaN or infinity).
        constexpr void validate() const noexcept(false)
        {
            using std::isfinite;
            if (!isfinite(x) || !isfinite(y) || !isfinite(z))
                throw std::invalid_argument("XYZ coordinates (x, y, z) must be finite");
        }

        std::ostream& print(std::ostream& out) const
        {
            using std::setw;
            out << std::fixed << std::setprecision(5);
            out << x << ", " << y << ", " << z;
            return out;
        }

        std::ostream& println(std::ostream& out) const { return print(out) << std::endl; }
    };

    /// @brief Calculates the absolute difference between two Cartesian coordinates.
    ///
    /// Computes the absolute difference for each corresponding component (x, y, z)
    /// of two xyz_coord objects.
    ///
    /// @tparam T The underlying numeric type of the coordinates (e.g., float, double).
    /// @param a The first Cartesian coordinate object.
    /// @param b The second Cartesian coordinate object.
    /// @return A std::tuple containing the absolute differences in the order:
    ///         (|a.x - b.x|, |a.y - b.y|, |a.z - b.z|).
    /// @note This function is guaranteed not to throw exceptions.
    template <typename T>
    std::tuple<T, T, T> diff(const xyz_coord<T>& a, const xyz_coord<T>& b) noexcept
    {
        using std::abs;
        return {abs(a.x - b.x), abs(a.y - b.y), abs(a.z - b.z)};
    }

    /// @brief Calculates the absolute difference between two geodetic coordinates.
    ///
    /// Computes the absolute difference for each corresponding component (latitude,
    /// longitude, altitude) of two geodetic_coord objects.
    ///
    /// @tparam T The underlying numeric type of the coordinates (e.g., float, double).
    /// @param a The first geodetic coordinate object.
    /// @param b The second geodetic coordinate object.
    /// @return A std::tuple containing the absolute differences in the order:
    ///         (|a.latitude - b.latitude|, |a.longitude - b.longitude|, |a.altitude - b.altitude|).
    /// @note This function is guaranteed not to throw exceptions.
    /// @warning This performs a simple arithmetic difference. For latitude/longitude,
    ///          this might not represent the true shortest angular distance, especially
    ///          near the poles or the +/-180 degree meridian. Consider using
    ///          specialized geodetic distance functions if needed.
    template <typename T>
    std::tuple<T, T, T> diff(const geodetic_coord<T>& a, const geodetic_coord<T>& b) noexcept
    {
        using std::abs;
        return {abs(a.latitude - b.latitude), abs(a.longitude - b.longitude), abs(a.altitude - b.altitude)};
    }

    /// @brief Alias for `xyz_coord<T>` representing geocentric (ECEF) coordinates.
    /// @tparam T The floating-point type for coordinate values.
    template <typename T>
    using geocentric_coord = xyz_coord<T>;

    /// @brief Alias for `xyz_coord<T>` representing projected (map) coordinates.
    ///
    /// Note: The meaning of X and Y depends on the specific Projected Coordinate Reference System (CRS).
    /// For EPSG:31700 (Stereo 70), X represents Northing and Y represents Easting.
    /// Z typically represents altitude/elevation above a reference surface.
    /// @tparam T The floating-point type for coordinate values.
    template <typename T>
    using projected_coord = xyz_coord<T>;

    /// @namespace kmx::gis::wgs84
    /// @brief Contains types and parameters specific to the WGS84 datum and coordinate system.
    namespace wgs84
    {
        /// @brief Alias for `geodetic_coord<T>` representing coordinates in the WGS84 datum.
        /// @tparam T The floating-point type for coordinate values.
        template <typename T>
        using coordinate = geodetic_coord<T>;
    } // namespace wgs84

} // namespace kmx::gis
