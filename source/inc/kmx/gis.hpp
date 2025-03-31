// Copyright (c) 2025 - present KMX Systems. All rights reserved.
/// @file gis_core.hpp
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

#ifndef PCH
    #include <algorithm> ///< For std::clamp
    #include <cmath> ///< For std::sqrt, std::sin, std::cos, std::atan, std::atan2, std::hypot, std::abs, std::copysign, std::pow, std::isfinite
    #include <limits>      ///< For std::numeric_limits
    #include <numbers>     ///< Requires C++20, for std::numbers::pi_v
    #include <stdexcept>   ///< For std::invalid_argument, std::runtime_error
    #include <string>      ///< For std::to_string in error messages
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

    /// @brief Represents a reference ellipsoid model.
    ///
    /// Stores the defining parameters of an ellipsoid (semi-major axis and inverse flattening)
    /// and calculates derived parameters like flattening, first and second eccentricity squared,
    /// first eccentricity, and semi-minor axis. These parameters are conceptually constant
    /// after the ellipsoid object is constructed.
    ///
    /// @tparam T The floating-point type (e.g., `float`, `double`, `long double`) for the ellipsoid parameters.
    ///           Must satisfy `std::is_floating_point_v<T>`.
    template <typename T>
    struct ellipsoid
    {
        /// @brief Static assertion to ensure T is a floating-point type.
        static_assert(std::is_floating_point_v<T>, "ellipsoid requires a floating-point type.");
        /// @brief Alias for constants of the corresponding floating-point type.
        using K = constants<T>;

        // Members are conceptually constant after construction
        /// @brief Semi-major axis (equatorial radius) in meters. Must be positive.
        const T a;
        /// @brief Inverse flattening (1/f). 0 for a perfect sphere. Cannot be negative.
        const T inv_f;
        /// @brief Flattening ((a-b)/a). Derived from `inv_f`.
        const T f;
        /// @brief First eccentricity squared (e^2 = (a^2 - b^2) / a^2 = 2f - f^2). Derived from `f`.
        const T e2;
        /// @brief First eccentricity (sqrt(e2)). Derived from `e2`.
        const T e;
        /// @brief Semi-minor axis (polar radius) in meters (b = a * (1 - f)). Derived from `a` and `f`.
        const T b;
        /// @brief Second eccentricity squared (e'^2 = (a^2 - b^2) / b^2 = e2 / (1 - e2)). Derived from `e2`. Can be infinity if e2
        /// == 1.
        const T ep2;

        /// @brief Constructs an ellipsoid from its semi-major axis and inverse flattening.
        ///
        /// Calculates all derived parameters upon construction. Throws an exception if the
        /// input parameters are invalid (e.g., non-positive semi-major axis, negative inverse flattening,
        /// or resulting parameters leading to mathematical instability like 1-e^2 <= 0 unless it's exactly 1).
        ///
        /// @param semi_major The semi-major axis (a) in meters. Must be positive.
        /// @param inverse_flattening The inverse flattening (1/f). Use 0 for a sphere. Must not be negative.
        /// @throws std::invalid_argument if `semi_major` is not positive.
        /// @throws std::invalid_argument if `inverse_flattening` is negative.
        /// @throws std::invalid_argument if derived parameters result in an invalid state (e.g., 1-e^2 <= 0 and e^2 != 1).
        constexpr ellipsoid(const T semi_major, const T inverse_flattening) noexcept(false):
            a(semi_major),
            inv_f(inverse_flattening),
            f(inv_f == K::values::zero ? K::values::zero : K::values::one / inv_f),  // Handle sphere case (inv_f = 0)
            e2(f == K::values::zero ? K::values::zero : K::values::two * f - f * f), // Handle sphere case (f = 0)
            e(std::sqrt(e2)),
            b(a * (K::values::one - f)),
            // Handle potential division by zero or invalid state for ep2 calculation
            ep2(e2 == K::values::one ? std::numeric_limits<T>::infinity() // Case where e=1 (degenerate ellipsoid)
                                       :
                                       (K::values::one - e2 <= K::tolerance::near_zero ? // Case where 1-e2 is near or below zero
                                            std::numeric_limits<T>::quiet_NaN()          // Or throw? NaN signals invalid intermediate state
                                            :
                                            e2 / (K::values::one - e2)) // Normal case
            )
        {
            // Input validation
            if (a <= K::values::zero)
                throw std::invalid_argument("Semi-major axis must be positive (got " + std::to_string(a) + ")");

            if (inv_f < K::values::zero)
                throw std::invalid_argument("Inverse flattening cannot be negative (got " + std::to_string(inv_f) + ")");

            // Derived parameter validation (check 1-e^2 explicitly)
            // Allow e2=1 (degenerate ellipsoid), but disallow 1-e2 <= 0 otherwise.
            if ((K::values::one - e2) <= K::values::zero && e2 != K::values::one)
                throw std::invalid_argument(
                    "Invalid ellipsoid parameters leading to non-positive (1-e^2) (value: " + std::to_string(K::values::one - e2) + ")");

            // Optional: Check if ep2 became NaN due to near-zero denominator in the calculation above
            if (!std::isfinite(ep2) && e2 != K::values::one) // Check finiteness unless e2=1 where infinity is expected
                throw std::invalid_argument("Invalid ellipsoid parameters leading to non-finite second eccentricity squared (ep2)");
        }

        /// @brief Checks if the ellipsoid represents a perfect sphere.
        /// @return `true` if the first eccentricity squared is zero (or numerically close enough), `false` otherwise.
        [[nodiscard]] constexpr bool is_sphere() const noexcept
        {
            // Check against zero or a very small tolerance related to epsilon if strict zero comparison is too fragile.
            // Using strict zero here as e2 is calculated directly.
            return e2 == K::values::zero;
        }
    };

    /// @brief Represents the parameters for a 7-parameter Helmert transformation (coordinate frame method).
    ///
    /// This structure holds the translation offsets (dx, dy, dz), rotation angles (rx, ry, rz),
    /// and scale factor difference (ds) used to transform coordinates between two geocentric
    /// Cartesian systems (e.g., WGS84 ECEF to Pulkovo 1942 ECEF).
    /// Parameters are typically mutable as they might be loaded from external sources.
    /// Includes methods to get rotation angles in radians and the absolute scale factor.
    ///
    /// @tparam T The floating-point type (e.g., `float`, `double`, `long double`) for the parameters.
    ///           Must satisfy `std::is_floating_point_v<T>`.
    /// @note The rotation convention used here is the Coordinate Frame Rotation. When transforming
    ///       from source to target: P_target = T + (1+ds) * R * P_source, where R contains small angle
    ///       approximations for rotations around X, Y, Z axes. The reverse transformation involves
    ///       inverting this process.
    template <typename T>
    struct helmert_params
    {
        /// @brief Static assertion to ensure T is a floating-point type.
        static_assert(std::is_floating_point_v<T>, "helmert_params requires a floating-point type.");
        /// @brief Alias for constants of the corresponding floating-point type.
        using K = constants<T>;

        // These parameters are mutable by design (e.g., loaded from external source)
        /// @brief Translation offset along the X-axis (meters).
        T dx {};
        /// @brief Translation offset along the Y-axis (meters).
        T dy {};
        /// @brief Translation offset along the Z-axis (meters).
        T dz {};
        /// @brief Rotation angle around the X-axis (arcseconds).
        T rx_sec {};
        /// @brief Rotation angle around the Y-axis (arcseconds).
        T ry_sec {};
        /// @brief Rotation angle around the Z-axis (arcseconds).
        T rz_sec {};
        /// @brief Scale factor difference (Scale = 1 + ds_ppm * 1e-6) in parts per million (ppm).
        T ds_ppm {};

        /// @brief Gets the X-axis rotation angle in radians.
        /// @return Rotation angle in radians.
        [[nodiscard]] constexpr T rx_rad() const noexcept { return rx_sec * K::conversions::arcsec_to_rad; }
        /// @brief Gets the Y-axis rotation angle in radians.
        /// @return Rotation angle in radians.
        [[nodiscard]] constexpr T ry_rad() const noexcept { return ry_sec * K::conversions::arcsec_to_rad; }
        /// @brief Gets the Z-axis rotation angle in radians.
        /// @return Rotation angle in radians.
        [[nodiscard]] constexpr T rz_rad() const noexcept { return rz_sec * K::conversions::arcsec_to_rad; }
        /// @brief Calculates the absolute scale factor (1 + ds).
        /// @return The dimensionless scale factor.
        [[nodiscard]] constexpr T scale_factor() const noexcept { return K::values::one + ds_ppm * K::conversions::ppm_to_scale; }

        /// @brief Validates that all Helmert parameters are finite numbers.
        ///
        /// Checks if translations, rotations (in arcseconds), and scale difference (in ppm)
        /// are finite. Does not modify the object's state.
        ///
        /// @throws std::invalid_argument if any parameter is not finite (NaN or infinity).
        constexpr void validate() const noexcept(false)
        {
            using std::isfinite;
            if (!isfinite(dx) || !isfinite(dy) || !isfinite(dz))
                throw std::invalid_argument("Helmert translation parameters (dx, dy, dz) must be finite");

            if (!isfinite(rx_sec) || !isfinite(ry_sec) || !isfinite(rz_sec))
                throw std::invalid_argument("Helmert rotation parameters (rx_sec, ry_sec, rz_sec) must be finite");

            if (!isfinite(ds_ppm))
                throw std::invalid_argument("Helmert scale parameter (ds_ppm) must be finite");

            // Optionally, validate the calculated scale factor is usable (e.g., not near zero if division is needed)
            if (std::abs(scale_factor()) < K::tolerance::near_zero)
                throw std::invalid_argument("Helmert scale factor is dangerously close to zero.");
        }
    };

    /// @brief Represents a geodetic coordinate (latitude, longitude, ellipsoidal height).
    ///
    /// Coordinates are stored in degrees for latitude and longitude, and meters for height.
    /// Latitude ranges from -90 to +90 degrees. Longitude is typically normalized to -180 to +180 degrees.
    /// Height is relative to the surface of the reference ellipsoid.
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

        // Coordinates are mutable (e.g., result of calculation, can be normalized)
        /// @brief Geodetic latitude in degrees. Expected range [-90, 90]. Positive north.
        T latitude {};
        /// @brief Geodetic longitude in degrees. Expected range [-180, 180]. Positive east.
        T longitude {};
        /// @brief Ellipsoidal height in meters. Height above the reference ellipsoid surface.
        T height {};

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
            // fmod approach: longitude = std::fmod(longitude + K::values::one_eighty, K::values::three_sixty) - K::values::one_eighty;
            // Handle the case where result is exactly -180, map it to +180? The loop does (-180, 180].
            while (longitude > K::values::one_eighty)
                longitude -= K::values::three_sixty;
            while (longitude <= -K::values::one_eighty) // Use <= to map -180 to range boundary handled by next loop
                longitude += K::values::three_sixty;
            // The interval is now (-180, 180]. If longitude was exactly -180, it becomes +180.
            // Let's refine to ensure (-180, 180] is strictly enforced if needed, or [-180, 180)
            if (longitude == -K::values::one_eighty)
            {
                // Some conventions prefer [-180, 180), others (-180, 180].
                // The loop above results in (-180, 180]. Stick with that unless required otherwise.
                // If [-180, 180) is needed:
                // longitude = std::fmod(std::fmod(longitude, K::values::three_sixty) + K::values::three_sixty, K::values::three_sixty);
                // // -> [0, 360) if (longitude >= K::values::one_eighty) longitude -= K::values::three_sixty; // -> [-180, 180)
            }
        }

        /// @brief Validates the geodetic coordinates.
        ///
        /// Checks if latitude, longitude, and height are finite numbers.
        /// Checks if latitude is within the valid range [-90, 90] degrees.
        /// Does not check longitude range as it can be normalized.
        /// Does not modify the object's state.
        ///
        /// @throws std::invalid_argument if any coordinate is not finite.
        /// @throws std::invalid_argument if latitude is outside the [-90, 90] degree range.
        constexpr void validate() const noexcept(false)
        {
            if (!std::isfinite(latitude) || !std::isfinite(longitude) || !std::isfinite(height))
                throw std::invalid_argument("Geodetic coordinates (latitude, longitude, height) must be finite");
            // Use a small tolerance for latitude range check if needed due to floating point inaccuracies?
            // Sticking to strict comparison for now.
            if (latitude < -K::values::ninety || latitude > K::values::ninety)
                throw std::invalid_argument("Latitude must be between -90 and 90 degrees (got " + std::to_string(latitude) + ")");
            // Longitude validation might be needed depending on context, but normalize() handles wrapping.
            // Height validation (e.g., reasonableness) is context-dependent and not done here.
        }
    };

    /// @brief Represents a 3D Cartesian coordinate (X, Y, Z).
    ///
    /// This structure is used for various Cartesian systems, such as Earth-Centered, Earth-Fixed (ECEF)
    /// geocentric coordinates or projected coordinates (often using X for Easting, Y for Northing,
    /// and Z optionally for height or elevation). Coordinates are typically in meters.
    /// Coordinates are mutable.
    ///
    /// @tparam T The floating-point type (e.g., `float`, `double`, `long double`) for the coordinate values.
    ///           Must satisfy `std::is_floating_point_v<T>`.
    template <typename T>
    struct xyz_coord
    {
        /// @brief Static assertion to ensure T is a floating-point type.
        static_assert(std::is_floating_point_v<T>);

        T x {}; ///< X coordinate value (e.g., meters in ECEF or Easting in projected).
        T y {}; ///< Y coordinate value (e.g., meters in ECEF or Northing in projected).
        T z {}; ///< Z coordinate value (e.g., meters in ECEF or Height/Elevation in projected).

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
    };

    /// @brief Alias for `xyz_coord<T>` representing geocentric (ECEF) coordinates.
    /// @tparam T The floating-point type for coordinate values.
    template <typename T>
    using geocentric_coord = xyz_coord<T>;

    /// @brief Alias for `xyz_coord<T>` representing projected (map) coordinates.
    ///
    /// Typically, X maps to Easting, Y maps to Northing, and Z represents height/elevation
    /// above a reference surface (often the geoid or ellipsoid, depending on the projection definition).
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

    /// @namespace kmx::gis::stereo70
    /// @brief Contains types and parameters specific to the Romanian Stereo70 projection system.
    ///
    /// Stereo70 is based on the Krasovsky 1940 ellipsoid and uses a Stereographic projection.
    namespace stereo70
    {
        /// @brief Alias for `projected_coord<T>` representing coordinates in the Stereo70 system.
        /// @tparam T The floating-point type for coordinate values.
        template <typename T>
        using coordinate = projected_coord<T>;

        /// @brief Defines the core projection parameters for the Stereo70 system.
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
            /// @brief Latitude of the projection origin (degrees).
            static constexpr T lat0_deg = T(46.0);
            /// @brief Longitude of the projection origin / Central Meridian (degrees).
            static constexpr T lon0_deg = T(25.0);
            /// @brief Scale factor at the projection origin.
            static constexpr T k0 = T(0.99975);
            /// @brief False Easting applied to the X coordinate (meters).
            static constexpr T fe = T(500000.0);
            /// @brief False Northing applied to the Y coordinate (meters).
            static constexpr T fn = T(500000.0);

            // Tolerances specific to Stereo70 calculations
            /// @brief Tolerance for checking if a point is at the projection origin (degrees).
            static constexpr T origin_tol_deg = T(1e-9);
            /// @brief Tolerance for checking proximity to the antipodal point of the projection origin.
            static constexpr T antipodal_tol = T(1e-10);
            /// @brief Tolerance for checking if a projected coordinate is at the origin (meters).
            static constexpr T origin_proj_tol = T(1e-6);
        };
    } // namespace stereo70

    /// @brief Provides static functions for performing GIS coordinate conversions.
    ///
    /// This struct contains the logic for transformations between different coordinate systems
    /// and datums, such as:
    /// - Geodetic to Geocentric (ECEF)
    /// - Geocentric (ECEF) to Geodetic
    /// - Helmert transformation (forward and reverse)
    /// - Krasovsky Geodetic to Stereo70 Projected
    /// - Stereo70 Projected to Krasovsky Geodetic
    /// - WGS84 Geodetic to Stereo70 Projected (via Helmert and Krasovsky)
    /// - Stereo70 Projected to WGS84 Geodetic (via Krasovsky and Helmert)
    ///
    /// It utilizes predefined ellipsoids (WGS84, Krasovsky 1940) and default Helmert parameters
    /// (Pulkovo 1958 to WGS84).
    ///
    /// @tparam T The floating-point type (e.g., `float`, `double`, `long double`) used for calculations.
    ///           Must satisfy `std::is_floating_point_v<T>`.
    template <typename T>
    struct conversion
    {
        /// @brief Alias for constants of the corresponding floating-point type.
        using K = constants<T>;
        /// @brief Alias for Stereo70 projection parameters.
        using stereo70_params = stereo70::projection_params<T>;

        // Ellipsoid definitions
        /// @brief Definition of the Krasovsky 1940 ellipsoid.
        static constexpr ellipsoid<T> krasovsky1940 {
            T(6378245.0), // a (meters)
            T(298.3)      // 1/f (inverse flattening)
        };

        /// @brief Definition of the WGS84 ellipsoid.
        static constexpr ellipsoid<T> wgs84 {
            T(6378137.0),    // a (meters)
            T(298.257223563) // 1/f (inverse flattening)
        };

        /// @brief Default Helmert parameters for transforming from Pulkovo 1942/58 datum (approximated via Krasovsky) to WGS84.
        /// @note These are approximate parameters often used for Romanian transformations. Specific applications might require
        /// different parameter sets. The sign convention matches the Coordinate Frame Rotation method used in
        /// `helmert_transformation_forward` and
        /// `_reverse`.
        static constexpr helmert_params<T> pulkovo58_to_wgs84 {
            .dx = T(33.4),       // meters
            .dy = T(-146.6),     // meters
            .dz = T(-76.3),      // meters
            .rx_sec = T(-0.359), // arcseconds
            .ry_sec = T(-0.053), // arcseconds
            .rz_sec = T(0.844),  // arcseconds
            .ds_ppm = T(-0.84)   // parts per million
        };

        /// @brief Pre-calculates and stores derived parameters needed for Stereo70 projection conversions.
        ///
        /// This nested struct computes values based on the Krasovsky ellipsoid and Stereo70
        /// projection parameters (`lat0_deg`, `lon0_deg`, `k0`) that are frequently reused
        /// in the projection and inverse projection calculations, improving efficiency and clarity.
        /// These parameters are constant after construction.
        ///
        /// @note This struct assumes the `krasovsky1940` ellipsoid and `stereo70::projection_params` are defined.
        struct derived_params
        {
            // Members are constant after construction
            const T lat0_rad;   ///< Latitude of origin in radians.
            const T lon0_rad;   ///< Longitude of origin (central meridian) in radians.
            const T e;          ///< First eccentricity of the Krasovsky ellipsoid.
            const T e2;         ///< First eccentricity squared of the Krasovsky ellipsoid.
            const T e4;         ///< e^4.
            const T e6;         ///< e^6.
            const T e8;         ///< e^8.
            const T sin_lat0;   ///< Sine of the latitude of origin.
            const T cos_lat0;   ///< Cosine of the latitude of origin.
            const T chi0;       ///< Conformal latitude of the origin.
            const T sin_chi0;   ///< Sine of the conformal latitude of the origin.
            const T cos_chi0;   ///< Cosine of the conformal latitude of the origin.
            const T n0;         ///< Radius of the projection sphere at the origin scaled by k0.
            const T inv_two_n0; ///< 1 / (2 * n0).
            const T c1;         ///< Coefficient C1 for inverse conformal latitude series.
            const T c2;         ///< Coefficient C2 for inverse conformal latitude series.
            const T c3;         ///< Coefficient C3 for inverse conformal latitude series.
            const T c4;         ///< Coefficient C4 for inverse conformal latitude series.

            /// @brief Constructs and calculates the derived projection parameters.
            /// @throws std::runtime_error if calculated parameters are invalid (e.g., NaN, division by zero).
            constexpr derived_params() noexcept(false):
                lat0_rad(stereo70_params::lat0_deg * K::conversions::deg_to_rad),
                lon0_rad(stereo70_params::lon0_deg * K::conversions::deg_to_rad),
                e(krasovsky1940.e),
                e2(krasovsky1940.e2),
                e4(e2 * e2),
                e6(e4 * e2),
                e8(e6 * e2),
                sin_lat0(std::sin(lat0_rad)),
                cos_lat0(std::cos(lat0_rad)),
                chi0(calculate_conformal_latitude(lat0_rad, e)), // Use static helper
                sin_chi0(std::sin(chi0)),
                cos_chi0(std::cos(chi0)),
                n0(calculate_n0()), // Use const helper method
                inv_two_n0(std::abs(n0) < K::tolerance::near_zero ? std::numeric_limits<T>::quiet_NaN() :
                                                                    K::values::half / n0), // Prevent division by zero
                c1(calculate_c1()),                                                        // Use const helper method
                c2(calculate_c2()),                                                        // Use const helper method
                c3(calculate_c3()),                                                        // Use const helper method
                c4(calculate_c4())                                                         // Use const helper method
            {
                // Validate calculated parameters after construction
                if (!is_valid())
                    throw std::runtime_error("Invalid derived projection parameters calculated (e.g., NaN or zero denominator)");
            }

        private:
            /// @brief Calculates N0, the radius of the projection sphere at the origin scaled by k0.
            /// Helper function for the constructor.
            /// @return Calculated N0 value, or NaN if inputs are invalid.
            [[nodiscard]] constexpr T calculate_n0() const noexcept
            {
                // W1^2 = 1 - e^2 * sin^2(lat0)
                const T w1_sq = K::values::one - e2 * sin_lat0 * sin_lat0;
                // Ensure W1 = sqrt(W1^2) is real and positive
                if (w1_sq <= K::values::zero)
                    return std::numeric_limits<T>::quiet_NaN(); // Invalid input latitude or eccentricity
                const T w1 = std::sqrt(w1_sq);

                // Denominator = W1 * cos(chi0)
                const T denominator = w1 * cos_chi0;
                // Ensure denominator is not too close to zero
                if (std::abs(denominator) < K::tolerance::near_zero)
                    return std::numeric_limits<T>::quiet_NaN(); // Origin near a pole or conformal latitude calculation issue

                // N0 = k0 * a * cos(lat0) / denominator
                // Use k0 from stereo70_params directly
                return stereo70_params::k0 * krasovsky1940.a * cos_lat0 / denominator;
            }

            /// @brief Calculates C1 coefficient. Helper for constructor. @return C1.
            [[nodiscard]] constexpr T calculate_c1() const noexcept
            {
                // C1 = e2/2 + 5e4/24 + e6/12 + 13e8/360
                return e2 * K::values::half + K::values::five * e4 / K::values::twenty_four + e6 / K::values::twelve +
                       K::values::thirteen * e8 / K::values::three_sixty;
            }

            /// @brief Calculates C2 coefficient. Helper for constructor. @return C2.
            [[nodiscard]] constexpr T calculate_c2() const noexcept
            {
                // C2 = 7e4/48 + 29e6/240 + 811e8/11520
                return (K::series::c2_term1_num / K::series::c2_term1_den) * e4 + (K::series::c2_term2_num / K::series::c2_term2_den) * e6 +
                       (K::series::c2_term3_num / K::series::c2_term3_den) * e8;
            }

            /// @brief Calculates C3 coefficient. Helper for constructor. @return C3.
            [[nodiscard]] constexpr T calculate_c3() const noexcept
            {
                // C3 = 7e6/120 + 81e8/1120
                return (K::series::c3_term1_num / K::series::c3_term1_den) * e6 + (K::series::c3_term2_num / K::series::c3_term2_den) * e8;
            }

            /// @brief Calculates C4 coefficient. Helper for constructor. @return C4.
            [[nodiscard]] constexpr T calculate_c4() const noexcept
            {
                // C4 = 4279e8 / 161280
                return (K::series::c4_term1_num / K::series::c4_term1_den) * e8;
            }

            /// @brief Validates that all calculated derived parameters are finite and usable.
            /// Helper function for the constructor.
            /// @return `true` if all parameters are valid, `false` otherwise.
            [[nodiscard]] constexpr bool is_valid() const noexcept
            {
                using std::isfinite;
                // Check essential parameters for finiteness and non-zero denominators where applicable
                return isfinite(lat0_rad) && isfinite(lon0_rad) && isfinite(e) && isfinite(e2) && // Assuming e, e2 are valid from ellipsoid
                       isfinite(sin_lat0) && isfinite(cos_lat0) && isfinite(chi0) && isfinite(sin_chi0) && isfinite(cos_chi0) &&
                       isfinite(n0) && std::abs(n0) > K::tolerance::near_zero &&     // Check n0 is finite and not near zero
                       isfinite(inv_two_n0) &&                                       // Check derived 1/(2*n0)
                       isfinite(c1) && isfinite(c2) && isfinite(c3) && isfinite(c4); // Check series coefficients
            }
        }; // end struct derived_params

        /// @brief Static constant instance of the derived projection parameters for Stereo70/Krasovsky.
        /// Initialized at compile time (or first use for static local).
        static constexpr derived_params proj_params {};

        // Static Helper Functions

        /// @brief Calculates the conformal latitude (χ) from the geodetic latitude (φ).
        ///
        /// Conformal latitude is used in some map projections (like Stereographic) to make them conformal.
        /// Formula: χ = 2 * atan( tan(π/4 + φ/2) * [(1 - e*sin(φ)) / (1 + e*sin(φ))]^(e/2) ) - π/2
        ///
        /// @param[in] lat_rad Geodetic latitude in radians.
        /// @param[in] e First eccentricity of the reference ellipsoid.
        /// @return Conformal latitude in radians, or +/- π/2 for poles, or NaN if calculation fails.
        /// @note Handles proximity to poles and potential calculation issues.
        [[nodiscard]] static constexpr T calculate_conformal_latitude(const T lat_rad, const T e) noexcept
        {
            using K = constants<T>; // Use constants specific to this function's scope if needed

            // Handle poles directly
            if (std::abs(std::abs(lat_rad) - K::math::half_pi) < K::tolerance::near_zero)
            {
                return std::copysign(K::math::half_pi, lat_rad);
            }

            const T sin_lat = std::sin(lat_rad);
            const T esin_lat = e * sin_lat;

            // Check for potential issue |e*sin(lat)| >= 1 (only possible if e >= 1, which ellipsoid validation should prevent)
            if (std::abs(esin_lat) >= K::values::one)
            {
                // This case implies e >= 1/|sin(lat)|. If e < 1, this shouldn't happen.
                // If it does, result is likely infinite or undefined. Return pole latitude.
                return std::copysign(K::math::half_pi, lat_rad);
            }

            // Term 1: tan(π/4 + φ/2)
            const T term1_arg = K::math::quarter_pi + lat_rad * K::values::half;
            const T term1 = std::tan(term1_arg);

            // Term 2: [(1 - e*sin(φ)) / (1 + e*sin(φ))]^(e/2)
            const T term2_num = K::values::one - esin_lat;
            const T term2_den = K::values::one + esin_lat;

            // Check for division by zero or negative base for power term
            if (std::abs(term2_den) < K::tolerance::near_zero)
                return std::numeric_limits<T>::quiet_NaN();                   // Division by zero
            if (term2_num <= K::values::zero || term2_den <= K::values::zero) // Base must be positive for real power e/2
                // This condition (term2_num <= 0 or term2_den <= 0) implies |e*sin(lat)| >= 1, handled above.
                // Re-check anyway for robustness.
                return std::numeric_limits<T>::quiet_NaN(); // Or return pole latitude?

            const T term2_base = term2_num / term2_den;
            const T term2 = std::pow(term2_base, e * K::values::half);

            // Argument for final atan: term1 * term2
            const T atan_arg = term1 * term2;

            // Check if atan argument became non-finite
            if (!std::isfinite(atan_arg))
            {
                // If argument is infinite, atan approaches +/- pi/2
                return std::copysign(K::math::half_pi, atan_arg);
            }

            // Final calculation: χ = 2 * atan(atan_arg) - π/2
            return K::values::two * std::atan(atan_arg) - K::math::half_pi;
        }

        // Core Transformation Functions

        /// @brief Converts Geodetic coordinates (lat, lon, height) to Geocentric Cartesian (X, Y, Z) coordinates.
        ///
        /// Uses the provided ellipsoid definition for the transformation.
        /// Formulae:
        ///   N = a / sqrt(1 - e^2 * sin^2(lat))
        ///   X = (N + h) * cos(lat) * cos(lon)
        ///   Y = (N + h) * cos(lat) * sin(lon)
        ///   Z = (N * (1 - e^2) + h) * sin(lat)
        ///
        /// @param[in] geo The input geodetic coordinate (latitude/longitude in degrees, height in meters).
        /// @param[in] ellip The reference ellipsoid parameters.
        /// @return The corresponding geocentric (ECEF) coordinate (X, Y, Z in meters).
        /// @throws std::invalid_argument if input `geo` coordinates are invalid (via `geo.validate()`).
        /// @throws std::runtime_error if the calculation encounters invalid states (e.g., division by zero due to N denominator).
        [[nodiscard]] static geocentric_coord<T> geodetic_to_geocentric(const geodetic_coord<T>& geo,
                                                                        const ellipsoid<T>& ellip) noexcept(false)
        {
            // Validate input coordinate first
            geo.validate(); // Throws std::invalid_argument if invalid

            // Convert lat/lon to radians
            const T lat_rad = geo.latitude * K::conversions::deg_to_rad;
            const T lon_rad = geo.longitude * K::conversions::deg_to_rad;
            const T sin_lat = std::sin(lat_rad);
            const T cos_lat = std::cos(lat_rad);

            // Calculate N, the radius of curvature in the prime vertical
            const T n_denominator_sq = K::values::one - ellip.e2 * sin_lat * sin_lat;
            // Check if denominator is positive before sqrt
            if (n_denominator_sq <= K::values::zero) // Should not happen for valid latitudes and e < 1
                throw std::runtime_error(
                    "Invalid geodetic coordinate or ellipsoid parameter leading to non-positive N denominator squared: " +
                    std::to_string(n_denominator_sq));

            const T n_denominator = std::sqrt(n_denominator_sq);
            if (std::abs(n_denominator) < K::tolerance::near_zero) // Should be covered by above check, but for safety
                throw std::runtime_error("Invalid geodetic coordinate or ellipsoid parameter leading to near-zero N denominator: " +
                                         std::to_string(n_denominator));

            const T N = ellip.a / n_denominator;
            const T h = geo.height; // Ellipsoidal height

            // Calculate X, Y, Z
            const T N_plus_h = N + h;
            const T xy_factor = N_plus_h * cos_lat; // Common factor for X and Y

            const T x = xy_factor * std::cos(lon_rad);
            const T y = xy_factor * std::sin(lon_rad);
            // Z = (N * (1 - e^2) + h) * sin(lat) = (N * b^2/a^2 + h) * sin(lat)
            const T z = (N * (K::values::one - ellip.e2) + h) * sin_lat;

            // Construct and return result
            geocentric_coord<T> result {.x = x, .y = y, .z = z};
            // Post-validation (optional, primarily checks for NaN/inf results)
            result.validate(); // Throws std::invalid_argument if non-finite
            return result;
        }

        /// @brief Converts Geocentric Cartesian (X, Y, Z) coordinates to Geodetic coordinates (lat, lon, height).
        ///
        /// Uses an iterative algorithm (e.g., Bowring's or similar) to find the latitude and height.
        /// Longitude is calculated directly from atan2(Y, X).
        ///
        /// @param[in] ecef The input geocentric (ECEF) coordinate (X, Y, Z in meters).
        /// @param[in] ellip The reference ellipsoid parameters.
        /// @param[in] tolerance The convergence tolerance for the latitude iteration (radians). Defaults to
        /// `K::tolerance::default_geodetic`.
        /// @param[in] max_iterations The maximum number of iterations allowed. Defaults to `K::tolerance::max_iterations`.
        /// @return The corresponding geodetic coordinate (latitude/longitude in degrees, height in meters). The coordinates are
        /// normalized.
        /// @throws std::invalid_argument if input `ecef` coordinates are invalid (via `ecef.validate()`).
        /// @throws std::runtime_error if the iteration fails to converge or encounters numerical instability.
        /// @note Handles special cases like points on the Z-axis (poles) and near the Earth's center.
        [[nodiscard]] static geodetic_coord<T> geocentric_to_geodetic(
            const geocentric_coord<T>& ecef, const ellipsoid<T>& ellip, const T tolerance = K::tolerance::default_geodetic,
            const int max_iterations = K::tolerance::max_iterations) noexcept(false)
        {
            // Validate input coordinate first
            ecef.validate(); // Throws std::invalid_argument if invalid

            const T x = ecef.x;
            const T y = ecef.y;
            const T z = ecef.z;

            // Calculate distance from Z-axis
            const T p = std::hypot(x, y);

            geodetic_coord<T> result; // Result coordinate (will be modified)

            // Handle singularity at poles (p is near zero)
            // Define a threshold based on ellipsoid size and machine precision
            const T pole_threshold = K::tolerance::near_zero * ellip.a; // Scale near_zero by semi-major axis
            if (p < pole_threshold)
            {
                result.longitude = K::values::zero; // Longitude is undefined at poles, conventionally set to 0.

                // Check if also near the center
                if (std::abs(z) < pole_threshold)
                {
                    // Point is very close to the origin
                    result.latitude = K::values::zero;
                    result.height =
                        -ellip.b; // Height relative to ellipsoid surface at origin (closest point is pole b) - or use -a? Using -b.
                    // Alternatively: throw std::runtime_error("Point is near the center of the ellipsoid, geodetic coordinates are
                    // ill-defined.");
                }
                else
                {
                    // Point is on or very near the Z-axis (pole)
                    result.latitude = std::copysign(K::values::ninety, z); // +/- 90 degrees
                    // Height is distance along Z from polar surface point
                    result.height = std::abs(z) - ellip.b; // Height above the semi-minor axis endpoint
                }
                // No iteration needed, result is determined. Already normalized implicitly.
                result.validate(); // Final check
                return result;
            }

            // General Case (Iterative Solution)

            // Calculate longitude directly
            result.longitude = std::atan2(y, x) * K::conversions::rad_to_deg;

            // Initial guess for latitude (based on spherical approximation or reduced latitude)
            // Using atan2(z, p * (1 - e^2)) is a common starting point (related to reduced latitude)
            T lat_rad = std::atan2(z, p * (K::values::one - ellip.e2));
            T lat_rad_prev = lat_rad + tolerance * K::values::two; // Ensure first loop runs

            T N = ellip.a;         // Radius of curvature in prime vertical
            T h = K::values::zero; // Ellipsoidal height

            int i = 0;
            for (; i < max_iterations; ++i)
            {
                const T sin_lat = std::sin(lat_rad);
                const T cos_lat = std::cos(lat_rad);

                // Calculate N for the current latitude estimate
                const T N_denominator_sq = K::values::one - ellip.e2 * sin_lat * sin_lat;
                if (N_denominator_sq <= K::tolerance::near_zero) // Use tolerance here, check could be <= 0 strictly
                    throw std::runtime_error("Numerical instability: N denominator near or below zero during iteration.");

                N = ellip.a / std::sqrt(N_denominator_sq);

                // Calculate height h
                // Avoid division by cos_lat near poles
                if (std::abs(cos_lat) < K::tolerance::near_zero)
                {
                    // Near pole: h = |z| / |sin(lat)| - N * (1 - e^2)
                    // Since sin(lat) will be close to +/- 1, this is safer.
                    h = std::abs(z) / std::abs(sin_lat) - N * (K::values::one - ellip.e2);
                    // Latitude is already effectively +/- pi/2, refinement below won't help much
                    // We can perhaps break here if lat_rad is already very close to +/- pi/2?
                }
                else
                {
                    h = p / cos_lat - N;
                }

                // Update latitude estimate
                // Avoid singularity if point is near the center (N+h is small)
                const T N_plus_h = N + h;
                const T center_threshold = K::tolerance::near_zero * ellip.a; // Threshold for being near center

                T lat_rad_next;
                if (std::abs(N_plus_h) < center_threshold)
                {
                    // Point is very close to the center, latitude is ill-defined.
                    // Revert to a simple approximation or stick to pole behavior.
                    // If p > 0 (checked earlier), maybe use atan2(z, p)?
                    // Or signal failure/stick to previous value? Sticking to pole approx.
                    lat_rad_next = (std::abs(z) < center_threshold) ? K::values::zero : // Assume equatorial plane if z is also small
                                                                      std::copysign(K::math::half_pi, z);              // Stick to pole otherwise
                    // Convergence check might fail here, potentially leading to max iterations.
                    // Consider breaking? For now, let it continue.
                }
                else
                {
                    // Standard iteration formula: lat = atan( Z / (P * (1 - e^2 * N / (N + h))) )
                    const T lat_denom = p * (K::values::one - ellip.e2 * N / N_plus_h);
                    if (std::abs(lat_denom) < center_threshold) // Check denominator too (though N+h check helps)
                    {
                        // Denominator near zero implies point is near Z-axis (pole)
                        lat_rad_next = std::copysign(K::math::half_pi, z);
                    }
                    else
                    {
                        lat_rad_next = std::atan2(z, lat_denom);
                    }
                }

                lat_rad_prev = lat_rad;
                lat_rad = lat_rad_next;

                // Check for convergence
                if (std::abs(lat_rad - lat_rad_prev) < tolerance)
                    break; // Converged
            }

            // Check if max iterations were reached without convergence
            if (i >= max_iterations)
                throw std::runtime_error("Geocentric to Geodetic conversion failed to converge within maximum iterations.");

            // Store final results
            result.latitude = lat_rad * K::conversions::rad_to_deg;
            result.height = h;

            // Normalize the final coordinates (especially longitude)
            result.normalize();
            // Final validation of calculated results
            result.validate(); // Throws if result is invalid (e.g., NaN)
            return result;
        }

        /// @brief Converts Krasovsky 1940 Geodetic coordinates to Stereo70 projected coordinates.
        ///
        /// Implements the forward Stereographic projection formula as defined for Stereo70.
        ///
        /// @param[in] geo The input geodetic coordinate (lat/lon in degrees, height in meters) on the Krasovsky 1940 ellipsoid.
        /// @return The corresponding Stereo70 projected coordinate (X=Easting, Y=Northing, Z=Height) in meters.
        /// @throws std::invalid_argument if input `geo` coordinates are invalid (via `geo.validate()`).
        /// @throws std::runtime_error if the coordinate is too close to the antipodal point of the projection origin or if calculation
        /// fails.
        /// @sa stereo70_to_krasovsky_geodetic()
        [[nodiscard]] static projected_coord<T> krasovsky_geodetic_to_stereo70(const geodetic_coord<T>& geo) noexcept(false)
        {
            // Validate input coordinate first
            geo.validate(); // Throws std::invalid_argument if invalid

            // Handle the origin case directly for precision and speed
            // Use tolerances defined in stereo70_params
            if (std::abs(geo.latitude - stereo70_params::lat0_deg) < stereo70_params::origin_tol_deg &&
                std::abs(geo.longitude - stereo70_params::lon0_deg) < stereo70_params::origin_tol_deg)
            {
                // At the origin, result is False Easting, False Northing
                return {.x = stereo70_params::fe, .y = stereo70_params::fn, .z = geo.height};
            }

            // Convert lat/lon to radians
            const T lat_rad = geo.latitude * K::conversions::deg_to_rad;
            const T lon_rad = geo.longitude * K::conversions::deg_to_rad;

            // Calculate conformal latitude (chi) for the input point
            const T chi = calculate_conformal_latitude(lat_rad, proj_params.e);
            if (!std::isfinite(chi)) // Check if conformal latitude calculation was successful
                throw std::runtime_error("Failed to calculate valid conformal latitude for input geodetic coordinates.");

            const T sin_chi = std::sin(chi);
            const T cos_chi = std::cos(chi);

            // Calculate longitude difference (lambda) relative to the central meridian, normalized to (-pi, pi]
            T lambda_diff = lon_rad - proj_params.lon0_rad;
            // Normalize lambda_diff to (-pi, pi] range
            while (lambda_diff > K::math::pi)
                lambda_diff -= K::math::two_pi;
            while (lambda_diff <= -K::math::pi)
                lambda_diff += K::math::two_pi;

            const T sin_lambda = std::sin(lambda_diff);
            const T cos_lambda = std::cos(lambda_diff);

            // Calculate the denominator D = 1 + sin(chi0)*sin(chi) + cos(chi0)*cos(chi)*cos(lambda)
            const T denominator = K::values::one + proj_params.sin_chi0 * sin_chi + proj_params.cos_chi0 * cos_chi * cos_lambda;

            // Check for projection singularity (antipodal point)
            // Use tolerance defined in stereo70_params
            if (std::abs(denominator) < stereo70_params::antipodal_tol)
                throw std::runtime_error("Input coordinate is too close to the antipodal point of the projection origin.");

            // Calculate the projection scale factor k' = 2 * N0 / D
            const T k_prime = K::values::two * proj_params.n0 / denominator;
            if (!std::isfinite(k_prime)) // Check if k_prime is valid
                throw std::runtime_error("Projection scale factor (k_prime) calculation resulted in non-finite value.");

            // Calculate Easting (X) and Northing (Y)
            // X = FE + k' * cos(chi) * sin(lambda)
            // Y = FN + k' * (cos(chi0)*sin(chi) - sin(chi0)*cos(chi)*cos(lambda))
            // Use fe, fn from stereo70_params
            const T x_proj = stereo70_params::fe + k_prime * cos_chi * sin_lambda;
            const T y_proj = stereo70_params::fn + k_prime * (proj_params.cos_chi0 * sin_chi - proj_params.sin_chi0 * cos_chi * cos_lambda);

            projected_coord<T> result {.x = x_proj, .y = y_proj, .z = geo.height}; // Z coordinate is typically the ellipsoidal height

            // Validate the result (check for NaN/inf)
            // Use the validate method of xyz_coord
            if (!std::isfinite(result.x) || !std::isfinite(result.y)) // Basic check, validate() is more robust
                throw std::runtime_error("Non-finite projection result (X or Y).");
            result.validate(); // Throws std::invalid_argument if non-finite

            return result;
        }

        /// @brief Converts WGS84 Geodetic coordinates to Stereo70 projected coordinates.
        ///
        /// This is a composite transformation:
        /// 1. WGS84 Geodetic -> WGS84 ECEF (`geodetic_to_geocentric`)
        /// 2. WGS84 ECEF -> Pulkovo ECEF (using reverse Helmert transformation) (`helmert_transformation_reverse`)
        /// 3. Pulkovo ECEF -> Krasovsky Geodetic (`geocentric_to_geodetic`)
        /// 4. Krasovsky Geodetic -> Stereo70 Projected (`krasovsky_geodetic_to_stereo70`)
        ///
        /// @param[in] wgs84_coord The input WGS84 geodetic coordinate.
        /// @param[in] params The Helmert transformation parameters to go from WGS84 to the system
        ///                   underlying Krasovsky/Stereo70 (typically Pulkovo 1942/58). Defaults to
        ///                   `conversion<T>::pulkovo58_to_wgs84` (which needs reversal).
        /// @return The corresponding Stereo70 projected coordinate.
        /// @throws std::invalid_argument if input `wgs84_coord` or `params` are invalid.
        /// @throws std::runtime_error if any intermediate conversion step fails.
        /// @sa stereo70_to_wgs84()
        [[nodiscard]] static stereo70::coordinate<T> wgs84_to_stereo70(
            const wgs84::coordinate<T>& wgs84_coord, const helmert_params<T>& params = conversion<T>::pulkovo58_to_wgs84) noexcept(false)
        {
            // Validate inputs
            params.validate();      // Throws if Helmert params invalid
            wgs84_coord.validate(); // Throws if WGS84 coords invalid

            // Step 1: WGS84 Geodetic -> WGS84 ECEF
            const auto wgs84_ecef = geodetic_to_geocentric(wgs84_coord, wgs84); // Can throw

            // Step 2: WGS84 ECEF -> Pulkovo ECEF (using REVERSE Helmert, as params are defined Pulkovo->WGS84)
            const auto pulkovo_ecef = helmert_transformation_reverse(wgs84_ecef, params); // Can throw

            // Step 3: Pulkovo ECEF -> Krasovsky Geodetic
            // Note: Assumes Pulkovo ECEF is compatible with Krasovsky ellipsoid for geodetic conversion
            const auto krasovsky_geo = geocentric_to_geodetic(pulkovo_ecef, krasovsky1940); // Can throw

            // Step 4: Krasovsky Geodetic -> Stereo70 Projected
            return krasovsky_geodetic_to_stereo70(krasovsky_geo); // Can throw
        }

        /// @brief Converts Stereo70 projected coordinates to WGS84 Geodetic coordinates.
        ///
        /// This is the inverse composite transformation:
        /// 1. Stereo70 Projected -> Krasovsky Geodetic (`stereo70_to_krasovsky_geodetic`)
        /// 2. Krasovsky Geodetic -> Pulkovo ECEF (`geodetic_to_geocentric`)
        /// 3. Pulkovo ECEF -> WGS84 ECEF (using forward Helmert transformation) (`helmert_transformation_forward`)
        /// 4. WGS84 ECEF -> WGS84 Geodetic (`geocentric_to_geodetic`)
        ///
        /// @param[in] stereo_coord The input Stereo70 projected coordinate.
        /// @param[in] params The Helmert transformation parameters to go from the system
        ///                   underlying Krasovsky/Stereo70 (typically Pulkovo 1942/58) to WGS84.
        ///                   Defaults to `conversion<T>::pulkovo58_to_wgs84`.
        /// @return The corresponding WGS84 geodetic coordinate. The coordinates are normalized.
        /// @throws std::invalid_argument if input `stereo_coord` or `params` are invalid.
        /// @throws std::runtime_error if any intermediate conversion step fails.
        /// @sa wgs84_to_stereo70()
        [[nodiscard]] static wgs84::coordinate<T> stereo70_to_wgs84(
            const stereo70::coordinate<T>& stereo_coord,
            const helmert_params<T>& params = conversion<T>::pulkovo58_to_wgs84) noexcept(false)
        {
            // Validate inputs
            params.validate();       // Throws if Helmert params invalid
            stereo_coord.validate(); // Throws if Stereo70 coords invalid

            // Step 1: Stereo70 Projected -> Krasovsky Geodetic
            const auto krasovsky_geo = stereo70_to_krasovsky_geodetic(stereo_coord); // Can throw

            // Step 2: Krasovsky Geodetic -> Pulkovo ECEF
            // Note: Assumes Pulkovo ECEF is compatible with Krasovsky ellipsoid for geodetic conversion
            const auto pulkovo_ecef = geodetic_to_geocentric(krasovsky_geo, krasovsky1940); // Can throw

            // Step 3: Pulkovo ECEF -> WGS84 ECEF (using FORWARD Helmert, as params are defined Pulkovo->WGS84)
            const auto wgs84_ecef = helmert_transformation_forward(pulkovo_ecef, params); // No throw expected here, just calculation

            // Step 4: WGS84 ECEF -> WGS84 Geodetic
            return geocentric_to_geodetic(wgs84_ecef, wgs84); // Can throw
        }

        /// @brief Converts Stereo70 projected coordinates (X, Y, Z=Height) to Krasovsky 1940 Geodetic coordinates.
        ///
        /// Implements the inverse Stereographic projection formula as defined for Stereo70, using
        /// a series expansion to convert conformal latitude back to geodetic latitude.
        ///
        /// @param[in] proj The input Stereo70 projected coordinate (X=Easting, Y=Northing, Z=Height) in meters.
        /// @return The corresponding geodetic coordinate (lat/lon in degrees, height in meters) on the Krasovsky 1940 ellipsoid. The
        /// coordinates are normalized.
        /// @throws std::invalid_argument if input `proj` coordinates are invalid (via `proj.validate()`).
        /// @throws std::runtime_error if the calculation encounters numerical instability or invalid intermediate values.
        /// @sa krasovsky_geodetic_to_stereo70()
        [[nodiscard]] static geodetic_coord<T> stereo70_to_krasovsky_geodetic(const projected_coord<T>& proj) noexcept(false)
        {
            // Validate input projected coordinate
            proj.validate(); // Throws std::invalid_argument if invalid

            // Calculate coordinates relative to false origin
            // Use fe, fn from stereo70_params
            const T x_rel = proj.x - stereo70_params::fe;
            const T y_rel = proj.y - stereo70_params::fn;

            // Calculate distance rho from the false origin
            const T rho = std::hypot(x_rel, y_rel);

            // Handle the origin case
            // Use tolerance defined in stereo70_params
            if (rho < stereo70_params::origin_proj_tol)
            {
                // At the origin, return the projection origin coordinates
                // Use lat0_deg, lon0_deg from stereo70_params
                return {.latitude = stereo70_params::lat0_deg, .longitude = stereo70_params::lon0_deg, .height = proj.z};
            }

            // Calculate intermediate angle 'c'
            // tan(c/2) = rho / (2 * N0)
            // Use pre-calculated inv_two_n0 = 1 / (2 * N0)
            if (!std::isfinite(proj_params.inv_two_n0)) // Check if precalculated value is valid
                throw std::runtime_error("Invalid derived parameter inv_two_n0 for inverse projection.");

            const T tan_c_half = rho * proj_params.inv_two_n0;
            if (!std::isfinite(tan_c_half)) // Check calculation result
                throw std::runtime_error("Invalid intermediate value (tan_c_half) in inverse projection calculation.");

            const T c = K::values::two * std::atan(tan_c_half);
            const T sin_c = std::sin(c);
            const T cos_c = std::cos(c);

            // Calculate conformal latitude (chi)
            // sin(chi) = cos(c)*sin(chi0) + (y' * sin(c) * cos(chi0) / rho)
            const T chi_arg_num = cos_c * proj_params.sin_chi0 + (y_rel * sin_c * proj_params.cos_chi0 / rho);
            // Clamp the argument to [-1, 1] before asin to handle potential floating point inaccuracies
            const T chi_arg = std::clamp(chi_arg_num, -K::values::one, K::values::one);
            const T chi = std::asin(chi_arg);
            if (!std::isfinite(chi)) // Check result of asin
                throw std::runtime_error("Failed to calculate valid conformal latitude (chi) in inverse projection.");

            // Calculate longitude difference (lambda)
            // tan(lambda) = (x' * sin(c)) / (rho * cos(chi0) * cos(c) - y' * sin(chi0) * sin(c))
            const T lambda_diff_num = x_rel * sin_c;
            const T lambda_diff_den = rho * proj_params.cos_chi0 * cos_c - y_rel * proj_params.sin_chi0 * sin_c;
            // Use atan2 for correct quadrant
            const T lambda_diff = std::atan2(lambda_diff_num, lambda_diff_den);
            if (!std::isfinite(lambda_diff)) // Check result of atan2
                throw std::runtime_error("Failed to calculate valid longitude difference (lambda_diff) in inverse projection.");

            // Calculate geodetic latitude (lat_rad) from conformal latitude (chi) using series expansion
            // lat = chi + C1*sin(2chi) + C2*sin(4chi) + C3*sin(6chi) + C4*sin(8chi)
            // Use pre-calculated coefficients c1, c2, c3, c4
            const T lat_rad = calculate_geodetic_latitude_series(chi, proj_params.c1, proj_params.c2, proj_params.c3, proj_params.c4);
            if (!std::isfinite(lat_rad)) // Check result of series calculation
                throw std::runtime_error("Failed to calculate valid geodetic latitude (lat_rad) from series expansion.");

            // Construct the result geodetic coordinate
            geodetic_coord<T> result {
                .latitude = lat_rad * K::conversions::rad_to_deg, // Convert back to degrees
                .longitude =
                    (proj_params.lon0_rad + lambda_diff) * K::conversions::rad_to_deg, // Add central meridian back, convert to degrees
                .height = proj.z                                                       // Height is preserved
            };

            // Normalize and validate the final result
            result.normalize(); // Ensure lat/lon are in standard ranges
            result.validate();  // Throws if result is invalid (e.g., NaN, lat out of range)
            return result;
        }

    private:
        // Private Static Helper Functions

        /// @brief Calculates geodetic latitude from conformal latitude using a series expansion.
        ///
        /// Formula: φ = χ + C1*sin(2χ) + C2*sin(4χ) + C3*sin(6χ) + C4*sin(8χ)
        /// Used in the inverse Stereographic projection for Stereo70.
        ///
        /// @param[in] chi Conformal latitude in radians.
        /// @param[in] c1 Series coefficient C1.
        /// @param[in] c2 Series coefficient C2.
        /// @param[in] c3 Series coefficient C3.
        /// @param[in] c4 Series coefficient C4.
        /// @return Geodetic latitude in radians.
        /// @note This is a private helper, assuming valid inputs.
        [[nodiscard]] static constexpr T calculate_geodetic_latitude_series(const T chi, const T c1, const T c2, const T c3,
                                                                            const T c4) noexcept
        {
            // Pre-calculate multiples of chi for sin arguments
            const T chi2 = K::values::two * chi;
            const T chi4 = K::values::four * chi;
            const T chi6 = K::values::six * chi;
            const T chi8 = K::values::eight * chi;

            // Calculate the series sum
            return chi + c1 * std::sin(chi2) + c2 * std::sin(chi4) + c3 * std::sin(chi6) + c4 * std::sin(chi8);
        }

        /// @brief Applies a forward 7-parameter Helmert transformation (Coordinate Frame method).
        ///
        /// Transforms coordinates FROM the source system defined by `params` TO the target system (e.g., Pulkovo ECEF -> WGS84 ECEF).
        /// Formula (small angle approximation):
        ///   X_target = dx + (1+ds) * (X_source - rz*Y_source + ry*Z_source)
        ///   Y_target = dy + (1+ds) * (rz*X_source + Y_source - rx*Z_source)
        ///   Z_target = dz + (1+ds) * (-ry*X_source + rx*Y_source + Z_source)
        /// where rx, ry, rz are rotation angles in radians.
        ///
        /// @param[in] source The source geocentric coordinate (X, Y, Z).
        /// @param[in] params The Helmert parameters defining the transformation FROM source TO target.
        /// @return The transformed geocentric coordinate in the target system.
        /// @note Assumes `params` have been validated. Does not throw exceptions itself but relies on finite inputs.
        /// @sa helmert_transformation_reverse()
        [[nodiscard]] static geocentric_coord<T> helmert_transformation_forward(const geocentric_coord<T>& source,
                                                                                const helmert_params<T>& params) noexcept
        {
            // Get parameters in needed units (radians for rotations, scale factor)
            const T rx = params.rx_rad();          // Rotation around X in radians
            const T ry = params.ry_rad();          // Rotation around Y in radians
            const T rz = params.rz_rad();          // Rotation around Z in radians
            const T scale = params.scale_factor(); // (1 + ds)

            // Apply rotation (using small angle approximation for Coordinate Frame rotation)
            const T x_rot = source.x - rz * source.y + ry * source.z;
            const T y_rot = rz * source.x + source.y - rx * source.z;
            const T z_rot = -ry * source.x + rx * source.y + source.z;

            // Apply scale and translation
            const T x_target = params.dx + scale * x_rot;
            const T y_target = params.dy + scale * y_rot;
            const T z_target = params.dz + scale * z_rot;

            return {.x = x_target, .y = y_target, .z = z_target};
        }

        /// @brief Applies a reverse 7-parameter Helmert transformation (Coordinate Frame method).
        ///
        /// Transforms coordinates FROM the target system TO the source system defined by `params` (e.g., WGS84 ECEF -> Pulkovo ECEF).
        /// Inverts the forward transformation.
        /// Formula (derived from forward, small angle approximation):
        ///   Let X' = (X_target - dx) / (1+ds), Y' = (Y_target - dy) / (1+ds), Z' = (Z_target - dz) / (1+ds)
        ///   X_source = X' + rz*Y' - ry*Z'
        ///   Y_source = -rz*X' + Y' + rx*Z'
        ///   Z_source = ry*X' - rx*Y' + Z'
        /// where rx, ry, rz are rotation angles in radians.
        ///
        /// @param[in] target The target geocentric coordinate (X, Y, Z) - the coordinate we are transforming FROM.
        /// @param[in] params The Helmert parameters defining the transformation FROM source TO target.
        /// @return The transformed geocentric coordinate in the source system.
        /// @throws std::runtime_error if the scale factor is too close to zero, preventing inversion.
        /// @note Assumes `params` have been validated for finiteness, but checks scale factor for inversion.
        /// @sa helmert_transformation_forward()
        [[nodiscard]] static geocentric_coord<T> helmert_transformation_reverse(const geocentric_coord<T>& target,
                                                                                const helmert_params<T>& params) noexcept(false)
        {
            // Get parameters in needed units
            const T rx = params.rx_rad();          // Rotation around X in radians
            const T ry = params.ry_rad();          // Rotation around Y in radians
            const T rz = params.rz_rad();          // Rotation around Z in radians
            const T scale = params.scale_factor(); // (1 + ds)

            // Check if scale factor is valid for division
            if (std::abs(scale) < K::tolerance::near_zero)
                throw std::runtime_error("Invalid scale factor (near zero: " + std::to_string(scale) +
                                         ") in reverse Helmert transformation, cannot invert.");

            const T inv_scale = K::values::one / scale;

            // Reverse translation and scale
            const T x_temp = target.x - params.dx;
            const T y_temp = target.y - params.dy;
            const T z_temp = target.z - params.dz;

            const T x_scaled = x_temp * inv_scale;
            const T y_scaled = y_temp * inv_scale;
            const T z_scaled = z_temp * inv_scale;

            // Apply inverse rotation (transpose of rotation matrix for small angles)
            const T x_source = x_scaled + rz * y_scaled - ry * z_scaled;
            const T y_source = -rz * x_scaled + y_scaled + rx * z_scaled;
            const T z_source = ry * x_scaled - rx * y_scaled + z_scaled;

            return {.x = x_source, .y = y_source, .z = z_source};
        }
    }; // end struct conversion

} // namespace kmx::gis
