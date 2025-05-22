/// Copyright (c) 2025 - present KMX Systems. All rights reserved.
/// @file gis.hpp
/// @brief Defines generic ellipsoid and Helmert parameters for GIS operations.
///
/// @note This file assumes a C++20 compliant compiler.
/// @warning Does not use PCH if PCH is not defined.
#pragma once
#ifndef PCH
    #include <kmx/gis/constants.hpp>
    #include <stdexcept> ///< For std::invalid_argument, std::runtime_error
    #include <string>    ///< For std::to_string in error messages
#endif

/// @namespace kmx::gis
/// @brief Namespace for Geographic Information System (GIS) related functionalities.
///
/// This namespace contains structures, constants, and algorithms commonly used
/// in geodetic and cartographic computations.
namespace kmx::gis
{
    /// @brief Represents a reference ellipsoid model.
    ///
    /// Stores the defining parameters of an ellipsoid (semi-major axis and inverse flattening)
    /// and calculates derived parameters like flattening, first and second eccentricity squared,
    /// first eccentricity, and semi-minor axis. These parameters are conceptually constant
    /// after the ellipsoid objectusing std::abs; is constructed.
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
        const T dx {};             ///< Translation offset along the X-axis (meters).
        const T dy {};             ///< Translation offset along the Y-axis (meters).
        const T dz {};             ///< Translation offset along the Z-axis (meters).
        const T rx_sec {};         ///< Rotation angle around the X-axis (arcseconds).
        const T ry_sec {};         ///< Rotation angle around the Y-axis (arcseconds).
        const T rz_sec {};         ///< Rotation angle around the Z-axis (arcseconds).
        const T ds_ppm {};         ///< Scale factor difference (Scale = 1 + ds_ppm * 1e-6) in parts per million (ppm).
        const char* const name {}; ///< Optional name.

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
            const T scale = scale_factor();
            if (!std::isfinite(scale))
                throw std::invalid_argument("Helmert scale factor calculates to non-finite value.");

            if (std::abs(scale) < K::tolerance::near_zero)
                throw std::invalid_argument("Helmert scale factor is dangerously close to zero.");
        }
    };
} // namespace gis
