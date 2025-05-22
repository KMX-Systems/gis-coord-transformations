/// Copyright (c) 2025 - present KMX Systems. All rights reserved.
/// @file gis/coordiante/conversion.hpp
/// @brief Defines conversion functions for GIS coordinates.
///
/// @note This file assumes a C++20 compliant compiler.
/// @warning Does not use PCH if PCH is not defined.
#pragma once
#ifndef PCH
    #include <kmx/gis/coordinate/stereo70.hpp>
    #include <kmx/gis/ellipsoid.hpp>
#endif

/// @namespace kmx::gis::stereo70
/// @brief Contains types and parameters specific to the Romanian Stereo70 projection system (EPSG:31700).
///
/// Stereo70 is based on the Krasovsky 1940 ellipsoid and uses an Oblique Stereographic projection.
/// The coordinate system definition (EPSG:31700) uses X for Northing and Y for Easting.
namespace kmx::gis::coordinate
{
    /// @brief Provides static functions for performing GIS coordinate conversions.
    ///
    /// This struct contains the logic for transformations between different coordinate systems
    /// and datums, such as:
    /// - Geodetic to Geocentric (ECEF)
    /// - Geocentric (ECEF) to Geodetic
    /// - Helmert transformation (forward and reverse)
    /// - Krasovsky Geodetic to Stereo70 Projected (EPSG:31700 convention: X=Northing, Y=Easting)
    /// - Stereo70 Projected (EPSG:31700 convention) to Krasovsky Geodetic
    /// - WGS84 Geodetic to Stereo70 Projected (via Helmert and Krasovsky)
    /// - Stereo70 Projected to WGS84 Geodetic (via Krasovsky and Helmert)
    ///
    /// It utilizes predefined ellipsoids (WGS84, Krasovsky 1940) and default Helmert parameters
    /// (Pulkovo 1942/58 to WGS84).
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
        /// @brief Definition of the Krasovsky 1940 ellipsoid (EPSG:7024).
        static constexpr ellipsoid<T> krasovsky1940 {
            T(6378245.0), // a (meters)
            T(298.3)      // 1/f (inverse flattening)
        };

        /// @brief Definition of the WGS84 ellipsoid (EPSG:7030).
        static constexpr ellipsoid<T> wgs84 {
            T(6378137.0),    // a (meters)
            T(298.257223563) // 1/f (inverse flattening)
        };

        // Helmert Parameter Sets
        struct transformation
        {
            /// @brief Helmert parameters based on EPSG Transformation Code 1241.
            /// @name Pulkovo 1942 to WGS 84 (3)
            /// @epsg_code 1241
            /// @source_crs EPSG:4284 - Pulkovo 1942
            /// @target_crs EPSG:4326 - WGS 84
            /// @accuracy ~3 meters.
            /// @area Albania; Bulgaria; Czech Republic; Germany - East Germany; Hungary; Poland; Romania; Slovakia; Former Soviet Union -
            /// onshore west of 87°E.
            /// @method Coordinate Frame rotation (EPSG:9607). Reverses CFr Bursa-Wolf (geocentric). Sign convention: P1 = P2 + T + D*R*P2.
            /// @remarks Valid transformation. Uses parameters defining the transformation FROM Pulkovo 1942 TO WGS 84. Note zero rotations
            /// for Rx, Ry.
            /// @reference https://epsg.org/transformation_1241/Pulkovo-1942-to-WGS-84-3.html
            static constexpr helmert_params<T> epsg1241 {.dx = T(24.0),      // meters
                                                         .dy = T(-124.0),    // meters
                                                         .dz = T(-79.0),     // meters
                                                         .rx_sec = T(0.0),   // arcseconds
                                                         .ry_sec = T(0.0),   // arcseconds
                                                         .rz_sec = T(-0.15), // arcseconds
                                                         .ds_ppm = T(0.09),  // parts per million
                                                         .name = "EPSG:1241 Pulkovo 1942 to WGS 84 (3)"};

            /// @brief Helmert parameters based on EPSG Transformation Code 1838.
            /// @name Dealul Piscului 1970 to WGS 84 (1)
            /// @epsg_code 1838
            /// @source_crs EPSG:4317 - Dealul Piscului 1970
            /// @target_crs EPSG:4326 - WGS 84
            /// @accuracy ~2 meters.
            /// @area Romania - onshore.
            /// @method Coordinate Frame rotation (EPSG:9607). Reverses CFr Bursa-Wolf (geocentric). Sign convention: P1 = P2 + T + D*R*P2.
            /// @remarks Valid transformation. Uses parameters defining the transformation FROM Dealul Piscului 1970 TO WGS 84. Note zero
            /// rotations (rx, ry, rz). Often used for Stereo 70 (EPSG:31700) workflows, which is based on this CRS.
            /// @reference https://epsg.org/transformation_1838/Dealul-Piscului-1970-to-WGS-84-1.html
            static constexpr helmert_params<T> epsg1838 {.dx = T(33.51),     // meters
                                                         .dy = T(-145.13),   // meters
                                                         .dz = T(-75.91),    // meters
                                                         .rx_sec = T(0.0),   // arcseconds
                                                         .ry_sec = T(0.0),   // arcseconds
                                                         .rz_sec = T(0.0),   // arcseconds
                                                         .ds_ppm = T(-0.82), // parts per million
                                                         .name = "EPSG:1838 Dealul Piscului 1970 to WGS 84 (1)"};

            /// @brief **[OBSOLETE]** Helmert parameters based on **DEPRECATED** EPSG Transformation Code 1839.
            /// @name Dealul Piscului 1970 to WGS 84 (2)
            /// @epsg_code 1839
            /// @source_crs EPSG:4317 - Dealul Piscului 1970
            /// @target_crs EPSG:4326 - WGS 84
            /// @accuracy ~5 meters.
            /// @area Romania - onshore.
            /// @method Geocentric translations (EPSG:9603). Reverses Geoc B-W GeoTr (geocentric). Sign convention: P1 = P2 + T.
            /// @remarks **DEPRECATED by EPSG.** 3-parameter transformation (translations only). Use EPSG:1838 for better accuracy if
            /// rotations/scale are significant. This variable was previously mislabeled as epsg31700 (which is the Stereo 70 *projected
            /// CRS* code).
            /// @reference https://epsg.org/transformation_1839/Dealul-Piscului-1970-to-WGS-84-2.html (May require login or show as
            /// deprecated)
            static constexpr helmert_params<T> epsg1839_obsolete {.dx = T(28.0),    // meters
                                                                  .dy = T(-121.0),  // meters
                                                                  .dz = T(-77.0),   // meters
                                                                  .rx_sec = T(0.0), // arcseconds (Implicitly zero in method 9603)
                                                                  .ry_sec = T(0.0), // arcseconds (Implicitly zero in method 9603)
                                                                  .rz_sec = T(0.0), // arcseconds (Implicitly zero in method 9603)
                                                                  .ds_ppm = T(0.0), // parts per million (Implicitly zero in method 9603)
                                                                  .name = "EPSG:1839 [DEPRECATED] Dealul Piscului 1970 to WGS 84 (2)"};

            /// @brief Commonly cited Helmert parameters for Pulkovo 1942(58) / Stereo70 to WGS 84 / ETRS89 transformation in Romania.
            /// @name Pulkovo 1942(58) to WGS 84 (Romania approx.)
            /// @epsg_code N/A - Not a standard EPSG transformation code.
            /// @source_crs Likely intended for EPSG:4179 - Pulkovo 1942(58) (used by Stereo 70)
            /// @target_crs Likely intended for EPSG:4326 - WGS 84 or EPSG:4258 - ETRS89
            /// @accuracy Undefined / Varies. ~1-2 meters expected, but area of validity is specific to where parameters were derived.
            /// @area Romania - intended for Stereo 70 context.
            /// @method Assumed Coordinate Frame rotation (EPSG:9607 convention).
            /// @remarks **Non-standard parameters.** Origin unclear (potentially older standard, software default, or local fit). Similar
            /// parameters appear in various sources related to Romanian Stereo 70 transformations prior to official grid adoption. Use with
            /// caution. This variable was previously mislabeled as epsg3844 (which is a *projected CRS* code). The identical parameters
            /// previously listed under `epsg4178` (also a CRS code) are omitted here as redundant.
            /// @reference See discussion in FOSS4G 2014 slides (slide 18) -
            /// https://europe.foss4g.org/2014/slides/Daniel%20Urda%20-%20PROJ4%20Issues.The%20case%20of%20Romania2.pdf
            static constexpr helmert_params<T> pulkovo58_wgs84_ro_approx {
                .dx = T(33.4),       // meters
                .dy = T(-146.6),     // meters
                .dz = T(-76.3),      // meters
                .rx_sec = T(-0.359), // arcseconds
                .ry_sec = T(-0.053), // arcseconds
                .rz_sec = T(0.844),  // arcseconds
                .ds_ppm = T(-0.84),  // parts per million
                .name = "Non-EPSG Pulkovo 1942(58) / Stereo70 to WGS84 (Romania Approx)"};

            /// @brief Helmert parameters based on EPSG Transformation Code 1188.
            /// @name Pulkovo 1942 to WGS 84 (7)
            /// @epsg_code 1188
            /// @source_crs EPSG:4284 - Pulkovo 1942
            /// @target_crs EPSG:4326 - WGS 84
            /// @accuracy ~3 meters.
            /// @area Europe - FSU onshore W of 60°E; Albania; Bulgaria; Czechia; Germany - East; Hungary; Poland; Romania; Slovakia.
            /// @method Coordinate Frame rotation (EPSG:9607). Reverses CFr Bursa-Wolf (geocentric). Sign convention: P1 = P2 + T + D*R*P2.
            /// @remarks Valid transformation. Uses parameters defining the transformation FROM Pulkovo 1942 TO WGS 84. This variable was
            /// previously mislabeled as epsg4284 (which is the Pulkovo 1942 *geographic CRS* code).
            /// @reference https://epsg.org/transformation_1188/Pulkovo-1942-to-WGS-84-7.html
            static constexpr helmert_params<T> epsg1188 {
                .dx = T(23.92),     // meters
                .dy = T(-141.27),   // meters
                .dz = T(-80.9),     // meters
                .rx_sec = T(0.0),   // arcseconds (Note: Official EPSG has 0 here, check source if -0 intended)
                .ry_sec = T(-0.35), // arcseconds
                .rz_sec = T(0.82),  // arcseconds
                .ds_ppm = T(-0.12), // parts per million
                .name = "EPSG:1188 Pulkovo 1942 to WGS 84 (7)"};

            /// @brief Helmert parameters based on EPSG Transformation Code 15861.
            /// @name Pulkovo 1942(58) to WGS 84 (1)
            /// @epsg_code 15861
            /// @source_crs EPSG:4179 - Pulkovo 1942(58)
            /// @target_crs EPSG:4326 - WGS 84
            /// @accuracy ~1 meter.
            /// @area Europe - Former Soviet Union onshore; Afghanistan; Albania; Bulgaria; Czech Republic; Germany - East Germany; Hungary;
            /// Mongolia; Poland; Romania; Slovakia.
            /// @method Coordinate Frame rotation (EPSG:9607). Reverses CFr Bursa-Wolf (geocentric). Sign convention: P1 = P2 + T + D*R*P2.
            /// @remarks Valid transformation. Uses parameters defining the transformation FROM Pulkovo 1942(58) TO WGS 84. This
            /// transformation is relevant for workflows involving Stereo 70 in Romania, as Stereo 70 is based on Pulkovo 1942(58)
            /// (EPSG:4179).
            /// @reference https://epsg.org/transformation_15861/Pulkovo-1942-58-to-WGS-84-1.html
            static constexpr helmert_params<T> epsg15861 {.dx = T(24.986),     // meters
                                                          .dy = T(-128.922),   // meters
                                                          .dz = T(-83.764),    // meters
                                                          .rx_sec = T(-0.018), // arcseconds
                                                          .ry_sec = T(0.027),  // arcseconds
                                                          .rz_sec = T(0.854),  // arcseconds
                                                          .ds_ppm = T(0.09),   // parts per million
                                                          .name = "EPSG:15861 Pulkovo 1942(58) to WGS 84 (1)"};

            /// @brief Helmert parameters defined by ANCPI (Romanian Cadastre Agency) often used as an approximation for Stereo 70 <->
            /// ETRS89/WGS 84.
            /// @name Stereo 70 (Pulkovo 1942(58)) to ETRS89/WGS 84 (ANCPI approximation)
            /// @epsg_code N/A - Not an EPSG transformation code. Defined by national authority.
            /// @source_crs Implicitly EPSG:4179 - Pulkovo 1942(58) (used by Stereo 70 - EPSG:31700)
            /// @target_crs Implicitly EPSG:4258 - ETRS89 (often treated as equivalent to WGS 84 - EPSG:4326 for this accuracy)
            /// @accuracy Variable, represents an average fit. Better accuracy typically achieved using the official RONET grid
            /// transformation (NTv2 format). Likely ~0.5-1 meter average.
            /// @area Romania - specific to ANCPI regulations.
            /// @method Assumed Coordinate Frame rotation (EPSG:9607 convention). Parameters define Stereo70 -> ETRS89.
            /// @remarks **Non-standard EPSG parameters.** These values originate from Romanian national geodetic regulations (e.g., ANCPI
            /// Order 79/2010 or similar) and represent a Helmert approximation of the official grid-based transformation between Stereo 70
            /// and ETRS89. Use the official grid (`.gsb` file) for highest accuracy cadastral work in Romania.
            /// @reference ANCPI regulations (e.g., Order 79/2010). Also discussed in FOSS4G 2014 slides (slide 20) -
            /// https://europe.foss4g.org/2014/slides/Daniel%20Urda%20-%20PROJ4%20Issues.The%20case%20of%20Romania2.pdf
            static constexpr helmert_params<T> ancpi_stereo70_etrs89_approx {.dx = T(2.3283),          // meters
                                                                             .dy = T(-147.0416),       // meters
                                                                             .dz = T(-92.0802),        // meters
                                                                             .rx_sec = T(-0.30924979), // arcseconds
                                                                             .ry_sec = T(0.32482188),  // arcseconds
                                                                             .rz_sec = T(0.49730012),  // arcseconds
                                                                             .ds_ppm = T(5.68907711),  // parts per million
                                                                             .name = "ANCPI Stereo70 (Pulkovo58) to ETRS89/WGS84 Approx"};
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
        /// Initialized at compile time (or first use for static local). Error during init will throw.
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
                // Consider returning NaN or throwing if this is truly exceptional
                // return std::numeric_limits<T>::quiet_NaN();
            }

            // Term 1: tan(π/4 + φ/2)
            const T term1_arg = K::math::quarter_pi + lat_rad * K::values::half;
            const T term1 = std::tan(term1_arg);

            // Term 2: [(1 - e*sin(φ)) / (1 + e*sin(φ))]^(e/2)
            const T term2_num = K::values::one - esin_lat;
            const T term2_den = K::values::one + esin_lat;

            // Check for division by zero or negative base for power term
            // Denominator check: |1 + e*sin(lat)| near zero implies e*sin(lat) near -1
            if (std::abs(term2_den) < K::tolerance::near_zero)
                // This case is highly unlikely for valid e < 1, implies e is near 1 and lat is near -pi/2
                return std::numeric_limits<T>::quiet_NaN(); // Indicate calculation failure

            // Numerator check: (1 - e*sin(lat)) <= 0 implies e*sin(lat) >= 1, handled above.
            // Base = num/den. If base <= 0, power term is problematic.
            const T term2_base = term2_num / term2_den;
            if (term2_base <= K::tolerance::near_zero) // Check base is positive
                // This case is also highly unlikely for valid e < 1.
                return std::numeric_limits<T>::quiet_NaN(); // Indicate calculation failure

            const T term2 = std::pow(term2_base, e * K::values::half);

            // Argument for final atan: term1 * term2
            const T atan_arg = term1 * term2;

            // Check if atan argument became non-finite
            if (!std::isfinite(atan_arg))
                // If argument is infinite, atan approaches +/- pi/2
                return std::copysign(K::math::half_pi, atan_arg);

            // Final calculation: χ = 2 * atan(atan_arg) - π/2
            return K::values::two * std::atan(atan_arg) - K::math::half_pi;
        }

        // Core Transformation Functions

        /// @brief Converts Geodetic coordinates (lat, lon, altitude) to Geocentric Cartesian (X, Y, Z) coordinates.
        ///
        /// Uses the provided ellipsoid definition for the transformation.
        /// Formulae:
        ///   N = a / sqrt(1 - e^2 * sin^2(lat))
        ///   X = (N + h) * cos(lat) * cos(lon)
        ///   Y = (N + h) * cos(lat) * sin(lon)
        ///   Z = (N * (1 - e^2) + h) * sin(lat)
        ///
        /// @param[in] geo The input geodetic coordinate (latitude/longitude in degrees, altitude in meters).
        /// @param[in] ellip The reference ellipsoid parameters.
        /// @return The corresponding geocentric (ECEF) coordinate (X, Y, Z in meters).
        /// @throws std::invalid_argument if input `geo` coordinates are invalid (via `geo.validate()`).
        /// @throws std::runtime_error if the calculation encounters invalid states (e.g., division by zero due to N denominator).
        [[nodiscard]] static geocentric_coord<T> geodetic_to_geocentric(const geodetic<T>& geo, const ellipsoid<T>& ellip) noexcept(false)
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
            if (n_denominator_sq <= K::tolerance::near_zero) // Use tolerance, though strict <= 0 is often sufficient
                throw std::runtime_error(
                    "Invalid geodetic coordinate or ellipsoid parameter leading to non-positive N denominator squared: " +
                    std::to_string(n_denominator_sq));

            const T n_denominator = std::sqrt(n_denominator_sq);
            // No need to check n_denominator near zero again, covered by the check above.

            const T N = ellip.a / n_denominator;
            const T h = geo.altitude; // Ellipsoidal altitude

            // Calculate X, Y, Z
            const T N_plus_h = N + h;
            const T xy_factor = N_plus_h * cos_lat; // Common factor for X and Y

            const T x = xy_factor * std::cos(lon_rad);
            const T y = xy_factor * std::sin(lon_rad);
            // Z = (N * (1 - e^2) + h) * sin(lat) = (N * b^2/a^2 + h) * sin(lat)
            const T z = (N * (K::values::one - ellip.e2) + h) * sin_lat;

            // Construct and return result
            geocentric_coord<T> result {x, y, z};
            // Post-validation (optional, primarily checks for NaN/inf results)
            result.validate(); // Throws std::invalid_argument if non-finite
            return result;
        }

        /// @brief Converts Geocentric Cartesian (X, Y, Z) coordinates to Geodetic coordinates (lat, lon, altitude).
        ///
        /// Uses an iterative algorithm (e.g., Bowring's or similar) to find the latitude and altitude.
        /// Longitude is calculated directly from atan2(Y, X).
        ///
        /// @param[in] ecef The input geocentric (ECEF) coordinate (X, Y, Z in meters).
        /// @param[in] ellip The reference ellipsoid parameters.
        /// @param[in] tolerance The convergence tolerance for the latitude iteration (radians). Defaults to
        /// `K::tolerance::default_geodetic`.
        /// @param[in] max_iterations The maximum number of iterations allowed. Defaults to `K::tolerance::max_iterations`.
        /// @return The corresponding geodetic coordinate (latitude/longitude in degrees, altitude in meters). The coordinates are
        /// normalized.
        /// @throws std::invalid_argument if input `ecef` coordinates are invalid (via `ecef.validate()`).
        /// @throws std::runtime_error if the iteration fails to converge or encounters numerical instability.
        /// @note Handles special cases like points on the Z-axis (poles) and near the Earth's center.
        [[nodiscard]] static geodetic<T> geocentric_to_geodetic(const geocentric_coord<T>& ecef, const ellipsoid<T>& ellip,
                                                                const T tolerance = K::tolerance::default_geodetic,
                                                                const int max_iterations = K::tolerance::max_iterations) noexcept(false)
        {
            // Validate input coordinate first
            ecef.validate(); // Throws std::invalid_argument if invalid

            const T x = ecef.x;
            const T y = ecef.y;
            const T z = ecef.z;

            // Calculate distance from Z-axis
            const T p = std::hypot(x, y);

            geodetic<T> result; // Result coordinate (will be modified)

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
                    // Altitude relative to ellipsoid surface at origin. Closest point is pole at distance b.
                    // If we consider the origin itself projected onto the surface, altitude is negative b.
                    result.altitude = -ellip.b;
                    // Or maybe altitude is distance from origin? result.altitude = -std::hypot(p, z); is near zero.
                    // Using -b seems conventional for altitude *relative to surface*.
                }
                else
                {
                    // Point is on or very near the Z-axis (pole)
                    result.latitude = std::copysign(K::values::ninety, z); // +/- 90 degrees
                    // Altitude is distance along Z from polar surface point (distance |z| from origin)
                    result.altitude = std::abs(z) - ellip.b; // Altitude above the semi-minor axis endpoint
                }

                // No iteration needed, result is determined. Already normalized implicitly (lon=0).
                // Need to validate the calculated altitude, could be large negative.
                if (!std::isfinite(result.altitude)) // Should not happen if inputs are finite
                    throw std::runtime_error("Non-finite altitude calculated in pole case.");
                // result.validate(); // Optional final check (redundant with finite check)
                return result;
            }

            // General Case (Iterative Solution)

            // Calculate longitude directly
            result.longitude = std::atan2(y, x) * K::conversions::rad_to_deg;

            // Initial guess for latitude (based on spherical approximation or reduced latitude)
            // Using atan2(z, p * (1 - e^2)) is a common starting point (related to reduced latitude)
            // Avoid division by zero if p is extremely small but not caught by pole check (shouldn't happen)
            T lat_rad = std::atan2(z, p * (K::values::one - ellip.e2));
            T lat_rad_prev; // = lat_rad + tolerance * K::values::two; // Ensure first loop runs

            T N = ellip.a;         // Radius of curvature in prime vertical
            T h = K::values::zero; // Ellipsoidal altitude

            int i = 0;
            do
            {
                lat_rad_prev = lat_rad; // Store previous latitude for convergence check

                const T sin_lat = std::sin(lat_rad);
                const T cos_lat = std::cos(lat_rad);

                // Calculate N for the current latitude estimate
                const T N_denominator_sq = K::values::one - ellip.e2 * sin_lat * sin_lat;
                if (N_denominator_sq <= K::tolerance::near_zero) // Use tolerance here
                    throw std::runtime_error("Numerical instability: N denominator near or below zero during iteration.");

                N = ellip.a / std::sqrt(N_denominator_sq);

                // Calculate altitude h
                // Avoid division by cos_lat near poles (where p is not near zero)
                if (std::abs(cos_lat) < K::tolerance::near_zero)
                {
                    // Near pole: h = |z| / |sin(lat)| - N * (1 - e^2)
                    // Since sin(lat) will be close to +/- 1, this is safer. Check sin_lat != 0.
                    if (std::abs(sin_lat) < K::tolerance::near_zero)
                        // Should not happen if cos_lat is near zero unless lat is 0 (equator)
                        throw std::runtime_error("Numerical instability: Both sin(lat) and cos(lat) near zero.");

                    h = std::abs(z) / std::abs(sin_lat) - N * (K::values::one - ellip.e2);
                }
                else
                    h = p / cos_lat - N;

                // Update latitude estimate
                // Formula: lat = atan( Z / (P * (1 - e^2 * N / (N + h))) )
                // Avoid division by zero if (N+h) is near zero (point near center)
                const T N_plus_h = N + h;
                // Use a threshold relative to semi-major axis for near-center check
                const T center_threshold = K::tolerance::near_zero * ellip.a;

                if (std::abs(N_plus_h) < center_threshold)
                {
                    // Point is very close to the center, latitude is ill-defined.
                    // This might happen if iteration goes wrong or input is truly near center.
                    // Revert to a simple approximation based on Z and P.
                    // Let previous loop check handle convergence based on the last valid lat_rad.
                    // Or signal failure / stick to previous value?
                    throw std::runtime_error("Numerical instability: Point projects near ellipsoid center (N+h near zero).");
                    // Or maybe break and use lat_rad_prev? For now, throw.
                }

                const T lat_denom_factor = ellip.e2 * N / N_plus_h;
                // Check if |lat_denom_factor| >= 1 (should not happen if h > -N*(1-e^2))
                if (std::abs(lat_denom_factor) >= K::values::one)
                    throw std::runtime_error("Numerical instability: Invalid factor e^2*N/(N+h) >= 1.");

                const T lat_denom = p * (K::values::one - lat_denom_factor);

                if (std::abs(lat_denom) < center_threshold) // Check denominator too (implies point is near Z-axis/pole)
                {
                    // Denominator near zero implies point is near Z-axis (pole)
                    lat_rad = std::copysign(K::math::half_pi, z);
                    // Let convergence check handle if it's stable here.
                }
                else
                    lat_rad = std::atan2(z, lat_denom);

                // Check for convergence
                if (std::abs(lat_rad - lat_rad_prev) < tolerance)
                    break; // Converged

                i++;

            } while (i < max_iterations);

            // Check if max iterations were reached without convergence
            if (i >= max_iterations)
                throw std::runtime_error("Geocentric to Geodetic conversion failed to converge within maximum iterations.");

            // Store final results
            result.latitude = lat_rad * K::conversions::rad_to_deg;
            result.altitude = h;

            // Normalize the final coordinates (especially longitude, clamp latitude)
            result.normalize();
            // Final validation of calculated results
            result.validate(); // Throws if result is invalid (e.g., NaN, lat out of range)
            return result;
        }

        /// @brief Converts Krasovsky 1940 Geodetic coordinates to Stereo70 projected coordinates.
        ///
        /// Implements the forward Oblique Stereographic projection formula as defined for Stereo70 (EPSG:9809).
        /// Outputs coordinates according to EPSG:31700 axis definition: X = Northing, Y = Easting.
        ///
        /// @param[in] geo The input geodetic coordinate (lat/lon in degrees, altitude in meters) on the Krasovsky 1940 ellipsoid.
        /// @return The corresponding Stereo70 projected coordinate (X=Northing, Y=Easting, Z=Altitude) in meters.
        /// @throws std::invalid_argument if input `geo` coordinates are invalid (via `geo.validate()`).
        /// @throws std::runtime_error if the coordinate is too close to the antipodal point of the projection origin or if calculation
        /// fails.
        /// @sa stereo70_to_krasovsky_geodetic()
        [[nodiscard]] static stereo70::coordinate<T> krasovsky_geodetic_to_stereo70(const geodetic<T>& geo) noexcept(false)
        {
            // Validate input coordinate first
            geo.validate(); // Throws std::invalid_argument if invalid

            // Handle the origin case directly for precision and speed
            // Use tolerances defined in stereo70_params
            if (std::abs(geo.latitude - stereo70_params::lat0_deg) < stereo70_params::origin_tol_deg &&
                std::abs(geo.longitude - stereo70_params::lon0_deg) < stereo70_params::origin_tol_deg)
                // At the origin, result is False Northing (X), False Easting (Y)
                return {stereo70_params::fn, stereo70_params::fe, geo.altitude};

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
            // Normalize lambda_diff to (-pi, pi] range more robustly
            if (std::isfinite(lambda_diff))
            {
                lambda_diff = std::fmod(lambda_diff + K::math::pi, K::math::two_pi); // -> [0, 2pi) or (-2pi, 0)
                if (lambda_diff < K::values::zero)
                    lambda_diff += K::math::two_pi; // -> [0, 2pi)
                lambda_diff -= K::math::pi; // -> [-pi, pi) -- adjusted to match common convention, should be (-pi, pi]? Check usage.
                // Let's stick to atan2 range standard (-pi, pi]
                while (lambda_diff > K::math::pi)
                    lambda_diff -= K::math::two_pi;
                while (lambda_diff <= -K::math::pi)
                    lambda_diff += K::math::two_pi;
            }
            else
            {
                throw std::runtime_error("Non-finite longitude difference calculated.");
            }

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

            // Calculate standard Easting (x_std) and Northing (y_std)
            // x_std = k' * cos(chi) * sin(lambda)
            // y_std = k' * (cos(chi0)*sin(chi) - sin(chi0)*cos(chi)*cos(lambda))
            const T std_easting = k_prime * cos_chi * sin_lambda;
            const T std_northing = k_prime * (proj_params.cos_chi0 * sin_chi - proj_params.sin_chi0 * cos_chi * cos_lambda);

            // Apply false easting/northing according to EPSG:31700 convention
            // X = False Northing + standard Northing
            // Y = False Easting + standard Easting
            const T final_northing_X = stereo70_params::fn + std_northing;
            const T final_easting_Y = stereo70_params::fe + std_easting;

            // Assign to result struct according to EPSG:31700 (X=Northing, Y=Easting)
            stereo70::coordinate<T> result {final_northing_X, final_easting_Y,
                                            geo.altitude}; // Z coordinate is typically the ellipsoidal altitude

            // Validate the result (check for NaN/inf)
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
        /// Resulting Stereo70 coordinates follow EPSG:31700 convention (X=Northing, Y=Easting).
        ///
        /// @param[in] wgs84_coord The input WGS84 geodetic coordinate.
        /// @param[in] params The Helmert transformation parameters defining the transformation FROM the system
        ///                   underlying Krasovsky/Stereo70 (typically Pulkovo 1942/58) TO WGS84.
        ///                   Defaults to `conversion<T>::pulkovo58_to_wgs84`. The reverse transformation is applied.
        /// @return The corresponding Stereo70 projected coordinate (X=Northing, Y=Easting, Z=Altitude).
        /// @throws std::invalid_argument if input `wgs84_coord` or `params` are invalid.
        /// @throws std::runtime_error if any intermediate conversion step fails.
        /// @sa stereo70_to_wgs84()
        [[nodiscard]] static stereo70::coordinate<T> wgs84_to_stereo70(
            const wgs84::coordinate<T>& coord,
            const helmert_params<T>& params = conversion<T>::transformation::ancpi_stereo70_etrs89_approx) noexcept(false)
        {
            // Validate inputs
            params.validate(); // Throws if Helmert params invalid
            coord.validate();  // Throws if WGS84 coords invalid

            // Step 1: WGS84 Geodetic -> WGS84 ECEF
            const auto wgs84_ecef = geodetic_to_geocentric(coord, wgs84); // Can throw

            // Step 2: WGS84 ECEF -> Pulkovo ECEF (using REVERSE Helmert)
            const auto pulkovo_ecef = helmert_transformation_reverse(wgs84_ecef, params); // Can throw

            // Step 3: Pulkovo ECEF -> Krasovsky Geodetic
            // Note: Assumes Pulkovo ECEF is compatible with Krasovsky ellipsoid for geodetic conversion
            const auto krasovsky_geo = geocentric_to_geodetic(pulkovo_ecef, krasovsky1940); // Can throw

            // Step 4: Krasovsky Geodetic -> Stereo70 Projected
            // This function now correctly returns X=Northing, Y=Easting
            return krasovsky_geodetic_to_stereo70(krasovsky_geo); // Can throw
        }

        /// @brief Convenience operator.
        [[nodiscard]] stereo70::coordinate<T> operator()(
            const wgs84::coordinate<T>& coord,
            const helmert_params<T>& params = conversion<T>::transformation::dealul_piscului_1970_to_wgs84_epsg1838) noexcept(false)
        {
            return wgs84_to_stereo70(coord, params);
        }

        /// @brief Converts Stereo70 projected coordinates to WGS84 Geodetic coordinates.
        ///
        /// Assumes input Stereo70 coordinates follow EPSG:31700 convention (X=Northing, Y=Easting).
        /// This is the inverse composite transformation:
        /// 1. Stereo70 Projected -> Krasovsky Geodetic (`stereo70_to_krasovsky_geodetic`)
        /// 2. Krasovsky Geodetic -> Pulkovo ECEF (`geodetic_to_geocentric`)
        /// 3. Pulkovo ECEF -> WGS84 ECEF (using forward Helmert transformation) (`helmert_transformation_forward`)
        /// 4. WGS84 ECEF -> WGS84 Geodetic (`geocentric_to_geodetic`)
        ///
        /// @param[in] coord The input Stereo70 projected coordinate (X=Northing, Y=Easting, Z=Altitude).
        /// @param[in] params The Helmert transformation parameters defining the transformation FROM the system
        ///                   underlying Krasovsky/Stereo70 (typically Pulkovo 1942/58) TO WGS84.
        ///                   Defaults to `conversion<T>::pulkovo58_to_wgs84`. The forward transformation is applied.
        /// @return The corresponding WGS84 geodetic coordinate. The coordinates are normalized.
        /// @throws std::invalid_argument if input `stereo_coord` or `params` are invalid.
        /// @throws std::runtime_error if any intermediate conversion step fails.
        /// @sa wgs84_to_stereo70()
        [[nodiscard]] static wgs84::coordinate<T> stereo70_to_wgs84(
            const stereo70::coordinate<T>& coord,
            const helmert_params<T>& params = conversion<T>::transformation::ancpi_stereo70_etrs89_approx) noexcept(false)
        {
            // Validate inputs
            params.validate(); // Throws if Helmert params invalid
            coord.validate();  // Throws if Stereo70 coords invalid (finite check)

            // Step 1: Stereo70 Projected -> Krasovsky Geodetic
            // This function now correctly interprets input X=Northing, Y=Easting
            const auto krasovsky_geo = stereo70_to_krasovsky_geodetic(coord); // Can throw

            // Step 2: Krasovsky Geodetic -> Pulkovo ECEF
            // Note: Assumes Pulkovo ECEF is compatible with Krasovsky ellipsoid for geodetic conversion
            const auto pulkovo_ecef = geodetic_to_geocentric(krasovsky_geo, krasovsky1940); // Can throw

            // Step 3: Pulkovo ECEF -> WGS84 ECEF (using FORWARD Helmert)
            const auto wgs84_ecef = helmert_transformation_forward(pulkovo_ecef, params); // Calculation, should not throw if inputs valid

            // Step 4: WGS84 ECEF -> WGS84 Geodetic
            return geocentric_to_geodetic(wgs84_ecef, wgs84); // Can throw
        }

        /// @brief Convenience operator.
        [[nodiscard]] wgs84::coordinate<T> operator()(
            const stereo70::coordinate<T>& coord,
            const helmert_params<T>& params = conversion<T>::transformation::ancpi_stereo70_etrs89_approx) noexcept(false)
        {
            return stereo70_to_wgs84(coord, params);
        }

        /// @brief Converts Stereo70 projected coordinates (X=Northing, Y=Easting, Z=Altitude) to Krasovsky 1940 Geodetic coordinates.
        ///
        /// Implements the inverse Oblique Stereographic projection formula as defined for Stereo70 (EPSG:9809),
        /// using a series expansion to convert conformal latitude back to geodetic latitude.
        /// Expects input coordinates according to EPSG:31700 axis definition: X = Northing, Y = Easting.
        ///
        /// @param[in] proj The input Stereo70 projected coordinate (X=Northing, Y=Easting, Z=Altitude) in meters.
        /// @return The corresponding geodetic coordinate (lat/lon in degrees, altitude in meters) on the Krasovsky 1940 ellipsoid. The
        /// coordinates are normalized.
        /// @throws std::invalid_argument if input `proj` coordinates are invalid (via `proj.validate()`).
        /// @throws std::runtime_error if the calculation encounters numerical instability or invalid intermediate values.
        /// @sa krasovsky_geodetic_to_stereo70()
        [[nodiscard]] static geodetic<T> stereo70_to_krasovsky_geodetic(const projected_coord<T>& proj) noexcept(false)
        {
            // Validate input projected coordinate
            proj.validate(); // Throws std::invalid_argument if invalid

            // Calculate coordinates relative to false origin, respecting EPSG:31700 axis convention.
            // Standard formulas use x' (relative easting) and y' (relative northing).
            // Input proj.x is Northing, proj.y is Easting.
            const T y_prime = proj.x - stereo70_params::fn; // Relative Northing (matches y' in standard formulas)
            const T x_prime = proj.y - stereo70_params::fe; // Relative Easting (matches x' in standard formulas)

            // Calculate distance rho from the projection center (using relative easting/northing)
            const T rho = std::hypot(x_prime, y_prime);

            // Handle the origin case
            // Use tolerance defined in stereo70_params
            if (rho < stereo70_params::origin_proj_tol)
                // At the origin, return the projection origin coordinates
                // Use lat0_deg, lon0_deg from stereo70_params
                return {.latitude = stereo70_params::lat0_deg, .longitude = stereo70_params::lon0_deg, .altitude = proj.z};

            // Calculate intermediate angle 'c' (angular distance on conformal sphere)
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

            // Calculate conformal latitude (chi) of the point
            // Formula: sin(chi) = cos(c)*sin(chi0) + (y' * sin(c) * cos(chi0) / rho)
            // Check rho is not zero before dividing (handled by origin check, but safety)
            if (std::abs(rho) < K::tolerance::near_zero)
                // Should be caught by origin_proj_tol check
                throw std::runtime_error("Numerical instability: rho is near zero after origin check.");

            const T chi_arg_num = cos_c * proj_params.sin_chi0 + (y_prime * sin_c * proj_params.cos_chi0 / rho);
            // Clamp the argument to [-1, 1] before asin to handle potential floating point inaccuracies
            const T chi_arg = std::clamp(chi_arg_num, -K::values::one, K::values::one);
            const T chi = std::asin(chi_arg);
            if (!std::isfinite(chi)) // Check result of asin
                throw std::runtime_error("Failed to calculate valid conformal latitude (chi) in inverse projection.");

            // Calculate longitude difference (lambda) relative to the central meridian
            // Formula: tan(lambda) = (x' * sin(c)) / (rho * cos(chi0) * cos(c) - y' * sin(chi0) * sin(c))
            const T lambda_diff_num = x_prime * sin_c;
            const T lambda_diff_den = rho * proj_params.cos_chi0 * cos_c - y_prime * proj_params.sin_chi0 * sin_c;
            // Use atan2 for correct quadrant
            const T lambda_diff = std::atan2(lambda_diff_num, lambda_diff_den);
            if (!std::isfinite(lambda_diff)) // Check result of atan2
                throw std::runtime_error("Failed to calculate valid longitude difference (lambda_diff) in inverse projection.");

            // Calculate geodetic latitude (lat_rad) from conformal latitude (chi) using series expansion
            // lat = chi + C1*sin(2chi) + C2*sin(4chi) + C3*sin(6chi) + C4*sin(8chi)
            // Use pre-calculated coefficients c1, c2, c3, c4 from proj_params
            const T lat_rad = calculate_geodetic_latitude_series(chi, proj_params.c1, proj_params.c2, proj_params.c3, proj_params.c4);
            if (!std::isfinite(lat_rad)) // Check result of series calculation
                throw std::runtime_error("Failed to calculate valid geodetic latitude (lat_rad) from series expansion.");

            // Construct the result geodetic coordinate
            geodetic<T> result {
                .latitude = lat_rad * K::conversions::rad_to_deg, // Convert back to degrees
                .longitude =
                    (proj_params.lon0_rad + lambda_diff) * K::conversions::rad_to_deg, // Add central meridian back, convert to degrees
                .altitude = proj.z                                                     // Altitude is preserved
            };

            // Normalize and validate the final result
            result.normalize(); // Ensure lat/lon are in standard ranges
            result.validate();  // Throws if result is invalid (e.g., NaN, lat out of range)
            return result;
        }

    private:
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
        /// @param[in] params The Helmert parameters defining the transformation FROM source TO target. Must be validated beforehand.
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

            return {x_target, y_target, z_target};
        }

        /// @brief Applies a reverse 7-parameter Helmert transformation (Coordinate Frame method).
        ///
        /// Transforms coordinates FROM the target system TO the source system defined by `params` (e.g., WGS84 ECEF -> Pulkovo ECEF).
        /// Inverts the forward transformation.
        /// Formula (derived from forward, small angle approximation):
        ///   Let X' = (X_target - dx), Y' = (Y_target - dy), Z' = (Z_target - dz)
        ///   X_source = (X' + rz*Y' - ry*Z') / (1+ds) -> Incorrect, rotation needs to be applied after scaling
        ///   Correct approach:
        ///   X_s = (X_t - dx), Y_s = (Y_t - dy), Z_s = (Z_t - dz) ; remove translation
        ///   X_r = X_s / scale, Y_r = Y_s / scale, Z_r = Z_s / scale ; remove scale
        ///   X_source = X_r + rz*Y_r - ry*Z_r ; apply inverse rotation (transposed matrix for small angles)
        ///   Y_source = -rz*X_r + Y_r + rx*Z_r
        ///   Z_source = ry*X_r - rx*Y_r + Z_r
        ///
        /// @param[in] target The target geocentric coordinate (X, Y, Z) - the coordinate we are transforming FROM.
        /// @param[in] params The Helmert parameters defining the transformation FROM source TO target. Must be validated beforehand.
        /// @return The transformed geocentric coordinate in the source system.
        /// @throws std::runtime_error if the scale factor is too close to zero, preventing inversion (should be caught by param
        /// validation).
        /// @note Assumes `params` have been validated for finiteness and non-zero scale factor.
        /// @sa helmert_transformation_forward()
        [[nodiscard]] static geocentric_coord<T> helmert_transformation_reverse(
            const geocentric_coord<T>& target,
            const helmert_params<T>& params) noexcept(false) // Can technically still fail if scale is bad despite validation
        {
            // Get parameters in needed units
            const T rx = params.rx_rad();          // Rotation around X in radians
            const T ry = params.ry_rad();          // Rotation around Y in radians
            const T rz = params.rz_rad();          // Rotation around Z in radians
            const T scale = params.scale_factor(); // (1 + ds)

            // Check if scale factor is valid for division (redundant if params.validate() was called, but safe)
            if (std::abs(scale) < K::tolerance::near_zero)
                throw std::runtime_error("Invalid scale factor (near zero: " + std::to_string(scale) +
                                         ") in reverse Helmert transformation, cannot invert.");

            const T inv_scale = K::values::one / scale;

            // Reverse translation
            const T x_temp = target.x - params.dx;
            const T y_temp = target.y - params.dy;
            const T z_temp = target.z - params.dz;

            // Reverse scale
            const T x_scaled = x_temp * inv_scale;
            const T y_scaled = y_temp * inv_scale;
            const T z_scaled = z_temp * inv_scale;

            // Apply inverse rotation (transpose of rotation matrix for small angles approximation)
            // Note: This is the inverse of the Coordinate Frame rotation matrix used in forward.
            const T x_source = x_scaled + rz * y_scaled - ry * z_scaled;
            const T y_source = -rz * x_scaled + y_scaled + rx * z_scaled;
            const T z_source = ry * x_scaled - rx * y_scaled + z_scaled;

            // Result validation (check for NaN/Inf which might occur if intermediate values are extreme)
            geocentric_coord<T> result {x_source, y_source, z_source};

            // No validate method on geocentric_coord alias directly, but it's an xyz
            // result.validate(); // Add this check for robustness
            if (!std::isfinite(result.x) || !std::isfinite(result.y) || !std::isfinite(result.z))
                throw std::runtime_error("Non-finite result in reverse Helmert transformation.");

            return result;
        }
    }; // end struct conversion
}
