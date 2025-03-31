// Copyright (c) 2025 - present KMX Systems. All rights reserved.
#ifndef PCH
    #include <algorithm>
    #include <cmath>
    #include <limits>
    #include <numbers> // requires C++2a
    #include <stdexcept>
    #include <type_traits>
#endif

namespace kmx::gis
{
    // Core Constants
    template <typename T = double>
    struct constants
    {
        static_assert(std::is_floating_point_v<T>, "constants requires a floating-point type.");

        // Avoid instantiation
        constants() = delete;

        // Common Values
        static constexpr T zero = static_cast<T>(0.0);
        static constexpr T point_five = static_cast<T>(0.5);
        static constexpr T one = static_cast<T>(1.0);
        static constexpr T two = static_cast<T>(2.0);
        static constexpr T three = static_cast<T>(3.0);
        static constexpr T four = static_cast<T>(4.0);
        static constexpr T five = static_cast<T>(5.0);
        static constexpr T six = static_cast<T>(6.0);
        static constexpr T seven = static_cast<T>(7.0);
        static constexpr T eight = static_cast<T>(8.0);
        static constexpr T nine = static_cast<T>(9.0);
        static constexpr T ten = static_cast<T>(10.0);
        static constexpr T twelve = static_cast<T>(12.0);
        static constexpr T thirteen = static_cast<T>(13.0);
        static constexpr T fifteen = static_cast<T>(15.0);
        static constexpr T twenty_four = static_cast<T>(24.0);
        static constexpr T twenty_nine = static_cast<T>(29.0);
        static constexpr T forty_eight = static_cast<T>(48.0);
        static constexpr T eighty_one = static_cast<T>(81.0);

        static constexpr T one_twenty = static_cast<T>(120.0);
        static constexpr T one_eighty = static_cast<T>(180.0);
        static constexpr T two_forty = static_cast<T>(240.0);
        static constexpr T three_sixty = static_cast<T>(360.0);
        static constexpr T one_one_two_zero = static_cast<T>(1120.0);             // Series
        static constexpr T one_one_five_two_zero = static_cast<T>(11520.0);       // Series
        static constexpr T one_six_one_two_eight_zero = static_cast<T>(161280.0); // Series

        static constexpr T million = static_cast<T>(1e6);

        // Mathematical Constants
        static constexpr T pi = std::numbers::pi_v<T>;
        static constexpr T two_pi = two * pi;
        static constexpr T half_pi = pi * point_five;
        static constexpr T quarter_pi = pi * point_five * point_five;

        // Specific small numbers for Helmert params
        static constexpr T point_three_five_nine = static_cast<T>(0.359);
        static constexpr T point_zero_five_three = static_cast<T>(0.053);
        static constexpr T point_eight_four_four = static_cast<T>(0.844);
        static constexpr T point_eight_four = static_cast<T>(0.84);

        static constexpr T sec_per_degree = static_cast<T>(3600.0); // Literal unavoidable unless defined elsewhere

        // Conversion Factors
        static constexpr T deg_to_rad = pi / one_eighty;
        static constexpr T rad_to_deg = one_eighty / pi;
        static constexpr T arcsec_to_rad = deg_to_rad / sec_per_degree;
        static constexpr T ppm_to_scale = one / million;

        // Tolerances & Limits
        static constexpr T default_geodetic_tolerance = static_cast<T>(1e-11);              // Literal unavoidable unless defined elsewhere
        static constexpr int default_max_iterations = 10;                                   // Integer literal is standard
        static constexpr T default_near_zero_threshold = std::numeric_limits<T>::epsilon(); // Use epsilon directly for robustness
    };

    // Core Structures
    template <typename T>
    struct ellipsoid
    {
        static_assert(std::is_floating_point_v<T>, "ellipsoid requires a floating-point type.");
        using K = constants<T>;

        T a;
        T inv_f;
        T f;
        T e2;
        T e;
        T b;
        T ep2;

        constexpr ellipsoid(const T semi_major, const T inverse_flattening):
            a(semi_major),
            inv_f(inverse_flattening),
            f(inv_f == K::zero ? K::zero : K::one / inv_f),
            e2(f == K::zero ? K::zero : K::two * f - f * f),
            e(std::sqrt(e2)),
            b(a * (K::one - f)),
            ep2(e2 == K::one ? std::numeric_limits<T>::infinity() : e2 / (K::one - e2))
        {
            if (semi_major <= K::zero)
                throw std::invalid_argument("Semi-major axis must be positive.");

            if (inverse_flattening < K::zero)
                throw std::invalid_argument("Inverse flattening cannot be negative.");

            if (K::one - e2 <= K::zero && e2 != K::one)
                throw std::runtime_error("Invalid ellipsoid parameters leading to non-positive (1-e^2).");
        }
    };

    template <typename T>
    struct helmert_params
    {
        static_assert(std::is_floating_point_v<T>, "helmert_params requires a floating-point type.");

        using K = constants<T>;

        T dx {};
        T dy {};
        T dz {};
        T rx_sec {};
        T ry_sec {};
        T rz_sec {};
        T ds_ppm {};

        constexpr T rx_pv_rad() const noexcept { return rx_sec * K::arcsec_to_rad; }
        constexpr T ry_pv_rad() const noexcept { return ry_sec * K::arcsec_to_rad; }
        constexpr T rz_pv_rad() const noexcept { return rz_sec * K::arcsec_to_rad; }
        constexpr T scale_factor() const noexcept { return K::one + ds_ppm * K::ppm_to_scale; }
    };

    template <typename T>
    struct geodetic_coord
    {
        static_assert(std::is_floating_point_v<T>);

        T latitude {};  // degrees
        T longitude {}; // degrees
        T height {};    // meters
    };

    template <typename T>
    struct xyz_coord // Base structure for Cartesian-like coordinates
    {
        static_assert(std::is_floating_point_v<T>);

        T x {}; // meters (Meaning depends on context: X, Easting)
        T y {}; // meters (Meaning depends on context: Y, Northing)
        T z {}; // meters (Meaning depends on context: Z, Height/Elevation)
    };

    // Specific Coordinate Type Aliases
    template <typename T>
    using geocentric_coord = xyz_coord<T>; // Represents ECEF (X, Y, Z)

    template <typename T>
    using projected_coord = xyz_coord<T>; // Represents projected (Easting, Northing, Height/Elevation)

    namespace wgs84
    {
        template <typename T = double>
        using coordinate = geodetic_coord<T>;
    }

    namespace stereo70
    {
        template <typename T = double>
        using coordinate = projected_coord<T>; // Stereo70 uses projected coords
    }

    // Conversion Logic
    template <typename T>
    struct conversion
    {
        using K = constants<T>;

        // Define Specific Parameters as Named Constants
        static constexpr T krasovsky_semi_major_axis = static_cast<T>(6378245.0);
        static constexpr T krasovsky_inverse_flattening = static_cast<T>(298.3);
        static constexpr T wgs84_semi_major_axis = static_cast<T>(6378137.0);
        static constexpr T wgs84_inverse_flattening = static_cast<T>(298.257223563);

        // Pulkovo 1942(58) -> WGS84 Helmert Parameters (Position Vector Convention assumed)
        static constexpr T pulkovo58_to_wgs84_dx = static_cast<T>(+33.4);
        static constexpr T pulkovo58_to_wgs84_dy = static_cast<T>(-146.6);
        static constexpr T pulkovo58_to_wgs84_dz = static_cast<T>(-76.3);
        static constexpr T pulkovo58_to_wgs84_rx_sec = static_cast<T>(-0.359);
        static constexpr T pulkovo58_to_wgs84_ry_sec = static_cast<T>(-0.053);
        static constexpr T pulkovo58_to_wgs84_rz_sec = static_cast<T>(+0.844);
        static constexpr T pulkovo58_to_wgs84_ds_ppm = static_cast<T>(-0.84);

        // Stereo70 Projection Parameters (based on Krasovsky 1940)
        static constexpr T stereo70_lat0_deg = static_cast<T>(46.0);
        static constexpr T stereo70_lon0_deg = static_cast<T>(25.0);
        static constexpr T stereo70_k0 = static_cast<T>(0.99975);
        static constexpr T stereo70_fe = static_cast<T>(500000.0);
        static constexpr T stereo70_fn = static_cast<T>(500000.0);

        // Tolerances specific to Stereo70 calculations
        static constexpr T stereo70_origin_coord_tolerance_deg = static_cast<T>(1e-9);
        static constexpr T stereo70_antipodal_threshold = static_cast<T>(1e-10);
        static constexpr T stereo70_origin_proj_tolerance_m = static_cast<T>(1e-6);

        // Ellipsoid Instances
        static constexpr ellipsoid<T> krasovsky1940_ellipsoid {krasovsky_semi_major_axis, krasovsky_inverse_flattening};
        static constexpr ellipsoid<T> wgs84_ellipsoid {wgs84_semi_major_axis, wgs84_inverse_flattening};

        // Helmert Parameter Instance
        static constexpr helmert_params<T> pulkovo58_to_wgs84_ellipsoid_params {.dx = pulkovo58_to_wgs84_dx,
                                                                                .dy = pulkovo58_to_wgs84_dy,
                                                                                .dz = pulkovo58_to_wgs84_dz,
                                                                                .rx_sec = pulkovo58_to_wgs84_rx_sec,
                                                                                .ry_sec = pulkovo58_to_wgs84_ry_sec,
                                                                                .rz_sec = pulkovo58_to_wgs84_rz_sec,
                                                                                .ds_ppm = pulkovo58_to_wgs84_ds_ppm};

        // Derived Constants for Stereo70 Projection
        struct derived_params
        {
            // Numerical constants specific to series expansion
            static constexpr T c2_term1_num = K::seven;
            static constexpr T c2_term1_den = K::forty_eight;
            static constexpr T c2_term2_num = K::twenty_nine;
            static constexpr T c2_term2_den = K::two_forty;
            static constexpr T c2_term3_num = static_cast<T>(811.0);
            static constexpr T c2_term3_den = K::one_one_five_two_zero;
            static constexpr T c3_term1_num = K::seven;
            static constexpr T c3_term1_den = K::one_twenty;
            static constexpr T c3_term2_num = K::eighty_one;
            static constexpr T c3_term2_den = K::one_one_two_zero;
            static constexpr T c4_term1_num = static_cast<T>(4279.0);
            static constexpr T c4_term1_den = K::one_six_one_two_eight_zero;

            T stereo70_lat0_rad;
            T stereo70_lon0_rad;
            T krasovsky_e;
            T krasovsky_e2;
            T krasovsky_e4;
            T krasovsky_e6;
            T krasovsky_e8;
            T sin_lat0;
            T cos_lat0;
            T chi0;
            T sin_chi0;
            T cos_chi0;
            T n0;
            T inv_two_n0;
            T c1, c2, c3, c4;
            bool valid;

            constexpr derived_params():
                stereo70_lat0_rad(stereo70_lat0_deg * K::deg_to_rad),
                stereo70_lon0_rad(stereo70_lon0_deg * K::deg_to_rad),
                krasovsky_e(krasovsky1940_ellipsoid.e),
                krasovsky_e2(krasovsky1940_ellipsoid.e2),
                krasovsky_e4(krasovsky_e2 * krasovsky_e2),
                krasovsky_e6(krasovsky_e4 * krasovsky_e2),
                krasovsky_e8(krasovsky_e6 * krasovsky_e2),
                sin_lat0(std::sin(stereo70_lat0_rad)),
                cos_lat0(std::cos(stereo70_lat0_rad)),
                chi0(calculate_conformal_latitude(stereo70_lat0_rad, krasovsky_e)),
                sin_chi0(std::sin(chi0)),
                cos_chi0(std::cos(chi0)),
                n0(calculate_n0(krasovsky1940_ellipsoid, chi0, sin_lat0, cos_lat0, cos_chi0)),
                inv_two_n0(calculate_inv_two_n0(n0)),
                c1(krasovsky_e2 * K::point_five + (K::five / K::twenty_four) * krasovsky_e4 + krasovsky_e6 / K::twelve +
                   (K::thirteen / K::three_sixty) * krasovsky_e8),
                c2((c2_term1_num / c2_term1_den) * krasovsky_e4 + (c2_term2_num / c2_term2_den) * krasovsky_e6 +
                   (c2_term3_num / c2_term3_den) * krasovsky_e8),
                c3((c3_term1_num / c3_term1_den) * krasovsky_e6 + (c3_term2_num / c3_term2_den) * krasovsky_e8),
                c4((c4_term1_num / c4_term1_den) * krasovsky_e8),
                valid(check())
            {
            }

            static constexpr T calculate_conformal_latitude(const T lat_rad, const T ellipsoid_e) noexcept
            {
                const T sin_lat = std::sin(lat_rad);
                const T esin_lat = ellipsoid_e * sin_lat;
                if (std::abs(esin_lat) >= K::one)
                    return std::copysign(K::half_pi, lat_rad);

                const T term1_arg = K::quarter_pi + lat_rad * K::point_five;
                const T term1 = std::tan(term1_arg);
                const T term2_num = K::one - esin_lat;
                const T term2_den = K::one + esin_lat;
                if (std::abs(term2_den) < K::default_near_zero_threshold || term2_num / term2_den <= K::zero)
                    return std::numeric_limits<T>::quiet_NaN();

                const T term2_base = term2_num / term2_den;
                const T term3_pow_arg = ellipsoid_e * K::point_five;
                const T term3 = std::pow(term2_base, term3_pow_arg);
                const T atan_arg = term1 * term3;
                if (!std::isfinite(atan_arg))
                    return std::copysign(K::half_pi, atan_arg);

                return K::two * std::atan(atan_arg) - K::half_pi;
            }

            static constexpr T calculate_n0(const ellipsoid<T>& ellip, const T chi0_val, const T sin_lat0_val, const T cos_lat0_val,
                                            const T cos_chi0_val) noexcept
            {
                if (!std::isfinite(chi0_val) || !std::isfinite(cos_chi0_val))
                    return std::numeric_limits<T>::quiet_NaN();

                const T e2_ = ellip.e2;
                const T w1_sq = K::one - e2_ * sin_lat0_val * sin_lat0_val;
                if (w1_sq <= K::zero)
                    return std::numeric_limits<T>::quiet_NaN();

                const T sqrt_w1 = std::sqrt(w1_sq);
                const T denominator = sqrt_w1 * cos_chi0_val;
                if (std::abs(denominator) < K::default_near_zero_threshold)
                    return std::numeric_limits<T>::quiet_NaN();

                return stereo70_k0 * ellip.a * cos_lat0_val / denominator;
            }

            static constexpr T calculate_inv_two_n0(const T n0_val) noexcept
            {
                if (!std::isfinite(n0_val) || std::abs(n0_val) < K::default_near_zero_threshold)
                    return std::numeric_limits<T>::quiet_NaN();

                return K::one / (K::two * n0_val);
            }

            constexpr bool check() const noexcept
            {
                using namespace std;
                return isfinite(stereo70_lat0_rad) && isfinite(stereo70_lon0_rad) && isfinite(krasovsky_e) && isfinite(krasovsky_e2) &&
                       isfinite(sin_lat0) && isfinite(cos_lat0) && isfinite(chi0) && isfinite(sin_chi0) && isfinite(cos_chi0) &&
                       isfinite(n0) && abs(n0) > K::default_near_zero_threshold && isfinite(inv_two_n0) && isfinite(c1) && isfinite(c2) &&
                       isfinite(c3) && isfinite(c4);
            }
        };

        static constexpr derived_params proj_params {};

        // Coordinate Transformation Functions

        // Geodetic (lat, lon, h) -> Geocentric Cartesian (X, Y, Z)
        static constexpr geocentric_coord<T> geodetic_to_geocentric(const geodetic_coord<T>& geo_coord, const ellipsoid<T>& ellip)
        {
            const T lat_rad = geo_coord.latitude * K::deg_to_rad;
            const T lon_rad = geo_coord.longitude * K::deg_to_rad;
            const T h = geo_coord.height;
            const T sin_lat = std::sin(lat_rad);
            const T cos_lat = std::cos(lat_rad);

            const T n_denominator_sq = K::one - ellip.e2 * sin_lat * sin_lat;
            if (n_denominator_sq <= K::zero)
                throw std::runtime_error("Invalid N denominator (geodetic_to_geocentric).");

            const T N = ellip.a / std::sqrt(n_denominator_sq);

            const T x = (N + h) * cos_lat * std::cos(lon_rad);
            const T y = (N + h) * cos_lat * std::sin(lon_rad);
            const T z = (N * (K::one - ellip.e2) + h) * sin_lat;
            return {.x = x, .y = y, .z = z}; // Returns geocentric_coord<T>
        }

        // Geocentric Cartesian (X, Y, Z) -> Geodetic (lat, lon, h)
        static constexpr geodetic_coord<T> geocentric_to_geodetic(const geocentric_coord<T>& geo_coord, const ellipsoid<T>& ellip,
                                                                  const T tolerance = K::default_geodetic_tolerance,
                                                                  const int max_iterations = K::default_max_iterations)
        {
            const T x = geo_coord.x;
            const T y = geo_coord.y;
            const T z = geo_coord.z;
            const T p = std::hypot(x, y);
            geodetic_coord<T> result_coord; // Renamed from geo_coord to avoid confusion with input type
            result_coord.longitude = K::zero;

            const T pole_threshold = K::default_near_zero_threshold * ellip.a;
            if (p < pole_threshold)
            {
                if (std::abs(z) < pole_threshold)
                {
                    result_coord.latitude = K::zero;
                    result_coord.height = -ellip.a;
                }
                else
                {
                    result_coord.latitude = std::copysign(K::half_pi, z) * K::rad_to_deg;
                    result_coord.height = std::abs(z) - ellip.b;
                }

                result_coord.longitude = K::zero;
                return result_coord;
            }

            result_coord.longitude = std::atan2(y, x) * K::rad_to_deg;
            T lat_rad = std::atan2(z, p * (K::one - ellip.e2));
            T lat_rad_prev = lat_rad;
            T N = ellip.a;

            for (int i = 0; i < max_iterations; ++i)
            {
                const T sin_lat = std::sin(lat_rad);
                const T cos_lat = std::cos(lat_rad);
                const T N_denominator_sq = K::one - ellip.e2 * sin_lat * sin_lat;
                if (N_denominator_sq <= K::zero)
                    throw std::runtime_error("Invalid N denominator during iteration (geocentric_to_geodetic).");

                N = ellip.a / std::sqrt(N_denominator_sq);

                if (std::abs(cos_lat) < K::default_near_zero_threshold)
                {
                    result_coord.height = std::abs(z) - ellip.b;
                    lat_rad = std::copysign(K::half_pi, z);
                    break;
                }
                else
                {
                    result_coord.height = (p / cos_lat) - N;
                }

                const T h_plus_N = N + result_coord.height;
                if (std::abs(h_plus_N) < K::default_near_zero_threshold)
                    lat_rad = (std::abs(z) < K::default_near_zero_threshold) ?
                                  K::zero :
                                  std::atan2(z, p * (K::one - ellip.e2 * N / h_plus_N)); // Let atan2 handle potential issues
                else
                {
                    const T lat_denom = p * (K::one - ellip.e2 * N / h_plus_N);
                    if (std::abs(lat_denom) < K::default_near_zero_threshold)
                    {
                        lat_rad = (std::abs(z) < K::default_near_zero_threshold) ? K::zero : std::copysign(K::half_pi, z);
                    }
                    else
                    {
                        lat_rad = std::atan2(z, lat_denom);
                    }
                }

                if (std::abs(lat_rad - lat_rad_prev) < tolerance)
                    break;

                lat_rad_prev = lat_rad;
                if (i == max_iterations - 1)
                    throw std::runtime_error("Geocentric to geodetic conversion did not converge.");
            }

            result_coord.latitude = lat_rad * K::rad_to_deg;
            while (result_coord.longitude > K::one_eighty)
                result_coord.longitude -= K::three_sixty;

            while (result_coord.longitude <= -K::one_eighty)
                result_coord.longitude += K::three_sixty;

            return result_coord;
        }

        // Helmert Transformations (operate on geocentric coordinates)
        static constexpr geocentric_coord<T> helmert_transformation_forward(const geocentric_coord<T>& source_coord,
                                                                            const helmert_params<T>& params)
        {
            const T rx_pv_rad = params.rx_pv_rad();
            const T ry_pv_rad = params.ry_pv_rad();
            const T rz_pv_rad = params.rz_pv_rad();
            const T scale_factor = params.scale_factor();
            const T x_src = source_coord.x;
            const T y_src = source_coord.y;
            const T z_src = source_coord.z;

            const T rotated_x = x_src - rz_pv_rad * y_src + ry_pv_rad * z_src;
            const T rotated_y = rz_pv_rad * x_src + y_src - rx_pv_rad * z_src;
            const T rotated_z = -ry_pv_rad * x_src + rx_pv_rad * y_src + z_src;

            return {.x = params.dx + scale_factor * rotated_x,
                    .y = params.dy + scale_factor * rotated_y,
                    .z = params.dz + scale_factor * rotated_z};
        }

        static constexpr geocentric_coord<T> helmert_transformation_reverse(const geocentric_coord<T>& target_coord,
                                                                            const helmert_params<T>& params)
        {
            const T rx_pv_rad = params.rx_pv_rad();
            const T ry_pv_rad = params.ry_pv_rad();
            const T rz_pv_rad = params.rz_pv_rad();
            const T forward_scale_factor = params.scale_factor();
            const T reverse_scale = (std::abs(forward_scale_factor) < K::default_near_zero_threshold) ?
                                        std::numeric_limits<T>::quiet_NaN() :
                                        K::one / forward_scale_factor;
            if (!std::isfinite(reverse_scale))
                throw std::runtime_error("Invalid scale factor in reverse Helmert transformation.");

            const T temp_x = target_coord.x - params.dx;
            const T temp_y = target_coord.y - params.dy;
            const T temp_z = target_coord.z - params.dz;

            const T rotated_x_corrected = temp_x + rz_pv_rad * temp_y - ry_pv_rad * temp_z;
            const T rotated_y_corrected = -rz_pv_rad * temp_x + temp_y + rx_pv_rad * temp_z;
            const T rotated_z_corrected = ry_pv_rad * temp_x - rx_pv_rad * temp_y + temp_z;

            return {.x = reverse_scale * rotated_x_corrected,
                    .y = reverse_scale * rotated_y_corrected,
                    .z = reverse_scale * rotated_z_corrected};
        }

        // Stereo70 Projection (operates between geodetic and projected coordinates)
        static constexpr projected_coord<T> krasovsky_geodetic_to_stereo70(
            const geodetic_coord<T>& krasovsky_coord) // Return type is projected_coord
        {
            if (!proj_params.valid)
                throw std::runtime_error("Projection setup error: Derived parameters invalid.");

            const T lat_deg_proj = krasovsky_coord.latitude;
            const T lon_deg_proj = krasovsky_coord.longitude;

            if (std::abs(lat_deg_proj - stereo70_lat0_deg) < stereo70_origin_coord_tolerance_deg &&
                std::abs(lon_deg_proj - stereo70_lon0_deg) < stereo70_origin_coord_tolerance_deg)
                return {.x = stereo70_fe, .y = stereo70_fn, .z = krasovsky_coord.height};

            const T lat_rad = lat_deg_proj * K::deg_to_rad;
            const T lon_rad = lon_deg_proj * K::deg_to_rad;
            const T chi = derived_params::calculate_conformal_latitude(lat_rad, proj_params.krasovsky_e);
            if (!std::isfinite(chi))
                throw std::runtime_error("Forward projection: calculate_conformal_latitude failed.");

            const T sin_chi = std::sin(chi);
            const T cos_chi = std::cos(chi);
            T lambda_diff = lon_rad - proj_params.stereo70_lon0_rad;
            while (lambda_diff > K::pi)
                lambda_diff -= K::two_pi;

            while (lambda_diff <= -K::pi)
                lambda_diff += K::two_pi;

            const T sin_lambda_diff = std::sin(lambda_diff);
            const T cos_lambda_diff = std::cos(lambda_diff);

            const T b_denominator = K::one + proj_params.sin_chi0 * sin_chi + proj_params.cos_chi0 * cos_chi * cos_lambda_diff;
            if (std::abs(b_denominator) < stereo70_antipodal_threshold)
                throw std::runtime_error("Forward projection: Input coordinate near antipodal point.");

            const T k_prime_factor = K::two * proj_params.n0 / b_denominator;

            const T x_proj_rel = k_prime_factor * cos_chi * sin_lambda_diff;
            const T y_proj_rel = k_prime_factor * (proj_params.cos_chi0 * sin_chi - proj_params.sin_chi0 * cos_chi * cos_lambda_diff);

            const projected_coord<T> proj_coord = {
                .x = stereo70_fe + x_proj_rel, .y = stereo70_fn + y_proj_rel, .z = krasovsky_coord.height};
            if (!std::isfinite(proj_coord.x) || !std::isfinite(proj_coord.y))
                throw std::runtime_error("Forward projection: Non-finite result calculated.");

            return proj_coord;
        }

        static constexpr T calculate_geodetic_latitude_series(const T chi, const T c1, const T c2, const T c3, const T c4) noexcept
        {
            const T sin_2chi = std::sin(K::two * chi);
            const T sin_4chi = std::sin(K::four * chi);
            const T sin_6chi = std::sin(K::six * chi);
            const T sin_8chi = std::sin(K::eight * chi);
            return chi + c1 * sin_2chi + c2 * sin_4chi + c3 * sin_6chi + c4 * sin_8chi;
        }

        static constexpr geodetic_coord<T> stereo70_to_krasovsky_geodetic(
            const projected_coord<T>& proj_coord) // Input type is projected_coord
        {
            if (!proj_params.valid)
                throw std::runtime_error("Projection setup error: Derived parameters invalid.");

            const T easting = proj_coord.x;
            const T northing = proj_coord.y;
            const T x_rel = easting - stereo70_fe;
            const T y_rel = northing - stereo70_fn;
            const T rho = std::hypot(x_rel, y_rel);

            if (rho < stereo70_origin_proj_tolerance_m)
                return {.latitude = stereo70_lat0_deg, .longitude = stereo70_lon0_deg, .height = proj_coord.z};

            if (!std::isfinite(proj_params.inv_two_n0))
                throw std::runtime_error("Inverse projection: inv_two_n0 is invalid.");

            const T tan_c_half = rho * proj_params.inv_two_n0;
            if (!std::isfinite(tan_c_half))
                throw std::runtime_error("Inverse projection: tan(c/2) calculation failed (non-finite).");

            const T c_angle = K::two * std::atan(tan_c_half);
            const T sin_c = std::sin(c_angle);
            const T cos_c = std::cos(c_angle);
            const T rho_inv = (rho > K::default_near_zero_threshold) ? K::one / rho : K::zero;
            T chi_inv_asin_arg_unclamped = cos_c * proj_params.sin_chi0 + (y_rel * sin_c * proj_params.cos_chi0 * rho_inv);
            const T chi_inv_asin_arg = std::clamp(chi_inv_asin_arg_unclamped, -K::one, K::one);
            const T chi_inv = std::asin(chi_inv_asin_arg);

            const T lambda_diff_inv_atan2_y = x_rel * sin_c;
            const T lambda_diff_inv_atan2_x = rho * proj_params.cos_chi0 * cos_c - y_rel * proj_params.sin_chi0 * sin_c;
            const T lambda_diff_inv = std::atan2(lambda_diff_inv_atan2_y, lambda_diff_inv_atan2_x);

            const T lat_rad = calculate_geodetic_latitude_series(chi_inv, proj_params.c1, proj_params.c2, proj_params.c3, proj_params.c4);
            if (!std::isfinite(lat_rad))
                throw std::runtime_error("Inverse projection: calculate_geodetic_latitude_series failed.");

            const T lon_rad = proj_params.stereo70_lon0_rad + lambda_diff_inv;
            geodetic_coord<T> geo_coord {.latitude = lat_rad * K::rad_to_deg, .longitude = lon_rad * K::rad_to_deg, .height = proj_coord.z};
            while (geo_coord.longitude > K::one_eighty)
                geo_coord.longitude -= K::three_sixty;

            while (geo_coord.longitude <= -K::one_eighty)
                geo_coord.longitude += K::three_sixty;

            if (!std::isfinite(geo_coord.latitude) || !std::isfinite(geo_coord.longitude))
                throw std::runtime_error("Inverse projection: Non-finite result after conversion/normalization.");

            return geo_coord;
        }
    }; // end struct conversion

    /// @brief Converts WGS84 geographic coordinate to Stereo70 projected coordinate.
    template <typename T = double>
    constexpr stereo70::coordinate<T> wgs84_to_stereo70( // Return type is stereo70::coordinate -> projected_coord
        const wgs84::coordinate<T>& wgs84_coord,         // Input type is wgs84::coordinate -> geodetic_coord
        const helmert_params<T>& pulkovo_to_wgs84_params = conversion<T>::pulkovo58_to_wgs84_ellipsoid_params)
    {
        using conv = conversion<T>;
        // 1. Geodetic -> Geocentric (WGS84)
        const geocentric_coord<T> wgs84_geocentric = conv::geodetic_to_geocentric(wgs84_coord, conv::wgs84_ellipsoid);
        // 2. Geocentric (WGS84) -> Geocentric (Pulkovo) via reverse Helmert
        const geocentric_coord<T> pulkovo_geocentric = conv::helmert_transformation_reverse(wgs84_geocentric, pulkovo_to_wgs84_params);
        // 3. Geocentric (Pulkovo) -> Geodetic (Krasovsky)
        const geodetic_coord<T> krasovsky_geodetic = conv::geocentric_to_geodetic(pulkovo_geocentric, conv::krasovsky1940_ellipsoid);
        // 4. Geodetic (Krasovsky) -> Projected (Stereo70)
        return conv::krasovsky_geodetic_to_stereo70(krasovsky_geodetic); // Returns projected_coord<T>
    }

    /// @brief Converts Stereo70 projected coordinate to WGS84 geographic coordinate.
    template <typename T = double>
    constexpr wgs84::coordinate<T> stereo70_to_wgs84( // Return type is wgs84::coordinate -> geodetic_coord
        const stereo70::coordinate<T>& stereo_coord,  // Input type is stereo70::coordinate -> projected_coord
        const helmert_params<T>& pulkovo_to_wgs84_params = conversion<T>::pulkovo58_to_wgs84_ellipsoid_params)
    {
        using conv = conversion<T>;
        // 1. Projected (Stereo70) -> Geodetic (Krasovsky)
        const geodetic_coord<T> krasovsky_geodetic = conv::stereo70_to_krasovsky_geodetic(stereo_coord);
        // 2. Geodetic (Krasovsky) -> Geocentric (Pulkovo)
        const geocentric_coord<T> pulkovo_geocentric = conv::geodetic_to_geocentric(krasovsky_geodetic, conv::krasovsky1940_ellipsoid);
        // 3. Geocentric (Pulkovo) -> Geocentric (WGS84) via forward Helmert
        const geocentric_coord<T> wgs84_geocentric = conv::helmert_transformation_forward(pulkovo_geocentric, pulkovo_to_wgs84_params);
        // 4. Geocentric (WGS84) -> Geodetic (WGS84)
        return conv::geocentric_to_geodetic(wgs84_geocentric, conv::wgs84_ellipsoid); // Returns geodetic_coord<T>
    }

} // namespace kmx::gis
