/// Copyright (c) 2025 - present KMX Systems. All rights reserved.
/// @file gis/coordinate/conversion.hpp
/// @brief Defines core constants, data structures, and conversion functions for GIS operations.
///
/// @note This file assumes a C++20 compliant compiler.
/// @warning Does not use PCH if PCH is not defined.
#pragma once
#ifndef PCH
    #include <cmath>    ///< Pentru std::isfinite
    #include <format>   ///< Pentru std::format (C++20, complet în C++23)
    #include <iomanip>  ///< For std::fixed, std::setprecision
    #include <iostream> ///< Pentru std::ostream, std::endl
    #include <iterator> ///< Pentru std::ostreambuf_iterator
    #include <kmx/gis/constants.hpp>
    #include <stdexcept> ///< For std::invalid_argument, std::runtime_error
    #include <stdexcept> ///< Pentru std::invalid_argument, std::runtime_error
    #include <string>    ///< Pentru std::string
#endif

/// @namespace kmx::gis
/// @brief Namespace for Geographic Information System (GIS) related functionalities.
///
/// This namespace contains structures, constants, and algorithms commonly used
/// in geodetic and cartographic computations.
namespace kmx::gis::coordinate
{
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
    struct geodetic
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

    // Forward declaration for the friend operator<<
    template <typename T>
    struct xy;

    template <typename T>
    std::ostream& operator<<(std::ostream& out, const xy<T>& coord);

    /// @brief Represents a generic 2D pair of values (X, Y).
    ///
    /// Can be used for spatial coordinates (if T is floating-point) or
    /// discrete indices (if T is integral).
    ///
    /// @tparam T The arithmetic type (e.g., `float`, `double`, `int`, `std::uint8_t`)
    ///           for the component values. Must satisfy `std::is_arithmetic_v<T>`.
    template <typename T>
    struct xy
    {
        /// @brief Static assertion to ensure T is an arithmetic type.
        static_assert(std::is_arithmetic_v<T>, "Template parameter T must be an arithmetic type for xy.");

        /// @brief The type of the values stored for each component.
        using value_type = T;

        value_type x {}; ///< The first component (e.g., X-coordinate or Column Index).
        value_type y {}; ///< The second component (e.g., Y-coordinate or Row Index).

        /// @brief Default constructor. Initializes x and y to zero.
        constexpr xy() noexcept = default;

        /// @brief Parameterized constructor to initialize components.
        /// @param x_val The value for the first component.
        /// @param y_val The value for the second component.
        constexpr xy(const value_type x_val, const value_type y_val) noexcept: x(x_val), y(y_val) {}

        /// @brief Validates that components are finite if T is a floating-point type.
        /// For integral types, this check is a no-op.
        /// @throws std::invalid_argument if any component (x, y) is not finite (for floating-point T).
        constexpr void validate() const noexcept(false)
        {
            if constexpr (std::is_floating_point_v<value_type>)
                if (!std::isfinite(x) || !std::isfinite(y))
                    throw std::invalid_argument("Floating-point xy components (x, y) must be finite");
        }

        /// @brief Formats the components as a string in the style "x, y".
        /// Adjusts formatting based on whether T is floating-point or integral.
        /// @param precision The number of decimal places for formatting floating-point values. Defaults to 5.
        ///                  Ignored for integral types.
        /// @return A std::string representing the components.
        [[nodiscard]] std::string to_string_simple(const int precision = 5) const
        {
            if constexpr (std::is_floating_point_v<value_type>)
            {
                const std::string fmt_string = std::format("{{:.{0}f}}, {{:.{0}f}}", precision);
                return std::vformat(fmt_string, std::make_format_args(x, y));
            }

            return std::format("{}, {}", x, y);
        }

        /// @brief Prints the components to the specified output stream.
        /// Adjusts formatting based on T.
        /// @param out The output stream (e.g., std::cout).
        /// @param precision The number of decimal places for floating-point T. Defaults to 5.
        /// @return A reference to the output stream.
        std::ostream& print(std::ostream& out, const int precision = 5) const
        {
            if constexpr (std::is_floating_point_v<value_type>)
            {
                const std::string fmt_string = std::format("{{:.{0}f}}, {{:.{0}f}}", precision);
                std::vformat_to(std::ostreambuf_iterator<char>(out), fmt_string, std::make_format_args(x, y));
            }
            else
            {
                if constexpr (sizeof(value_type) == 1 && (std::is_integral_v<value_type>) )
                    std::format_to(std::ostreambuf_iterator<char>(out), "{}, {}", static_cast<int>(x), static_cast<int>(y));
                else
                    std::format_to(std::ostreambuf_iterator<char>(out), "{}, {}", x, y);
            }

            return out;
        }

        /// @brief Prints the components to the stream, followed by a newline.
        /// @param out The output stream.
        /// @param precision The number of decimal places for formatting. Defaults to 5.
        /// @return A reference to the output stream.
        std::ostream& println(std::ostream& out, const int precision = 5) const { return print(out, precision) << std::endl; }

        /// @brief Defaulted spaceship operator for lexicographical comparisons.
        /// Compares first by x, then by y using their natural comparison.
        /// @param other The xy object to compare against.
        /// @return The result of the comparison (e.g., std::strong_ordering).
        constexpr auto operator<=>(const xy& other) const noexcept = default;

        /// @brief Explicit equality operator.
        /// Uses epsilon comparison for floating-point types, exact comparison for integral types.
        /// @param other The xy object to compare against.
        /// @return True if the objects are considered equal based on their type, false otherwise.
        constexpr bool operator==(const xy& other) const noexcept
        {
            if constexpr (std::is_floating_point_v<value_type>)
            {
                constexpr value_type type_specific_epsilon = []
                {
                    if constexpr (std::is_same_v<value_type, float>)
                        return 1e-6f;
                    else if constexpr (std::is_same_v<value_type, double>)
                        return 1e-9;
                    else if constexpr (std::is_same_v<value_type, long double>)
                        return 1e-12L;
                    else
                        return static_cast<value_type>(0); // Should not be reached for FP types
                }();
                return (std::abs(x - other.x) < type_specific_epsilon) && (std::abs(y - other.y) < type_specific_epsilon);
            }

            // Integral type
            return (x == other.x) && (y == other.y);
        }

        /// @brief Addition operator for two xy objects.
        /// @param other The xy object to add.
        /// @return A new xy object representing the sum of components.
        [[nodiscard]] constexpr xy operator+(const xy& other) const noexcept
        {
            return {static_cast<value_type>(x + other.x), static_cast<value_type>(y + other.y)};
        }

        /// @brief Addition assignment operator.
        /// @param other The xy object to add to this object.
        /// @return A reference to this modified object.
        constexpr xy& operator+=(const xy& other) noexcept
        {
            x = static_cast<value_type>(x + other.x);
            y = static_cast<value_type>(y + other.y);
            return *this;
        }

        /// @brief Subtraction operator for two xy objects.
        /// @param other The xy object to subtract.
        /// @return A new xy object representing the difference of components.
        [[nodiscard]] constexpr xy operator-(const xy& other) const noexcept
        {
            return {static_cast<value_type>(x - other.x), static_cast<value_type>(y - other.y)};
        }

        /// @brief Subtraction assignment operator.
        /// @param other The xy object to subtract from this object.
        /// @return A reference to this modified object.
        constexpr xy& operator-=(const xy& other) noexcept
        {
            x = static_cast<value_type>(x - other.x);
            y = static_cast<value_type>(y - other.y);
            return *this;
        }

        /// @brief Multiplication by a scalar operator.
        /// @param scalar The scalar value to multiply components by.
        /// @return A new xy object with scaled components.
        [[nodiscard]] constexpr xy operator*(const value_type scalar) const noexcept
        {
            return {static_cast<value_type>(x * scalar), static_cast<value_type>(y * scalar)};
        }

        /// @brief Multiplication by a scalar and assignment operator.
        /// @param scalar The scalar value to multiply components by.
        /// @return A reference to this modified object.
        constexpr xy& operator*=(const value_type scalar) noexcept
        {
            x = static_cast<value_type>(x * scalar);
            y = static_cast<value_type>(y * scalar);
            return *this;
        }

        /// @brief Division by a scalar operator.
        /// For integral types, this performs integer division.
        /// @param scalar The scalar value to divide components by.
        /// @throws std::runtime_error if scalar is zero.
        /// @return A new xy object with divided components.
        [[nodiscard]] constexpr xy operator/(const value_type scalar) const noexcept(false)
        {
            if (scalar == static_cast<value_type>(0))
                throw std::runtime_error("xy: Division by zero");
            return {static_cast<value_type>(x / scalar), static_cast<value_type>(y / scalar)};
        }
        /// @brief Division by a scalar and assignment operator.
        /// For integral types, this performs integer division.
        /// @param scalar The scalar value to divide components by.
        /// @throws std::runtime_error if scalar is zero.
        /// @return A reference to this modified object.
        constexpr xy& operator/=(const value_type scalar) noexcept(false)
        {
            if (scalar == static_cast<value_type>(0))
                throw std::runtime_error("xy: Division by zero (operator/=)");
            x = static_cast<value_type>(x / scalar);
            y = static_cast<value_type>(y / scalar);
            return *this;
        }
        /// @brief Unary minus (negation) operator.
        /// @return A new xy object with negated components.
        [[nodiscard]] constexpr xy operator-() const noexcept
        {
            // Note: Unary minus on unsigned types has defined behavior (modulo arithmetic),
            // but may not be intuitive for all use cases (e.g., indices).
            return {static_cast<value_type>(-x), static_cast<value_type>(-y)};
        }

        /// @brief Friend function for the stream insertion operator.
        /// Formats the output like "xy(x_val, y_val)" with default precision for floats.
        friend std::ostream& operator<< <T>(std::ostream& out, const xy<T>& coord);
    };

    /// @brief Overload of the stream insertion operator for xy coordinates.
    /// @tparam T The arithmetic type of the coordinates.
    /// @param out The output stream.
    /// @param coord The xy coordinate to print.
    /// @return A reference to the output stream.
    template <typename T>
    std::ostream& operator<<(std::ostream& out, const xy<T>& coord)
    {
        static_assert(std::is_arithmetic_v<T>, "Template parameter T must be an arithmetic type for xy operator<<.");
        const std::ostream::sentry sentry(out);
        if (sentry)
        {
            if constexpr (std::is_floating_point_v<T>)
            {
                const std::string fmt_string = std::format("xy({{:.{0}f}}, {{:.{0}f}})", 5); // Default precision 5
                std::vformat_to(std::ostreambuf_iterator<char>(out), fmt_string, std::make_format_args(coord.x, coord.y));
            }
            else
            { // Integral type
                if constexpr (sizeof(T) == 1 && (std::is_integral_v<T>) )
                    std::format_to(std::ostreambuf_iterator<char>(out), "xy({}, {})", static_cast<int>(coord.x), static_cast<int>(coord.y));
                else
                    std::format_to(std::ostreambuf_iterator<char>(out), "xy({}, {})", coord.x, coord.y);
            }
        }
        return out;
    }

    /// @brief Multiplication operator: scalar * coordinate (non-member form).
    /// @tparam T The floating-point type of the coordinates.
    /// @param scalar The scalar value.
    /// @param coord The xy coordinate.
    /// @return A new xy coordinate resulting from the multiplication.
    template <typename T>
    [[nodiscard]] constexpr xy<T> operator*(typename xy<T>::value_type scalar, const xy<T>& coord) noexcept
    {
        static_assert(std::is_arithmetic_v<T>, "Template parameter T must be an arithmetic type for xy operator*.");
        return coord * scalar; // Uses the member operator
    }

    // Forward declaration for the stream insertion operator
    template <typename T>
    struct xyz;

    template <typename T>
    std::ostream& operator<<(std::ostream& out, const xyz<T>& coord);

    /// @brief Represents a 3D Cartesian coordinate (X, Y, Z).
    ///
    /// This structure is used for various 3D Cartesian systems,
    /// such as ECEF or projected coordinates with altitude.
    /// Coordinates are typically in meters. Coordinates are mutable.
    /// Optimized for C++23.
    ///
    /// @tparam T The floating-point type (e.g., `float`, `double`, `long double`) for the coordinate values.
    ///           Must satisfy `std::is_floating_point_v<T>`.
    template <typename T>
    struct xyz
    {
        /// @brief Static assertion to ensure T is a floating-point type.
        static_assert(std::is_floating_point_v<T>, "Template parameter T must be a floating-point type for xyz.");

        /// @brief The type of the values stored for each coordinate component.
        using value_type = T;

        value_type x {}; ///< X-coordinate (e.g., meters in ECEF; Easting).
        value_type y {}; ///< Y-coordinate (e.g., meters in ECEF; Northing).
        value_type z {}; ///< Z-coordinate (e.g., meters in ECEF; Altitude/Elevation).

        /// @brief Default constructor. Initializes x, y, and z to zero (origin).
        constexpr xyz() noexcept = default;

        /// @brief Parameterized constructor to initialize coordinates.
        /// @param x_val The value for the X-coordinate.
        /// @param y_val The value for the Y-coordinate.
        /// @param z_val The value for the Z-coordinate.
        constexpr xyz(value_type x_val, value_type y_val, value_type z_val) noexcept: x(x_val), y(y_val), z(z_val) {}

        /// @brief Validates that the X, Y, and Z coordinates are finite numbers (not NaN or infinity).
        /// @throws std::invalid_argument if any coordinate (x, y, z) is not finite.
        constexpr void validate() const noexcept(false)
        {
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
                throw std::invalid_argument("XYZ coordinates (x, y, z) must be finite");
        }

        /// @brief Formats the coordinates as a string in the style "x, y, z".
        /// @param precision The number of decimal places for formatting the values. Defaults to 5.
        /// @return A std::string representing the coordinates.
        [[nodiscard]] std::string to_string_simple(int precision = 5) const
        {
            if (precision == 5)
                return std::format("{:.5f}, {:.5f}, {:.5f}", x, y, z);
            return std::vformat(std::format("{{:.{0}f}}, {{:.{0}f}}, {{:.{0}f}}", precision), std::make_format_args(x, y, z));
        }

        /// @brief Prints the coordinates to the specified output stream, in the style "x, y, z".
        /// Preserves the original format from the initial requirement.
        /// @param out The output stream (e.g., std::cout).
        /// @param precision The number of decimal places. Defaults to 5.
        /// @return A reference to the output stream, for chaining.
        std::ostream& print(std::ostream& out, int precision = 5) const
        {
            if (precision == 5)
                std::format_to(std::ostreambuf_iterator<char>(out), "{:.5f}, {:.5f}, {:.5f}", x, y, z);
            else
                std::vformat_to(std::ostreambuf_iterator<char>(out), std::format("{{:.{0}f}}, {{:.{0}f}}, {{:.{0}f}}", precision),
                                std::make_format_args(x, y, z));
            return out;
        }

        /// @brief Prints the coordinates to the specified output stream, followed by a newline.
        /// @param out The output stream.
        /// @param precision The number of decimal places. Defaults to 5.
        /// @return A reference to the output stream.
        std::ostream& println(std::ostream& out, int precision = 5) const { return print(out, precision) << std::endl; }

        /// @brief Spaceship operator for comparisons (generates ==, !=, <, <=, >, >=).
        /// Comparison is lexicographical: first by x, then by y, then by z.
        /// @param other The xyz coordinate to compare against.
        /// @return The result of the comparison (std::strong_ordering).
        constexpr auto operator<=>(const xyz& other) const noexcept = default;

        /// @brief Addition operator for two xyz coordinates.
        /// @param other The xyz coordinate to add.
        /// @return A new xyz coordinate representing the sum.
        [[nodiscard]] constexpr xyz operator+(const xyz& other) const noexcept { return {x + other.x, y + other.y, z + other.z}; }

        /// @brief Addition assignment operator.
        /// @param other The xyz coordinate to add.
        /// @return A reference to the current object (*this), modified.
        constexpr xyz& operator+=(const xyz& other) noexcept
        {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        /// @brief Subtraction operator for two xyz coordinates.
        /// @param other The xyz coordinate to subtract.
        /// @return A new xyz coordinate representing the difference.
        [[nodiscard]] constexpr xyz operator-(const xyz& other) const noexcept { return {x - other.x, y - other.y, z - other.z}; }

        /// @brief Subtraction assignment operator.
        /// @param other The xyz coordinate to subtract.
        /// @return A reference to the current object (*this), modified.
        constexpr xyz& operator-=(const xyz& other) noexcept
        {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        /// @brief Multiplication by a scalar operator.
        /// @param scalar The scalar value to multiply the components by.
        /// @return A new xyz coordinate with scaled components.
        [[nodiscard]] constexpr xyz operator*(value_type scalar) const noexcept { return {x * scalar, y * scalar, z * scalar}; }

        /// @brief Multiplication by a scalar and assignment operator.
        /// @param scalar The scalar value to multiply the components by.
        /// @return A reference to the current object (*this), modified.
        constexpr xyz& operator*=(value_type scalar) noexcept
        {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        /// @brief Division by a scalar operator.
        /// @param scalar The scalar value to divide the components by.
        /// @throws std::runtime_error if scalar is zero.
        /// @return A new xyz coordinate with divided components.
        [[nodiscard]] constexpr xyz operator/(value_type scalar) const noexcept(false)
        {
            if (scalar == static_cast<value_type>(0))
                throw std::runtime_error("xyz: Division by zero");

            return {x / scalar, y / scalar, z / scalar};
        }

        /// @brief Division by a scalar and assignment operator.
        /// @param scalar The scalar value to divide the components by.
        /// @throws std::runtime_error if scalar is zero.
        /// @return A reference to the current object (*this), modified.
        constexpr xyz& operator/=(value_type scalar) noexcept(false)
        {
            if (scalar == static_cast<value_type>(0))
                throw std::runtime_error("xyz: Division by zero (operator/=)");

            x /= scalar;
            y /= scalar;
            z /= scalar;
            return *this;
        }

        /// @brief Unary minus (negation) operator.
        /// @return A new xyz coordinate with negated components.
        [[nodiscard]] constexpr xyz operator-() const noexcept { return {-x, -y, -z}; }

        /// @brief Friend function for the stream insertion operator (e.g., std::cout << coord;).
        /// Formats the output as "xyz(x_val, y_val, z_val)".
        friend std::ostream& operator<< <T>(std::ostream& out, const xyz<T>& coord);
    };

    /// @brief Overload of the stream insertion operator for xyz coordinates.
    /// @tparam T The floating-point type of the coordinates.
    /// @param out The output stream.
    /// @param coord The xyz coordinate to print.
    /// @return A reference to the output stream.
    template <typename T>
    std::ostream& operator<<(std::ostream& out, const xyz<T>& coord)
    {
        static_assert(std::is_floating_point_v<T>, "Template parameter T must be a floating-point type for xyz operator<<.");
        std::ostream::sentry s(out); // Ensures the stream is in a valid state
        if (s)
        {
            // Descriptive formatting, directly to the stream
            std::format_to(std::ostreambuf_iterator<char>(out), "xyz({:.5f}, {:.5f}, {:.5f})", coord.x, coord.y, coord.z);
        }
        return out;
    }

    /// @brief Multiplication operator: scalar * coordinate (non-member form).
    /// @tparam T The floating-point type of the coordinates.
    /// @param scalar The scalar value.
    /// @param coord The xyz coordinate.
    /// @return A new xyz coordinate resulting from the multiplication.
    template <typename T>
    [[nodiscard]] constexpr xyz<T> operator*(typename xyz<T>::value_type scalar, const xyz<T>& coord) noexcept
    {
        static_assert(std::is_floating_point_v<T>, "Template parameter T must be a floating-point type for xyz scalar multiplication.");
        return coord * scalar; // Uses the member operator
    }

    /// @brief Calculates the absolute difference between two Cartesian coordinates.
    ///
    /// Computes the absolute difference for each corresponding component (x, y, z)
    /// of two xyz objects.
    ///
    /// @tparam T The underlying numeric type of the coordinates (e.g., float, double).
    /// @param a The first Cartesian coordinate object.
    /// @param b The second Cartesian coordinate object.
    /// @return A std::tuple containing the absolute differences in the order:
    ///         (|a.x - b.x|, |a.y - b.y|, |a.z - b.z|).
    /// @note This function is guaranteed not to throw exceptions.
    template <typename T>
    constexpr std::tuple<T, T, T> diff(const xyz<T>& a, const xyz<T>& b) noexcept
    {
        using std::abs;
        return {abs(a.x - b.x), abs(a.y - b.y), abs(a.z - b.z)};
    }

    /// @brief Calculates the absolute difference between two geodetic coordinates.
    ///
    /// Computes the absolute difference for each corresponding component (latitude,
    /// longitude, altitude) of two geodetic objects.
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
    constexpr std::tuple<T, T, T> diff(const geodetic<T>& a, const geodetic<T>& b) noexcept
    {
        using std::abs;
        return {abs(a.latitude - b.latitude), abs(a.longitude - b.longitude), abs(a.altitude - b.altitude)};
    }

    /// @brief Alias for `xyz<T>` representing geocentric (ECEF) coordinates.
    /// @tparam T The floating-point type for coordinate values.
    template <typename T>
    using geocentric_coord = xyz<T>;

    /// @brief Alias for `xyz<T>` representing projected (map) coordinates.
    ///
    /// Note: The meaning of X and Y depends on the specific Projected Coordinate Reference System (CRS).
    /// For EPSG:31700 (Stereo 70), X represents Northing and Y represents Easting.
    /// Z typically represents altitude/elevation above a reference surface.
    /// @tparam T The floating-point type for coordinate values.
    template <typename T>
    using projected_coord = xyz<T>;
} // namespace gis

/// @namespace kmx::gis::wgs84
/// @brief Contains types and parameters specific to the WGS84 datum and coordinate system.
namespace kmx::gis::wgs84
{
    /// @brief Alias for `geodetic<T>` representing coordinates in the WGS84 datum.
    /// @tparam T The floating-point type for coordinate values.
    template <typename T>
    using coordinate = coordinate::geodetic<T>;
} // namespace wgs84

// Specialization of std::hash for kmx::gis::coordinate::xy<T>
namespace std
{
    /// @brief Specialization of std::hash for the kmx::gis::coordinate::xy struct.
    /// Allows `xy` objects (with arithmetic T) to be used as keys in unordered containers.
    /// @tparam T The arithmetic type of the xy components.
    template <typename T>
    struct hash<kmx::gis::coordinate::xy<T>>
    {
        /// @brief Static assertion to ensure T is an arithmetic type.
        static_assert(std::is_arithmetic_v<T>, "T for hashing kmx::gis::coordinate::xy must be an arithmetic type.");

        /// @brief Calculates the hash value for a given xy object.
        /// @param p The xy object to hash.
        /// @return The calculated hash value as std::size_t.
        std::size_t operator()(const kmx::gis::coordinate::xy<T>& p) const noexcept
        {
            const std::size_t hx = std::hash<T> {}(p.x);
            const std::size_t hy = std::hash<T> {}(p.y);
            // Combine hashes
            return hx ^ (hy + 0x9e3779b9u + (hx << 6u) + (hx >> 2u));
        }
    };
} // namespace std
