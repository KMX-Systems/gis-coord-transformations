/// Copyright (c) 2025 - present KMX Systems. All rights reserved.
#pragma once
#ifndef PCH
    #include <array>       ///< For std::array
    #include <functional>  ///< For std::reference_wrapper, std::cref
    #include <iostream>    ///< For std::ostream, std::endl
    #include <kmx/gis.hpp> ///< Core GIS library header
    #include <stdexcept>   ///< For std::invalid_argument, std::out_of_range, std::runtime_error
    #include <string>      ///< For std::string, std::stod, std::to_string
    #include <string_view> ///< For std::string_view
#endif

// Define the application-specific components within the kmx::gis namespace
namespace kmx::gis
{
    /// @brief Performs coordinate conversions using multiple predefined Helmert transformations.
    ///
    /// This class takes command-line arguments representing an input coordinate, converts it
    /// to an output coordinate type using several different sets of standard Helmert parameters
    /// (defined within `gis::conversion<T>::transformation`), and prints the results along
    /// with the differences between the results obtained using different parameter sets.
    /// It handles argument parsing and orchestrates the conversion process for a fixed set of Helmert parameters.
    ///
    /// @tparam InputCoordT The type of the input coordinate (e.g., `wgs84::coordinate<double>`). Must have a `value_type` member and a
    /// constructor taking (x, y, z).
    /// @tparam OutputCoordT The type of the output coordinate (e.g., `stereo70::coordinate<double>`). Must be default constructible and
    /// assignable. Must have a `print` method.
    template <typename InputCoordT, typename OutputCoordT>
    class coordinate_converter
    {
    private:
        /// @brief The underlying floating-point type derived from the InputCoordT.
        using value_type = typename InputCoordT::value_type;
        /// @brief Alias for the input coordinate type template parameter.
        using input_coord_t = InputCoordT;
        /// @brief Alias for the output coordinate type template parameter.
        using output_coord_t = OutputCoordT;
        /// @brief Alias for the core GIS conversion utility struct templated on value_type.
        using converter_t = gis::conversion<value_type>;
        /// @brief Alias for the Helmert parameters struct templated on value_type.
        using helmert_params_t = gis::helmert_params<value_type>;

        /// @brief The number of predefined Helmert transformation parameter sets to use.
        static constexpr size_t num_helmert_params = 4;

    public:
        /// @brief Runs the coordinate conversion process using command-line arguments.
        ///
        /// Parses the command-line arguments to get the input coordinates, then performs
        /// the conversion from `InputCoordT` to `OutputCoordT` using each of the predefined
        /// Helmert parameter sets. Prints each result and the name of the parameter set used.
        /// Finally, calculates and prints the difference between the first result and subsequent results.
        ///
        /// @param argc The number of command-line arguments. Expected to be 3 or 4 (program name + 2 or 3 coordinate values).
        /// @param argv An array of C-style strings containing the command-line arguments. argv[1], argv[2], [argv[3]] should be the
        /// coordinate values.
        /// @param out The output stream to print results and differences to.
        /// @throws std::invalid_argument If the number of arguments is incorrect (not 3 or 4), or if arguments have an invalid numeric
        /// format, or are empty.
        /// @throws std::out_of_range If numeric arguments are outside the representable range for `value_type`.
        /// @throws std::runtime_error If an error occurs during the underlying GIS conversion process (e.g., convergence failure, invalid
        /// intermediate values).
        void operator()(const int argc, const char* const argv[], std::ostream& out) noexcept(false)
        {
            // Parse command line arguments into the input coordinate structure.
            // Throws on parsing errors (invalid count, format, range).
            const auto input_coord = parse_args(argc, argv);

            // Define the Helmert parameter sets to be used from the gis::conversion struct.
            using tf = converter_t::transformation;
            static constexpr std::array<std::reference_wrapper<const helmert_params_t>, num_helmert_params> helmert_params {
                std::cref(tf::pulkovo58_to_wgs84_non_standard),       ///< Non-standard Pulkovo 1942(58) -> WGS84 params.
                std::cref(tf::pulkovo42_to_wgs84_epsg1241),           ///< EPSG:1241 Pulkovo 1942 -> WGS84 params.
                std::cref(tf::pulkovo58_to_wgs84_epsg15861),          ///< EPSG:15861 Pulkovo 1942(58) -> WGS84 params.
                std::cref(tf::dealul_piscului_1970_to_wgs84_epsg1838) ///< EPSG:1838 Dealul Piscului 1970 -> WGS84 params.
            };

            // Array to store the results for each Helmert parameter set.
            std::array<output_coord_t, helmert_params.size()> resulted_coords {};
            auto dest = resulted_coords.begin(); // Iterator for storing results.

            // Create an instance of the converter utility.
            converter_t converter {};

            // Iterate through each Helmert parameter set, perform conversion, and print results.
            out << "Conversion results using different Helmert parameters:\n";
            for (const auto& params_ref: helmert_params)
            {
                const auto& params = params_ref.get(); // Get the actual Helmert parameters object.
                const std::string_view params_name = params.name != nullptr ? params.name : "Unnamed";

                // Perform the conversion using the converter instance, input coordinate, and current parameters.
                // This step can throw std::runtime_error if the conversion fails.
                *dest = converter(input_coord, params);

                // Print the resulting output coordinate and the name of the Helmert parameters used.
                dest->print(out) << " <- [" << params_name << "]" << std::endl;
                ++dest; // Move to the next position in the results array.
            }

            // Calculate and print the differences between the first result and the subsequent results.
            out << "\nDifferences relative to the first result ("
                << (helmert_params[0].get().name ? helmert_params[0].get().name : "Unnamed") << "):\n";
            if constexpr (num_helmert_params > 1) // Only calculate differences if there's more than one result.
            {
                for (size_t i = 1; i < num_helmert_params; ++i)
                {
                    // Use the generic diff function (requires appropriate overload for OutputCoordT).
                    // Assumes diff returns a tuple of differences (dx, dy, dz or dLat, dLon, dAlt).
                    const auto tuple_diff = diff(resulted_coords.front(), resulted_coords[i]);
                    using std::get;
                    using std::setw;
                    // Print differences with labels (assuming d0, d1, d2 correspond to x/y/z or lat/lon/alt diffs)
                    out << "Diff vs [" << setw(38) << (helmert_params[i].get().name ? helmert_params[i].get().name : "Unnamed") << "]: ";
                    out << "d0 = " << std::fixed << std::setprecision(5) << setw(10) << get<0>(tuple_diff)
                        << ", "; // Adjust precision as needed
                    out << "d1 = " << std::fixed << std::setprecision(5) << setw(10) << get<1>(tuple_diff) << ", ";
                    out << "d2 = " << std::fixed << std::setprecision(5) << setw(10) << get<2>(tuple_diff) << std::endl;
                }
            }
            else
            {
                out << "(No differences to calculate, only one parameter set used)\n";
            }
        }

    protected:
        /// @brief Parses a single command-line argument string into a numeric value.
        ///
        /// Attempts to convert the string at the given index in `argv` to type `T` using `std::stod`
        /// (implicitly converting `double` to `T` if necessary). Ensures the entire string is consumed.
        ///
        /// @tparam T The numeric type to parse the argument into (e.g., float, double).
        /// @param index The index of the argument in the `argv` array (e.g., 1 for the first value).
        /// @param argv The array of C-style strings containing the command-line arguments.
        /// @return The parsed numeric value of type T.
        /// @throws std::invalid_argument If the argument string is empty, does not represent a valid number, or contains extra characters
        /// after the number.
        /// @throws std::out_of_range If the parsed number is outside the range representable by `double` (or subsequently `T` if conversion
        /// narrows).
        template <typename T>
        static T parse_arg(const int index, const char* const argv[]) noexcept(false)
        {
            // Ensure index is valid? No, caller (parse_args) ensures argv[index] is accessed within bounds checked by argc.
            size_t processed_chars = 0;
            const std::string arg_str = argv[index]; // Copy to std::string for convenience

            if (arg_str.empty())
                throw std::invalid_argument("Argument at index " + std::to_string(index) + " cannot be empty.");

            // Use std::stod for parsing, as it handles various float formats and provides error checking.
            // Note: This parses as double first. Potential precision loss if T is float.
            // Consider std::stof if T is float, or more advanced parsing if needed.
            try
            {
                const double parsed_double = std::stod(arg_str, &processed_chars);

                // Ensure the entire string was consumed during parsing.
                if (processed_chars != arg_str.length())
                    throw std::invalid_argument("Invalid characters found after number in argument " + std::to_string(index) + ": '" +
                                                arg_str + "'");

                // Convert the parsed double to the target type T. Check for potential narrowing issues if necessary.
                // Static cast is usually sufficient for double -> float/double/long double.
                const T value = static_cast<T>(parsed_double);
                // Optional: Add range check specific to type T if std::stod's double range isn't sufficient.
                return value;
            }
            catch (const std::invalid_argument&)
            {
                // std::stod throws invalid_argument if no conversion could be performed.
                throw std::invalid_argument("Invalid number format in argument " + std::to_string(index) + ": '" + arg_str + "'");
            }
            catch (const std::out_of_range&)
            {
                // std::stod throws out_of_range if the value is out of double's range.
                throw std::out_of_range("Numeric value out of range in argument " + std::to_string(index) + ": '" + arg_str + "'");
            }
        }

        /// @brief Parses command-line arguments into an input coordinate object.
        ///
        /// Expects 2 or 3 numeric arguments after the program name (argc = 3 or 4).
        /// Parses argv[1], argv[2], and optionally argv[3] into the components of `input_coord_t`.
        /// If only 2 arguments are provided (argc=3), the third coordinate component (e.g., altitude/Z) is set to 0.
        ///
        /// @param argc The total number of command-line arguments.
        /// @param argv The array of C-style strings containing the command-line arguments.
        /// @return An `input_coord_t` object initialized with the parsed values.
        /// @throws std::invalid_argument If `argc` is not 3 or 4, or if any argument fails parsing via `parse_arg`.
        /// @throws std::out_of_range If any numeric argument is out of range via `parse_arg`.
        static input_coord_t parse_args(const int argc, const char* const argv[]) noexcept(false)
        {
            // Check if the number of arguments is correct (program name + 2 or 3 values).
            if (argc == 4 || argc == 3) // Allows 3 values (argc=4) or 2 values (argc=3)
            {
                try
                {
                    // Call parse_arg for each required argument.
                    const value_type val1 = parse_arg<value_type>(1, argv);
                    const value_type val2 = parse_arg<value_type>(2, argv);
                    // Parse the third argument only if provided, otherwise default to 0.
                    const value_type val3 = (argc == 4) ? parse_arg<value_type>(3, argv) : static_cast<value_type>(0);

                    // Construct and return the input coordinate object. Assumes InputCoordT has a constructor taking (val1, val2, val3).
                    return input_coord_t(val1, val2, val3);
                }
                catch (const std::invalid_argument& e)
                {
                    // Re-throw parsing errors with potentially more context.
                    throw std::invalid_argument("Error parsing command line arguments: " + std::string(e.what()));
                }
                catch (const std::out_of_range& e)
                {
                    // Re-throw range errors.
                    throw std::out_of_range("Error parsing command line arguments: " + std::string(e.what()));
                }
            }
            else // Incorrect number of arguments.
            {
                throw std::invalid_argument("Incorrect number of arguments provided. Expected 2 or 3 coordinate values (got " +
                                            std::to_string(argc - 1) + "). Usage: <prog> val1 val2 [val3]");
            }
            // Should not be reachable due to the throw above, but prevents compiler warning.
            // return input_coord_t{};
        }

    }; // class coordinate_converter

} // namespace kmx::gis
