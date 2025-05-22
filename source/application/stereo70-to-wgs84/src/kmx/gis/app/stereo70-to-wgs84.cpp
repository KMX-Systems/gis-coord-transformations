#include "kmx/gis/coordinate/custom_converter.hpp"
#include <exception>

/// @brief Main function for the Stereo70 to WGS84 converter.
/// Parses command line arguments, creates the application object, and runs it.
/// Catches standard exceptions and prints usage information on argument errors.
/// @param[in] argc Argument count passed from the operating system.
/// @param[in] argv Argument vector passed from the operating system.
/// @return 0 on successful execution, 1 on error (e.g., invalid arguments, conversion failure).
int main(int argc, char* argv[])
{
    using T = double;
    using namespace kmx::gis;
    using app_t = coordinate::custom_converter<stereo70::coordinate<T>, wgs84::coordinate<T>>;

    app_t app {};

    try
    {
        app(argc, argv, std::cout);
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    catch (const std::out_of_range& e)
    {
        std::cerr << "Error: Input value out of range. " << e.what() << '\n';
        return 2;
    }
    catch (const std::exception& e)
    {
        std::cerr << "An unexpected error occurred: " << e.what() << '\n';
        return 3;
    }

    return 0;
}
