import qbs

StaticLibrary {
    Depends { name: "cpp" }
    consoleApplication: true
    cpp.cxxLanguageVersion: "c++23"
    cpp.enableRtti: false
    install: true
    name: "gis-coord-conversion-lib"
    cpp.includePaths: [
        "inc",
        "inc_dep"
    ]
    files: [
        "inc/kmx/gis.hpp",
        "inc/kmx/gis/constants.hpp",
        "inc/kmx/gis/coordinate/base.hpp",
        "inc/kmx/gis/coordinate/conversion.hpp",
        "inc/kmx/gis/coordinate/stereo70.hpp",
        "inc/kmx/gis/ellipsoid.hpp",
        "inc/kmx/gis/uniform_grid_mapper.hpp",
    ]
}
