import qbs

StaticLibrary {
    Depends { name: "cpp" }
    consoleApplication: true
    cpp.cxxLanguageVersion: "c++23"
    cpp.enableRtti: false
    install: true
    name: "gis-coord-conversion-app-common"
    cpp.includePaths: [
        "inc",
        "inc_dep"
    ]
    files: [
        "inc/kmx/gis/coordinate_converter.hpp",
    ]
}
