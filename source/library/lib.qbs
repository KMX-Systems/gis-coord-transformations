import qbs

StaticLibrary {
    Depends { name: "cpp" }
    consoleApplication: true
    cpp.cxxLanguageVersion: "c++23"
    cpp.enableRtti: false
    install: true
    name: "gis-coord-conversion-lib"
    files: "inc/kmx/gis.hpp"
}
