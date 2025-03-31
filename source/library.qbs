import qbs

StaticLibrary {
    Depends { name: "cpp" }
    consoleApplication: true
    cpp.cxxLanguageVersion: "c++20"
    //cpp.cxxFlags: "-gdwarf-4"
    cpp.enableRtti: false
    install: true
    name: "gis-coord-conversion"
    files: "inc/kmx/gis.hpp"
}
