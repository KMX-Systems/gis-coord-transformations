import qbs

CppApplication {
    Depends
    {
        name: 'gis-coord-conversion-lib'
    }

    name: "gis-coord-conversion-app-stereo70-to-wgs84"
    consoleApplication: true
    cpp.debugInformation: true

    files: [
        "src/kmx/gis/app/stereo70-to-wgs84.cpp",
    ]
    cpp.cxxLanguageVersion: "c++23"
    cpp.enableRtti: false
    cpp.includePaths: [
        "inc",
        "inc_dep"
    ]
    cpp.systemIncludePaths: [
        "/usr/local/include",
        "/usr/include"
    ]
    cpp.staticLibraries: [
        "/usr/local/lib/libCatch2Main.a",
        "/usr/local/lib/libCatch2.a"
    ]
}
