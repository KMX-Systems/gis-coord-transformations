import qbs 1.0

Project {
    references: [
        "application/common/lib.qbs",
        "application/wgs84-to-stereo70/app.qbs",
        "application/stereo70-to-wgs84/app.qbs",
        "library/lib.qbs",
        "library-test/unit-test.qbs"
    ]
}

