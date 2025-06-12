add_library(disp
    secfunc.cc
    disp.cc
    toms748.cc
    utils.cc
    swegn96.cc
)

target_link_libraries(disp
    PRIVATE
    fmt::fmt-header-only
    segn96_fortran
)

add_executable(secfunc
    main_secfunc.cc
)

target_link_libraries(secfunc
    PRIVATE
    disp
    HighFive
    fmt::fmt-header-only
    toml11
)

add_executable(forward
    main_forward.cc
)

target_link_libraries(forward
    PRIVATE
    disp
    fmt::fmt-header-only
    toml11
)
