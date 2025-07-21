add_library(timer
    timer.cc
)

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
    HighFive
    fmt::fmt-header-only
    toml11
)

add_executable(inversion
    main_inversion.cc
    problem.cc
    model.cc
)

target_link_libraries(inversion
    PRIVATE
    disp
    fmt::fmt-header-only
    toml11
    HighFive
    lbfgspp
    timer
)

if(ENABLE_DEBUG)
    target_compile_definitions(inversion PRIVATE DEBUG)
endif()

if (USE_OPENMP)
target_link_libraries(inversion
    PRIVATE
    OpenMP::OpenMP_CXX
)
endif()