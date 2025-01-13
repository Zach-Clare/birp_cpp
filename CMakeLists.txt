cmake_minimum_required(VERSION 3.10)
project(birp VERSION 1.0) # create the project handle

# taken from https://stackoverflow.com/questions/41361631/optimize-in-cmake-by-default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_STANDARD 23) 

# let's add the elefits library so we can create .FITS files
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
find_package(EleFits REQUIRED)

# add files to be used as part of the executable
add_executable(birp
    ${CMAKE_SOURCE_DIR}/src/main.cpp
    ${CMAKE_SOURCE_DIR}/src/Camera.cpp
    ${CMAKE_SOURCE_DIR}/src/DataCube.cpp
    ${CMAKE_SOURCE_DIR}/src/Helper.cpp
    ${CMAKE_SOURCE_DIR}/src/Cmem.cpp
)

# set where available includes can come from
target_include_directories(birp
    PRIVATE ${CMAKE_SOURCE_DIR}/src # holds .cpp source files
    PRIVATE ${CMAKE_SOURCE_DIR}/include # holds the .h header files for the classes
)

# set(CMAKE_CXX_STANDARD 11) # use version 11
# set(CMAKE_CXX_STANDARD_REQUIRED True) # force
# set(CMAKE_CXX_EXTENSIONS OFF)
# target_compile_features(birp PRIVATE cxx_std_13)


target_link_libraries(birp
    PRIVATE EleFits EleCfitsioWrapper cfitsio EleFitsData # need ffppr
    # PRIVATE Classes
)

# add_subdirectory("${PROJECT_BINARY_DIR}../Classes")
# add_library(Classes Classes/Camera.cpp Classes/DataCube.cpp Classes/Volume.cpp Helper.cpp)
# add_subdirectory(Classes)

# configure_file(mainconfig.h.in mainconfig.h)
# target_include_directories(birp
#     PUBLIC "${PROJECT_BINARY_DIR}"
#     # PUBLIC "${PROJECT_BINARY_DIR}/Classes"
# )

# ################# Testing framework - can't build due to DataCube linking error?
# set(GOOGLETEST_VERSION 1.15.2)

# # set up test framework
# enable_testing()
# include(CTest)
# find_package(GTest REQUIRED)
# include(GoogleTest)

# # add_subdirectory(googletest)
# # include_directories(googletest/include) # this is so we can #include <gtest/gtest.h>
# add_executable(tests
#     ${CMAKE_SOURCE_DIR}/tests/Camera_test.cpp
#     ${CMAKE_SOURCE_DIR}/tests/CMEM_test.cpp
#     ${CMAKE_SOURCE_DIR}/tests/DataCube_test.cpp
#     ${CMAKE_SOURCE_DIR}/tests/Space_test.cpp
#     # ${CMAKE_SOURCE_DIR}/src/main.cpp
#     ${CMAKE_SOURCE_DIR}/src/Camera.cpp
#     ${CMAKE_SOURCE_DIR}/src/DataCube.cpp
#     ${CMAKE_SOURCE_DIR}/src/Helper.cpp
#     ${CMAKE_SOURCE_DIR}/src/Cmem.cpp
# )
# target_link_libraries(tests GTest::gtest GTest::Main EleFits EleCfitsioWrapper cfitsio EleFitsData)

# # set where available includes can come from
# target_include_directories(tests
#     PRIVATE ${CMAKE_SOURCE_DIR}/src # holds .cpp source files
#     PRIVATE ${CMAKE_SOURCE_DIR}/include # holds the .h header files for the classes
# )


# # include(CTest)
# # include(GoogleTest)
# gtest_discover_tests(tests DISCOVERY_TIMEOUT 600)
# # add_test(hellotest hellotest) # This is the crucial part. Without this, CTest has no idea the GTests exist