cmake_minimum_required(VERSION 3.10)
project(birp VERSION 1.0) # create the project handle

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

