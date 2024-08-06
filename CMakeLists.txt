cmake_minimum_required(VERSION 3.10)
project(birp VERSION 1.0)
add_executable(birp main.cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add_subdirectory("${PROJECT_BINARY_DIR}../Classes")
add_library(Classes Classes/Camera.cpp Classes/DataCube.cpp Classes/Volume.cpp Helper.cpp)
add_subdirectory(Classes)
# let's add the elefits library so we can create .FITS files
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
find_package(EleFits REQUIRED)
target_link_libraries(birp
    PUBLIC Classes
    PRIVATE EleFits
)

configure_file(mainconfig.h.in mainconfig.h)
target_include_directories(birp
    PUBLIC "${PROJECT_BINARY_DIR}"
    PUBLIC "${PROJECT_BINARY_DIR}/Classes"
)

