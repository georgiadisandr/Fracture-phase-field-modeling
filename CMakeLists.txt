cmake_minimum_required(VERSION 3.16)
project(PhaseFieldFEM1 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Build type (Release / Debug / RelWithDebInfo / MinSizeRel)" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ── Eigen3 ──────────────────────────────────────────────────────────────────
if(WIN32)
    set(EIGEN_DIR "" CACHE PATH "Path to Eigen headers directory (Windows only)")
    if(NOT EIGEN_DIR)
        message(FATAL_ERROR "Set EIGEN_DIR to your Eigen headers path, e.g. -DEIGEN_DIR=C:/path/to/eigen")
    endif()
else()
    find_package(Eigen3 REQUIRED
        HINTS /usr/share/eigen3/cmake)
endif()

# ── Gmsh ────────────────────────────────────────────────────────────────────
if(WIN32)
    set(GMSH_INCLUDE_DIR "" CACHE PATH "Path to Gmsh include dir (Windows only)")
    set(GMSH_LIB_DIR     "" CACHE PATH "Path to Gmsh lib dir (Windows only)")
    if(NOT GMSH_INCLUDE_DIR OR NOT GMSH_LIB_DIR)
        message(FATAL_ERROR "Set GMSH_INCLUDE_DIR and GMSH_LIB_DIR to your Gmsh SDK paths")
    endif()
    set(GMSH_LIBRARY "${GMSH_LIB_DIR}/gmsh-4.15.dll")
else()
    find_path(GMSH_INCLUDE_DIR NAMES gmsh.h REQUIRED)
    find_library(GMSH_LIBRARY  NAMES gmsh   REQUIRED)
endif()

message(STATUS "Gmsh include: ${GMSH_INCLUDE_DIR}")
message(STATUS "Gmsh library: ${GMSH_LIBRARY}")

add_executable(phasefield_1
    main.cpp
    MeshGeneration.cpp
    Gauss.cpp
    Jacob2.cpp
    ShapeFunc.cpp
    GmshReader.cpp
    ConstitutiveModel.cpp
    ElementUtils.cpp
    StiffnessPFM.cpp
    ResidualPFM.cpp
    History.cpp
    BoundaryConditions.cpp
    NeumannBC.cpp
    Solver.cpp
    Output.cpp
    Config.cpp
)

# The project root is on the include path so the vendored single-header
# tomlplusplus.hpp (sitting next to the .cpp files) is found by Config.cpp.
target_include_directories(phasefield_1 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${GMSH_INCLUDE_DIR})

if(WIN32)
    target_include_directories(phasefield_1 PRIVATE ${EIGEN_DIR})
    target_link_libraries(phasefield_1 PRIVATE "${GMSH_LIBRARY}")
else()
    target_link_libraries(phasefield_1 PRIVATE Eigen3::Eigen "${GMSH_LIBRARY}")
endif()

# On Windows, copy the gmsh DLL next to the executable after build so the .exe
# can find it without manual PATH juggling.
if(WIN32)
    find_file(GMSH_DLL NAMES gmsh-4.15.dll gmsh.dll
        HINTS "${GMSH_LIB_DIR}"
        PATH_SUFFIXES bin lib)
    if(GMSH_DLL)
        add_custom_command(TARGET phasefield_1 POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${GMSH_DLL}" "$<TARGET_FILE_DIR:phasefield_1>")
    else()
        message(WARNING "gmsh DLL not found; phasefield_1.exe will fail to start at runtime")
    endif()
endif()
