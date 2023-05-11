cmake_minimum_required(VERSION 3.25)
project(QuadRetopology)

set(CMAKE_CXX_STANDARD 14)

include_directories(extern extern/vcglib)

include_directories(extern/patterns
        extern/patterns/lp_solve_5.5
        extern/patterns/lp_solve_5.5/shared
        extern/patterns/lp_solve_5.5/bfp
        extern/patterns/lp_solve_5.5/bfp/bfp_LUSOL
        extern/patterns/lp_solve_5.5/bfp/bfp_LUSOL/LUSOL
        extern/patterns/lp_solve_5.5/colamd)


set(DIRS ${DIRS};src)
set(DIRS ${DIRS};extern)
set(DIRS ${DIRS};extern/Eigen)
set(DIRS ${DIRS};extern/patterns)
set(DIRS ${DIRS};extern/vcglib)

file(GLOB_RECURSE openmesh_src extern/patterns/OpenMesh/Core/*.cc)
file(GLOB_RECURSE patchgen_src extern/patterns/patchgen/*.cpp)
file(GLOB lpsolve_src extern/patterns/lp_solve_5.5/*.c)
file(GLOB pattern_src extern/patterns/*.cpp)

file(GLOB_RECURSE file_h *.h  extern/*.h)

foreach(DIR ${DIRS})
    string(REGEX REPLACE "\\+" "\\\\+" DIR ${DIR})
    source_group(${DIR} REGULAR_EXPRESSION ${DIR}/.*)
endforeach()

add_executable(${CMAKE_PROJECT_NAME}
        demo.cpp
        qr_internal.cpp
        extern/vcglib/wrap/ply/plylib.cpp
        extern/patterns/determine_geometry.cpp
        extern/patterns/generate_patch.cpp
        extern/patterns/Patch.cpp
        extern/patterns/lp_solve_5.5/shared/commonlib.c
        extern/patterns/lp_solve_5.5/shared/mmio.c
        extern/patterns/lp_solve_5.5/shared/myblas.c
        extern/patterns/lp_solve_5.5/colamd/colamd.c
        extern/patterns/lp_solve_5.5/bfp/bfp_LUSOL/lp_LUSOL.c
        extern/patterns/lp_solve_5.5/bfp/bfp_LUSOL/LUSOL/lusol.c
        ${openmesh_src}
        ${patchgen_src}
        ${lpsolve_src}
        ${file_h}
        )

target_compile_definitions(${CMAKE_PROJECT_NAME} PUBLIC
        YY_NEVER_INTERACTIVE
        PARSER_LP
        INVERSE_ACTIVE=INVERSE_LUSOL
        RoleIsExternalInvEngine
        )



if(CYGWIN OR MSYS)
    # detect if we use the MinGW compilers on Cygwin - if we do, handle the build
    # as a pure Windows build and make cmake find pure Windows import libraries
    # (.lib)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
            CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpmachine
                OUTPUT_VARIABLE CXX_COMPILER_MACHINE
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(CXX_COMPILER_MACHINE MATCHES "mingw")
            set(WIN32 1)
            add_definitions(-DWIN32)
        endif()
    endif()
endif()


if(MSVC)
    add_definitions(/bigobj -DNOMINMAX -D_CRT_SECURE_NO_DEPRECATE -D_CRT_SECURE_NO_WARNINGS
            -D_SCL_SECURE_NO_DEPRECATE)
    add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/source-charset:utf-8>")
    #    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O2")
    #    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MD")
    #    string(PREPEND "/RTC1" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endif()


if(WIN32)
    add_definitions(-D_USE_MATH_DEFINES)
endif()