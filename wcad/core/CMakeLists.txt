cmake_minimum_required(VERSION 3.12)
project(wcad)

file(GLOB_RECURSE WCAD_INCLUDE ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB_RECURSE WCAD_SRC ${PROJECT_SOURCE_DIR}/src/*.h ${PROJECT_SOURCE_DIR}/src/*.cpp)

if(${CMAKE_SYSTEM_NAME} MATCHES "windows")
    if (BUILD_SHARED_LIBS)
        set(OUT_PATH ${PROJECT_SOURCE_DIR}/build/shared)
    else()
        set(OUT_PATH ${PROJECT_SOURCE_DIR}/build/static)
    endif()
    set(GLEW_DIR ${PROJECT_SOURCE_DIR}/../../wgp/thirdparty/glew)
    if (CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(BIN_DIR ${OUT_PATH}/windows/bin/x64)
        set(LIB_DIR ${OUT_PATH}/windows/lib/x64)
        set(WGP_LIB_DIR ${PROJECT_SOURCE_DIR}/../../wgp/build/static/windows/lib/x64)
        set(GLEW_LIB_PATH ${GLEW_DIR}/lib/Release/x64)
    else()
        set(BIN_DIR ${OUT_PATH}/build/windows/bin/win32)
        set(LIB_DIR ${OUT_PATH}/build/windows/lib/win32)
        set(WGP_LIB_DIR ${PROJECT_SOURCE_DIR}/../../wgp/build/static/windows/lib/win32)
        set(GLEW_LIB_PATH ${GLEW_DIR}/lib/Release/Win32)
    endif()

    #wcad
    add_library(wcad 
        ${WCAD_INCLUDE} 
        ${WCAD_SRC}
    )
    set_target_properties(wcad PROPERTIES 
        ARCHIVE_OUTPUT_DIRECTORY "${LIB_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${BIN_DIR}" 
        RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
        OUTPUT_NAME "wcad"
        PREFIX "")
    target_include_directories(wcad PRIVATE ${PROJECT_SOURCE_DIR}/include ${GLEW_DIR}/include ${PROJECT_SOURCE_DIR}/../../wgp/include)
    target_link_directories(wcad PRIVATE ${WGP_LIB_DIR} ${GLEW_LIB_PATH})
    target_link_libraries(wcad PRIVATE wstd wgeo wscene glew32 opengl32)
    add_definitions(-DWCAD_EXPORTS)
    if (BUILD_SHARED_LIBS)
    else()
        add_definitions(-DWGP_STATIC)
        add_definitions(-DWCAD_STATIC)
    endif()
endif()

