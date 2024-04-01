cmake_minimum_required(VERSION 3.12)
project(wgp)

file(GLOB_RECURSE STD_INCLUDE ${PROJECT_SOURCE_DIR}/include/wstd/*.h)
file(GLOB_RECURSE STD_SRC ${PROJECT_SOURCE_DIR}/src/wstd/*.h ${PROJECT_SOURCE_DIR}/src/wstd/*.cpp)
file(GLOB_RECURSE GEO_INCLUDE ${PROJECT_SOURCE_DIR}/include/wgeo/*.h)
file(GLOB_RECURSE GEO_SRC ${PROJECT_SOURCE_DIR}/src/wgeo/*.h ${PROJECT_SOURCE_DIR}/src/wgeo/*.cpp)
file(GLOB_RECURSE SCENE_INCLUDE ${PROJECT_SOURCE_DIR}/include/wscene/*.h)
file(GLOB_RECURSE SCENE_SRC ${PROJECT_SOURCE_DIR}/src/wscene/*.h ${PROJECT_SOURCE_DIR}/src/wscene/*.cpp)

if(${CMAKE_SYSTEM_NAME} MATCHES "windows")
    if (BUILD_SHARED_LIBS)
        set(OUT_PATH ${PROJECT_SOURCE_DIR}/build/shared)
    else()
        set(OUT_PATH ${PROJECT_SOURCE_DIR}/build/static)
    endif()
    if (CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(BIN_DIR ${OUT_PATH}/windows/bin/x64)
        set(LIB_DIR ${OUT_PATH}/windows/lib/x64)
    else()
        set(BIN_DIR ${OUT_PATH}/build/windows/bin/win32)
        set(LIB_DIR ${OUT_PATH}/build/windows/lib/win32)
    endif()

    #wstd
    add_library(wstd ${PROJECT_SOURCE_DIR}/include/wbase.h 
        ${STD_INCLUDE} 
        ${STD_SRC}
    )
    set_target_properties(wstd PROPERTIES 
        ARCHIVE_OUTPUT_DIRECTORY "${LIB_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${BIN_DIR}" 
        RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
        OUTPUT_NAME "wstd"
        PREFIX "")
    target_include_directories(wstd PRIVATE ${PROJECT_SOURCE_DIR}/include)
    #wgeo
    add_library(wgeo ${PROJECT_SOURCE_DIR}/include/wbase.h 
        ${GEO_INCLUDE} 
        ${GEO_SRC}
    )
    set_target_properties(wgeo PROPERTIES 
        ARCHIVE_OUTPUT_DIRECTORY "${LIB_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${BIN_DIR}" 
        RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
        OUTPUT_NAME "wgeo"
        PREFIX "")
    target_include_directories(wgeo PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_link_directories(wgeo PRIVATE ${LIB_DIR})
    target_link_libraries(wgeo PRIVATE wstd)
    #wscene
    add_library(wscene ${PROJECT_SOURCE_DIR}/include/wbase.h 
        ${SCENE_INCLUDE} 
        ${SCENE_SRC}
    )
    set_target_properties(wscene PROPERTIES 
        ARCHIVE_OUTPUT_DIRECTORY "${LIB_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${BIN_DIR}" 
        RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
        OUTPUT_NAME "wscene"
        PREFIX "")
    target_include_directories(wscene PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_link_directories(wscene PRIVATE ${LIB_DIR})
    target_link_libraries(wscene PRIVATE wstd wgeo)
    add_definitions(-DWGP_EXPORTS)
    if (BUILD_SHARED_LIBS)
    else()
        add_definitions(-DWGP_STATIC)
    endif()
endif()

