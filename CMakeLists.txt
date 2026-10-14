cmake_minimum_required(VERSION 3.12)

project(3rd_party_libraries LANGUAGES CXX C VERSION 1.0.0)

# Create library variables
###############################################################################################
set(3rd_party_libraries_LIBRARY_COMPILE_DEFINITIONS 3rd_party_libraries_compile_definitions)
set(3rd_party_libraries_INSTALL_DIR ${PROJECT_NAME})
set(3rd_party_libraries_INSTALL_BINARY_DIR bin/)
set(3rd_party_libraries_INSTALL_INCLUDE_DIR include/)
set(3rd_party_libraries_INSTALL_LIB_DIR lib/)
set(3rd_party_libraries_INSTALL_ARCHIVE_DIR lib/)
set(3rd_party_libraries_INSTALL_DESTINATION lib/cmake/${PROJECT_NAME})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()
Message(STATUS "3rd_party_libraries Build configuration: " ${CMAKE_BUILD_TYPE})

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${3rd_party_libraries_INSTALL_DIR} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()
Message(STATUS "3rd_party_libraries Library will be installed in: " ${CMAKE_INSTALL_PREFIX})

# External libraries
###########################################################################################
set(ENABLE_TRIANGLE OFF CACHE BOOL "Enable Triangle extension")
set(ENABLE_TETGEN OFF CACHE BOOL "Enable Tetgen extension")
set(ENABLE_METIS OFF CACHE BOOL "Enable Metis extension")
set(ENABLE_VTK OFF CACHE BOOL "Enable VTK extension")
set(ENABLE_MKL OFF CACHE BOOL "Enable MKL extension")
set(ENABLE_VORO OFF CACHE BOOL "Enable VORO extension")
set(ENABLE_SUITESPARSE OFF CACHE BOOL "Enable SuiteSparse extension")
set(ENABLE_PETSC OFF CACHE BOOL "Enable PETSc extension: use MPI to active it")

# Add dependencies
###########################################################################################
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

## GeDiM
if (NOT TARGET GeDiM::GeDiM)
        find_package(GeDiM CONFIG REQUIRED)
        list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES PRIVATE GeDiM::GeDiM)
else()
        get_target_property(GEDIM_LIBRARY_HEADERS GeDiM::GeDiM PUBLIC_HEADER)
        get_target_property(GEDIM_LIBRARY_NAME GeDiM::GeDiM EXPORT_NAME)
        get_target_property(GEDIM_LIBRARY_LINKED_LIBRARIES GeDiM::GeDiM LINK_LIBRARIES)
        get_target_property(GEDIM_LIBRARY_INCLUDE_DIRECTORIES GeDiM::GeDiM INCLUDE_DIRECTORIES)
        get_target_property(GEDIM_LIBRARY_SOURCES GeDiM::GeDiM SOURCES)

        list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES PRIVATE ${GEDIM_LIBRARY_NAME} ${GEDIM_LIBRARY_LINKED_LIBRARIES})
        list(APPEND 3rd_party_libraries_LIBRARY_HEADERS ${GEDIM_LIBRARY_HEADERS})
        list(APPEND 3rd_party_libraries_LIBRARY_INCLUDE ${GEDIM_LIBRARY_INCLUDE_DIRECTORIES})
endif()

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Blas
find_package(BLAS REQUIRED)
list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES PUBLIC BLAS::BLAS)

## Lapack
find_package(LAPACK REQUIRED)
list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES PUBLIC LAPACK::LAPACK)

## MPI
if(${USE_MPI})
    find_package(MPI REQUIRED)
    list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES PUBLIC MPI::MPI_CXX)
endif()

## MKL
if(${ENABLE_MKL})
    set(MKL_INTERFACE lp64) # Eigen requires lp64 interface
    find_package(MKL CONFIG REQUIRED)
    list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES MKL::MKL)
endif()

## triangle
if(${ENABLE_TRIANGLE})
    find_package(triangle CONFIG REQUIRED)
    list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES PUBLIC triangle::triangle)
endif()

## tetgen
if(${ENABLE_TETGEN})
    find_package(tetgen CONFIG REQUIRED)
    list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES PUBLIC tetgen::tetgen)
endif()

## VORO
if(${ENABLE_VORO})
    find_package(VORO CONFIG REQUIRED)
    list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES PUBLIC VORO::voro++)
endif()

## VTK
if(${ENABLE_VTK})
    find_package(VTK CONFIG REQUIRED)
    list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES PUBLIC ${VTK_LIBRARIES})
endif()

## METIS
if(${ENABLE_METIS})
    find_package(metis CONFIG REQUIRED)
    list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES METIS::metis)
endif()

## SuiteSparse
if(${ENABLE_SUITESPARSE})
    find_package(SuiteSparse CONFIG REQUIRED)
    message("Found SuiteSparse: ${SUITESPARSE_LIBRARIES}")
    list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES ${SUITESPARSE_LIBRARIES})
    list(APPEND 3rd_party_libraries_LIBRARY_INCLUDE ${SUITESPARSE_INCLUDE_DIRS})
endif()

## PETSC
if(${ENABLE_PETSC})
    # PkgConfig
    find_package(PkgConfig)

    # PETSc
    if (PKG_CONFIG_FOUND)
        pkg_check_modules(PETSC PETSc)
    endif()

    if (PETSC_FOUND)
        message("Found PETSc: ${PETSC_LINK_LIBRARIES}")
        list(APPEND COMPILE_OPTIONS ${PETSC_CFLAGS})
        list(APPEND 3rd_party_libraries_LIBRARY_INCLUDE ${PETSC_INCLUDE_DIRS})
        list(APPEND 3rd_party_libraries_LIBRARY_LINKED_LIBRARIES ${PETSC_LINK_LIBRARIES})
        set(LINK_FLAGS "${LINK_FLAGS} ${PETSC_LDFLAGS}")
        set(CMAKE_REQUIRED_FLAGS ${PETSC_CFLAGS})
        set(CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDE_DIRS}")
    endif()
endif()

# Insert Sources
###############################################################################################
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
list(APPEND 3rd_party_libraries_LIBRARY_SOURCES ${3rd_party_libraries_src_sources})
list(APPEND 3rd_party_libraries_LIBRARY_HEADERS ${3rd_party_libraries_src_headers})
list(APPEND 3rd_party_libraries_LIBRARY_INCLUDE ${3rd_party_libraries_src_includes})

# Add variables for clang anc cppcheck
###############################################################################################
list(APPEND CLANGFORMAT_FILES ${3rd_party_libraries_src_sources})
list(APPEND CLANGFORMAT_FILES ${3rd_party_libraries_src_headers})
set(CLANGFORMAT_FILES ${CLANGFORMAT_FILES} PARENT_SCOPE)

list(APPEND CPPCHECK_FILES ${3rd_party_libraries_src_sources})
list(APPEND CPPCHECK_FILES ${3rd_party_libraries_src_headers})
set(CPPCHECK_FILES ${CPPCHECK_FILES} PARENT_SCOPE)

# Create public headers
###############################################################################################
set(3rd_party_libraries_LIBRARY_ADDITIONAL_INCLUDE_PATH ${CMAKE_CURRENT_BINARY_DIR}/3rd_party_libraries_additional_include)

foreach(header ${3rd_party_libraries_LIBRARY_HEADERS})
    get_filename_component(headerDirectory ${header} DIRECTORY)
    get_filename_component(headerFileName ${header} NAME_WE)
    get_filename_component(headerExtension ${header} EXT)

    if(${headerExtension} MATCHES ".in")
        set(headerExtension ".hpp")
        set(headerDirectory ${3rd_party_libraries_LIBRARY_ADDITIONAL_INCLUDE_PATH})

        configure_file(${header} ${headerDirectory}/${headerFileName}${headerExtension} @ONLY)
    elseif(NOT ${headerExtension} MATCHES ".hpp" AND NOT ${headerExtension} MATCHES ".h")
        message(SEND_ERROR "Header extension not recognized for file ${headerFileName}${headerExtension}")
    endif()

    list(APPEND 3rd_party_libraries_LIBRARY_PUBLIC_HEADERS ${headerDirectory}/${headerFileName}${headerExtension})
endforeach()

#foreach(header ${3rd_party_libraries_LIBRARY_SOURCES})
#    get_filename_component(headerDirectory ${header} DIRECTORY)
#    get_filename_component(headerFileName ${header} NAME_WE)
#    get_filename_component(headerExtension ${header} EXT)

#    if(${headerFileName} MATCHES "GeometryUtilities")
#        MESSAGE(WARNING "Add option to " ${header})
#        set_source_files_properties(${header} PROPERTIES COMPILE_OPTIONS -Weverything)
#    endif()
#
#endforeach()

list(APPEND 3rd_party_libraries_LIBRARY_INCLUDE ${3rd_party_libraries_LIBRARY_ADDITIONAL_INCLUDE_PATH})

# Create library
###############################################################################################
add_library(${PROJECT_NAME} STATIC ${3rd_party_libraries_LIBRARY_SOURCES} ${3rd_party_libraries_LIBRARY_HEADERS})
add_library(3rd_party_libraries::3rd_party_libraries ALIAS ${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} ${3rd_party_libraries_LIBRARY_LINKED_LIBRARIES})

target_include_directories(${PROJECT_NAME} PRIVATE ${3rd_party_libraries_LIBRARY_INCLUDE})
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:${3rd_party_libraries_INSTALL_INCLUDE_DIR}>)
target_compile_definitions(${PROJECT_NAME} PUBLIC ${${3rd_party_libraries_LIBRARY_COMPILE_DEFINITIONS}})

target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

# Create install target
###############################################################################################
set_target_properties(${PROJECT_NAME} PROPERTIES 
    PUBLIC_HEADER "${3rd_party_libraries_LIBRARY_PUBLIC_HEADERS}"
    EXPORT_NAME 3rd_party_libraries
    LINKER_LANGUAGE CXX
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF)

install(
    TARGETS ${PROJECT_NAME}
    EXPORT 3rd_party_librariesTargets
    RUNTIME DESTINATION ${3rd_party_libraries_INSTALL_BINARY_DIR} # runtime executable
    LIBRARY DESTINATION ${3rd_party_libraries_INSTALL_LIB_DIR} # dynamic libraries
    ARCHIVE DESTINATION ${3rd_party_libraries_INSTALL_ARCHIVE_DIR} # static libraries
    PUBLIC_HEADER DESTINATION ${3rd_party_libraries_INSTALL_INCLUDE_DIR} # headers
)


include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/3rd_party_librariesConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/3rd_party_librariesConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/${3rd_party_libraries_INSTALL_DESTINATION}
    PATH_VARS 3rd_party_libraries_INSTALL_INCLUDE_DIR
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/3rd_party_librariesConfigVersion.cmake
    VERSION ${3rd_party_libraries_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(
    TARGETS ${PROJECT_NAME}
    NAMESPACE 3rd_party_libraries::
    FILE ${CMAKE_CURRENT_BINARY_DIR}/3rd_party_librariesTargets.cmake
)

install(
    EXPORT 3rd_party_librariesTargets
    NAMESPACE 3rd_party_libraries::
    DESTINATION ${CMAKE_INSTALL_PREFIX}/${3rd_party_libraries_INSTALL_DESTINATION}
)

install(
    FILES
    ${CMAKE_CURRENT_BINARY_DIR}/3rd_party_librariesConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/3rd_party_librariesConfigVersion.cmake
    DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${3rd_party_libraries_INSTALL_DESTINATION}
)
