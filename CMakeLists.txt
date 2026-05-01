##################################################################################
#                                                                                #
# CMake file to be used to install KDSource distribution.                        #
#                                                                                #
# See README.md file for installation instructions.                              #
#                                                                                #
# Written 2021 by O. I. Abbate.                                                  #
# Major rewrite 2026 by T. Kittelmann.                                           #
#                                                                                #
##################################################################################

cmake_minimum_required(VERSION 3.26.1...3.31)

file(READ "${CMAKE_SOURCE_DIR}/VERSION" PROJECT_VERSION)
string(STRIP "${PROJECT_VERSION}" PROJECT_VERSION)

project(KDSource VERSION "${PROJECT_VERSION}" LANGUAGES C )

option( KDS_EMBED_LIBXML2 "Download and embed libxml2 statically" OFF)
option( KDS_STRICT "Stricter compilation " OFF)

###################################################################
#Detect sdist generation mode:
###################################################################

#NB: When generating the sdist, all we want to do is essentially download the
#    libxml2 sources for inclusion.

set( KDS_GENERATING_SDIST OFF )
if ( DEFINED SKBUILD_PROJECT_NAME AND "${SKBUILD_STATE}" STREQUAL "sdist" )
  set( KDS_GENERATING_SDIST ON )
endif()

###################################################################
#Find dependencies:
###################################################################

if ( NOT KDS_EMBED_LIBXML2 )
  if ( NOT KDS_GENERATING_SDIST )
    find_package(LibXml2 REQUIRED)
  endif()
else()
  # Configure FetchContent for libxml2 and ensure it is a static lib. If running
  # from an sdist, it might already have been fetched into third_party/libxml2.
  set( KDS_LIBXML2SRC "${CMAKE_CURRENT_BINARY_DIR}/kds_third_party_libxml2_src" )
  if ( DEFINED SKBUILD_PROJECT_NAME )
    set( KDS_LIBXML2SRC "${CMAKE_CURRENT_SOURCE_DIR}/third_party/libxml2" )
  endif()
  if ( NOT EXISTS "${KDS_LIBXML2SRC}/CMakeLists.txt" )
    include(FetchContent)
    FetchContent_Populate(
      libxml2
      #Use our own patched version, for a leaner build and to avoid installing
      #system libxml2 binaries into the wheel.
      URL https://github.com/KDSource/libxml2/archive/refs/tags/v2.15.3_kdsource_patch1.tar.gz
      URL_HASH SHA256=f422b40be3b99356bf9d397c9657a1a05ce223dc6c898af8857e827297ff036f
      SOURCE_DIR "${KDS_LIBXML2SRC}"
    )
  endif()
  if ( NOT KDS_GENERATING_SDIST )
    set(KDS_OLD_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
    set(BUILD_SHARED_LIBS OFF)
    add_subdirectory("${KDS_LIBXML2SRC}")
    set(BUILD_SHARED_LIBS ${KDS_OLD_BUILD_SHARED_LIBS})
  endif()
endif()

if ( KDS_GENERATING_SDIST )
  #sdist generation does not need anything else
  return()
endif()

if( NOT DEFINED "MCPL_DIR" )
  #Need to invoke "mcpl-config --show cmakedir" if we want to be able to
  #work with mcpl-core installed via python wheels:
  execute_process(
    COMMAND mcpl-config --show cmakedir
    OUTPUT_VARIABLE "MCPL_DIR" OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()
find_package(MCPL REQUIRED)

###################################################################
#Prepare autogen dir for processed files (possibly in-source if building via
#python tools)
###################################################################

set( kdsc_autogendir "${CMAKE_CURRENT_BINARY_DIR}/kdsc_autogen" )

#Always cleanup autogendir
if ( IS_DIRECTORY "${kdsc_autogendir}" )
  message( STATUS "Cleaning out existing ${kdsc_autogendir}" )
  file( REMOVE_RECURSE "${kdsc_autogendir}" )
endif()
file( MAKE_DIRECTORY "${kdsc_autogendir}" )

###################################################################
#Find source files and prepare in autogen dir, including
#processed .in files:
###################################################################

#On windows DLL's must be in %PATH% so we put them in bin (same as seems to be
#the way most packages do it on conda-forge). For wheels we will place the
#mcpl.dll into the scripts folder, since that seems to be the only way to get
#the dll onto the PATH. The following variable decides if ultimately placing
#.dll in wheel scripts (boolean value should here be "1" or "" due to
#configure_file usage):
set( kds_skbuild_shlib_in_wheel_scripts "" )
if ( WIN32 AND DEFINED SKBUILD_PROJECT_NAME )
  set( kds_skbuild_shlib_in_wheel_scripts "1" )
endif()

set( KDS_PYMODNAME "kdsource" )

file(GLOB KDS_C CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/clib/src/*.c")
file(GLOB KDS_H CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/clib/include/kdsource/*.h")
file(GLOB KDS_HIN CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/clib/include/kdsource/*.h.in")
file(GLOB KDS_PY CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/python/${KDS_PYMODNAME}/*.py")
file(GLOB KDS_PYIN CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/python/${KDS_PYMODNAME}/*.py.in")

set( t1 "1000000 * ${PROJECT_VERSION_MAJOR}" )
set( t2 "1000 * ${PROJECT_VERSION_MINOR}" )
set( t3 "${PROJECT_VERSION_PATCH}" )
math(EXPR "${PROJECT_NAME}_version_int" "(${t1})+(${t2})+(${t3})" )

file( COPY ${KDS_PY} DESTINATION "${kdsc_autogendir}/python/kdsource/" )
foreach( infile ${KDS_PYIN} )
  get_filename_component(infilebn "${infile}" NAME_WE)
  configure_file(
    "${infile}"
    "${kdsc_autogendir}/python/kdsource/${infilebn}.py"
    @ONLY
  )
endforeach()

file( COPY ${KDS_C} DESTINATION "${kdsc_autogendir}/clib/src" )
file( COPY ${KDS_H} DESTINATION "${kdsc_autogendir}/clib/include/kdsource" )
foreach( infile ${KDS_HIN} )
  get_filename_component( infilebn "${infile}" NAME_WE )
  configure_file(
    "${infile}"
    "${kdsc_autogendir}/clib/include/kdsource/${infilebn}.h"
    @ONLY
  )
endforeach()

###################################################################
# Installation destinations:
###################################################################

include(GNUInstallDirs)
if ( DEFINED SKBUILD_PROJECT_NAME )
  message( STATUS "scikit-build mode detected: clib will be placed in python module." )
  set( KDS_INSTALL_LIBDIR "${KDS_PYMODNAME}/data/lib" )
  set( KDS_INSTALL_INCLUDEDIR "${KDS_PYMODNAME}/data/include" )
else()
  set( KDS_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}" )
  set( KDS_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}" )
endif()

if ( kds_skbuild_shlib_in_wheel_scripts )
  set( KDS_INSTALL_SHLIBDIR "${SKBUILD_SCRIPTS_DIR}" )
else()
  #sane standard world
  set( KDS_INSTALL_SHLIBDIR "${KDS_INSTALL_LIBDIR}" )
endif()

###################################################################
# Build library (entirely from autogendir)
###################################################################

file(GLOB KDS_SOURCES CONFIGURE_DEPENDS "${kdsc_autogendir}/clib/src/*.c")
file(GLOB KDS_HEADERS CONFIGURE_DEPENDS "${kdsc_autogendir}/clib/include/kdsource/*.h")

add_library(kdsource SHARED ${KDS_SOURCES} )

set_target_properties(
  kdsource
  PROPERTIES C_STANDARD 99
  C_STANDARD_REQUIRED ON C_EXTENSIONS OFF
)

target_include_directories(
  kdsource
  PRIVATE "${kdsc_autogendir}/clib/src"
  PUBLIC
  $<BUILD_INTERFACE:${kdsc_autogendir}/clib/include>
  $<INSTALL_INTERFACE:${KDS_INCLUDEDIR}>
)

target_link_libraries(
  kdsource
  PUBLIC LibXml2::LibXml2 MCPL::MCPL
)

if ( WIN32 )
  set_target_properties( kdsource PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE )
endif()

install(
  FILES ${KDS_HEADERS}
  DESTINATION "${KDS_INSTALL_INCLUDEDIR}/kdsource/"
)

if ( KDS_STRICT )
  if ( MSVC )
    #/disable warning 4996 for now, since it will require rewriting the sscanf
    #usage.
    target_compile_options( kdsource PRIVATE /WX /W4 /wd4996 )
  else()
    target_compile_options(
      kdsource PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wformat=2
      -Wconversion -Wsign-conversion -Wnull-dereference -Werror
    )
  endif()
endif()

install(
  TARGETS kdsource
  RUNTIME DESTINATION "${KDS_INSTALL_SHLIBDIR}"
  ARCHIVE DESTINATION "${KDS_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${KDS_INSTALL_LIBDIR}"
)


###################################################################
# Install python modules (only if building a wheel)
###################################################################

if ( DEFINED SKBUILD_PROJECT_NAME )
  file( GLOB KDS_PYWHL CONFIGURE_DEPENDS "${kdsc_autogendir}/python/${KDS_PYMODNAME}/*.py" )
  install( FILES ${KDS_PYWHL} DESTINATION "${SKBUILD_PLATLIB_DIR}/${KDS_PYMODNAME}" )
endif()
