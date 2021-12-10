
##################################################################################
#                                                                                #
# CMake file to be used to install KDSource distribution.                        #
#                                                                                #
# See README.md file for installation instructions.                              #
#                                                                                #
# Written 2021 by O. I. Abbate.                                                  #
#                                                                                #
##################################################################################

cmake_minimum_required(VERSION 3.0.0)

project(KDSource VERSION 0.1.0)
configure_file(KDSourceConfig.h.in KDSourceConfig.h)

enable_testing()

set(INSTALL_PY OFF CACHE STRING "Whether to also install kdsource and mcpl python files.")


set(CMAKE_C_STANDARD 99)

set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")

find_package(LibXml2 REQUIRED)

set(BUILD_WITHG4 OFF CACHE STRING "Whether to build Geant4 plugins if Geant4 is available.")
add_subdirectory(mcpl)

add_library(kdsource SHARED "${SRC}/kdsource/kdsource.c"
                            "${SRC}/kdsource/plist.c"
                            "${SRC}/kdsource/geom.c"
                            "${SRC}/kdsource/utils.c")
target_include_directories(kdsource PUBLIC "${SRC}/kdsource" "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(kdsource PUBLIC m LibXml2::LibXml2 mcpl)
install(TARGETS kdsource DESTINATION lib)
install(FILES "${SRC}/kdsource/kdsource.h"
              "${SRC}/kdsource/plist.h"
              "${SRC}/kdsource/geom.h"
              "${SRC}/kdsource/utils.h"
              "${CMAKE_CURRENT_BINARY_DIR}/KDSourceConfig.h"
              DESTINATION include)

install(PROGRAMS "${SRC}/kdtool/kdtool.sh" DESTINATION bin RENAME kdtool)
install(PROGRAMS "${SRC}/kdtool/templates.sh" DESTINATION bin RENAME kdtool-templates)
add_executable(kdtool-resample "${SRC}/kdtool/resample.c")
add_executable(kdtool-beamtest "${SRC}/kdtool/beamtest.c")
target_link_libraries(kdtool-resample PUBLIC kdsource)
target_link_libraries(kdtool-beamtest PUBLIC kdsource)
install(TARGETS kdtool-resample kdtool-beamtest DESTINATION bin)

add_subdirectory(tests)

install(DIRECTORY docs mcstas templates DESTINATION .)

if (INSTALL_PY)
	install(DIRECTORY "python/kdsource" DESTINATION python)
endif()
