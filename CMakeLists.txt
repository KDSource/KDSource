
##################################################################################
#                                                                                #
# CMake file to be used to install KSource distribution.                         #
#                                                                                #
# See README.md file for installation instructions.                              #
#                                                                                #
# Written 2021 by O. I. Abbate.                                                  #
#                                                                                #
##################################################################################

cmake_minimum_required(VERSION 3.0.0)

project(KSource VERSION 1.0.0)
configure_file(KSourceConfig.h.in KSourceConfig.h)

enable_testing()

set(INSTALL_PY OFF CACHE STRING "Whether to also install ksource and mcpl python files.")


set(CMAKE_C_STANDARD 99)

set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")

find_package(LibXml2 REQUIRED)

set(BUILD_WITHG4 OFF CACHE STRING "Whether to build Geant4 plugins if Geant4 is available.")
add_subdirectory(mcpl)

add_library(ksource SHARED "${SRC}/ksource/ksource.c"
                           "${SRC}/ksource/plist.c"
                           "${SRC}/ksource/geom.c"
                           "${SRC}/ksource/utils.c")
target_include_directories(ksource PUBLIC "${SRC}/ksource" "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(ksource PUBLIC m LibXml2::LibXml2 mcpl)
install(TARGETS ksource DESTINATION lib)
install(FILES "${SRC}/ksource/ksource.h"
              "${SRC}/ksource/plist.h"
              "${SRC}/ksource/geom.h"
              "${SRC}/ksource/utils.h"
              "${CMAKE_CURRENT_BINARY_DIR}/KSourceConfig.h"
              DESTINATION include)

install(PROGRAMS "${SRC}/kstool/kstool.sh" DESTINATION bin RENAME kstool)
install(PROGRAMS "${SRC}/kstool/templates.sh" DESTINATION bin RENAME kstool-templates)
add_executable(kstool-resample "${SRC}/kstool/resample.c")
add_executable(kstool-beamtest "${SRC}/kstool/beamtest.c")
target_link_libraries(kstool-resample PUBLIC ksource)
target_link_libraries(kstool-beamtest PUBLIC ksource)
install(TARGETS kstool-resample kstool-beamtest DESTINATION bin)

add_subdirectory(tests)

install(DIRECTORY docs examples mcstas templates DESTINATION .)

if (INSTALL_PY)
	install(DIRECTORY "python/ksource" DESTINATION python)
endif()
