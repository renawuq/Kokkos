# Install script for directory: /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/external/kokkos

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/core/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/containers/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/algorithms/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/example/cmake_install.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos" TYPE FILE FILES
    "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/KokkosConfig.cmake"
    "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/KokkosConfigCommon.cmake"
    "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/KokkosConfigVersion.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos/KokkosTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos/KokkosTargets.cmake"
         "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/CMakeFiles/Export/lib/cmake/Kokkos/KokkosTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos/KokkosTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos/KokkosTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos" TYPE FILE FILES "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/CMakeFiles/Export/lib/cmake/Kokkos/KokkosTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos" TYPE FILE FILES "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/CMakeFiles/Export/lib/cmake/Kokkos/KokkosTargets-debug.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/kokkos" TYPE FILE FILES "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/KokkosCore_config.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES
    "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/external/kokkos/bin/nvcc_wrapper"
    "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/external/kokkos/bin/hpcbind"
    "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/external/kokkos/bin/kokkos_launch_compiler"
    "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/temp/kokkos_launch_compiler"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/kokkos" TYPE FILE FILES
    "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/KokkosCore_config.h"
    "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/KokkosCore_Config_FwdBackend.hpp"
    "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/KokkosCore_Config_SetupBackend.hpp"
    "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/KokkosCore_Config_DeclareBackend.hpp"
    "/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/KokkosCore_Config_PostInclude.hpp"
    )
endif()

