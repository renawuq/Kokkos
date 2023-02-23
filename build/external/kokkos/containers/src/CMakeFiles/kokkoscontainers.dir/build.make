# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.23.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.23.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build

# Include any dependencies generated for this target.
include external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/compiler_depend.make

# Include the progress variables for this target.
include external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/progress.make

# Include the compile flags for this target's objects.
include external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/flags.make

external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.o: external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/flags.make
external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.o: ../external/kokkos/containers/src/impl/Kokkos_UnorderedMap_impl.cpp
external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.o: external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.o"
	cd /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/containers/src && /usr/local/bin/g++-10 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.o -MF CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.o.d -o CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.o -c /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/external/kokkos/containers/src/impl/Kokkos_UnorderedMap_impl.cpp

external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.i"
	cd /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/containers/src && /usr/local/bin/g++-10 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/external/kokkos/containers/src/impl/Kokkos_UnorderedMap_impl.cpp > CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.i

external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.s"
	cd /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/containers/src && /usr/local/bin/g++-10 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/external/kokkos/containers/src/impl/Kokkos_UnorderedMap_impl.cpp -o CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.s

# Object files for target kokkoscontainers
kokkoscontainers_OBJECTS = \
"CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.o"

# External object files for target kokkoscontainers
kokkoscontainers_EXTERNAL_OBJECTS =

external/kokkos/containers/src/libkokkoscontainers.a: external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/impl/Kokkos_UnorderedMap_impl.cpp.o
external/kokkos/containers/src/libkokkoscontainers.a: external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/build.make
external/kokkos/containers/src/libkokkoscontainers.a: external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libkokkoscontainers.a"
	cd /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/containers/src && $(CMAKE_COMMAND) -P CMakeFiles/kokkoscontainers.dir/cmake_clean_target.cmake
	cd /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/containers/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/kokkoscontainers.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/build: external/kokkos/containers/src/libkokkoscontainers.a
.PHONY : external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/build

external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/clean:
	cd /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/containers/src && $(CMAKE_COMMAND) -P CMakeFiles/kokkoscontainers.dir/cmake_clean.cmake
.PHONY : external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/clean

external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/depend:
	cd /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/external/kokkos/containers/src /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/containers/src /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/kokkos/containers/src/CMakeFiles/kokkoscontainers.dir/depend
