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

# Utility rule file for NightlyStart.

# Include any custom commands dependencies for this target.
include external/kokkos/CMakeFiles/NightlyStart.dir/compiler_depend.make

# Include the progress variables for this target.
include external/kokkos/CMakeFiles/NightlyStart.dir/progress.make

external/kokkos/CMakeFiles/NightlyStart:
	cd /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos && /usr/local/Cellar/cmake/3.23.2/bin/ctest -D NightlyStart

NightlyStart: external/kokkos/CMakeFiles/NightlyStart
NightlyStart: external/kokkos/CMakeFiles/NightlyStart.dir/build.make
.PHONY : NightlyStart

# Rule to build all files generated by this target.
external/kokkos/CMakeFiles/NightlyStart.dir/build: NightlyStart
.PHONY : external/kokkos/CMakeFiles/NightlyStart.dir/build

external/kokkos/CMakeFiles/NightlyStart.dir/clean:
	cd /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos && $(CMAKE_COMMAND) -P CMakeFiles/NightlyStart.dir/cmake_clean.cmake
.PHONY : external/kokkos/CMakeFiles/NightlyStart.dir/clean

external/kokkos/CMakeFiles/NightlyStart.dir/depend:
	cd /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/external/kokkos /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos /Users/mengxiang/Desktop/prework/SCEC-Kokkos-Project/acoustic_kokkos/build/external/kokkos/CMakeFiles/NightlyStart.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/kokkos/CMakeFiles/NightlyStart.dir/depend

