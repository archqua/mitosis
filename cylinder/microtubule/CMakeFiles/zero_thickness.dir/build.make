# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mitosis/cylinder/zero_thickness

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mitosis/cylinder/zero_thickness

# Include any dependencies generated for this target.
include CMakeFiles/zero_thickness.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/zero_thickness.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/zero_thickness.dir/flags.make

CMakeFiles/zero_thickness.dir/main.cpp.o: CMakeFiles/zero_thickness.dir/flags.make
CMakeFiles/zero_thickness.dir/main.cpp.o: main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitosis/cylinder/zero_thickness/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/zero_thickness.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/zero_thickness.dir/main.cpp.o -c /home/mitosis/cylinder/zero_thickness/main.cpp

CMakeFiles/zero_thickness.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/zero_thickness.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitosis/cylinder/zero_thickness/main.cpp > CMakeFiles/zero_thickness.dir/main.cpp.i

CMakeFiles/zero_thickness.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/zero_thickness.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitosis/cylinder/zero_thickness/main.cpp -o CMakeFiles/zero_thickness.dir/main.cpp.s

# Object files for target zero_thickness
zero_thickness_OBJECTS = \
"CMakeFiles/zero_thickness.dir/main.cpp.o"

# External object files for target zero_thickness
zero_thickness_EXTERNAL_OBJECTS =

zero_thickness: CMakeFiles/zero_thickness.dir/main.cpp.o
zero_thickness: CMakeFiles/zero_thickness.dir/build.make
zero_thickness: /usr/lib/libarmadillo.so
zero_thickness: CMakeFiles/zero_thickness.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mitosis/cylinder/zero_thickness/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable zero_thickness"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/zero_thickness.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/zero_thickness.dir/build: zero_thickness

.PHONY : CMakeFiles/zero_thickness.dir/build

CMakeFiles/zero_thickness.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/zero_thickness.dir/cmake_clean.cmake
.PHONY : CMakeFiles/zero_thickness.dir/clean

CMakeFiles/zero_thickness.dir/depend:
	cd /home/mitosis/cylinder/zero_thickness && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mitosis/cylinder/zero_thickness /home/mitosis/cylinder/zero_thickness /home/mitosis/cylinder/zero_thickness /home/mitosis/cylinder/zero_thickness /home/mitosis/cylinder/zero_thickness/CMakeFiles/zero_thickness.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/zero_thickness.dir/depend
