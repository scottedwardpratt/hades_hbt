# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.20.5/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.20.5/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/mariastefaniak/Desktop/Work/femtoscopy_workshop/hades_hbt/scott_run

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/mariastefaniak/Desktop/Work/femtoscopy_workshop/hades_hbt/scott_run

# Utility rule file for extern_hades_hbt.

# Include any custom commands dependencies for this target.
include CMakeFiles/extern_hades_hbt.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/extern_hades_hbt.dir/progress.make

CMakeFiles/extern_hades_hbt:
	cd /Users/mariastefaniak/Desktop/Work/femtoscopy_workshop/hades_hbt/software && make

extern_hades_hbt: CMakeFiles/extern_hades_hbt
extern_hades_hbt: CMakeFiles/extern_hades_hbt.dir/build.make
.PHONY : extern_hades_hbt

# Rule to build all files generated by this target.
CMakeFiles/extern_hades_hbt.dir/build: extern_hades_hbt
.PHONY : CMakeFiles/extern_hades_hbt.dir/build

CMakeFiles/extern_hades_hbt.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/extern_hades_hbt.dir/cmake_clean.cmake
.PHONY : CMakeFiles/extern_hades_hbt.dir/clean

CMakeFiles/extern_hades_hbt.dir/depend:
	cd /Users/mariastefaniak/Desktop/Work/femtoscopy_workshop/hades_hbt/scott_run && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mariastefaniak/Desktop/Work/femtoscopy_workshop/hades_hbt/scott_run /Users/mariastefaniak/Desktop/Work/femtoscopy_workshop/hades_hbt/scott_run /Users/mariastefaniak/Desktop/Work/femtoscopy_workshop/hades_hbt/scott_run /Users/mariastefaniak/Desktop/Work/femtoscopy_workshop/hades_hbt/scott_run /Users/mariastefaniak/Desktop/Work/femtoscopy_workshop/hades_hbt/scott_run/CMakeFiles/extern_hades_hbt.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/extern_hades_hbt.dir/depend

