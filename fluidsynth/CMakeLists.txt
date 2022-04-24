# For more information about using CMake with Android Studio, read the
# documentation: https://d.android.com/studio/projects/add-native-code.html

# Sets the minimum version of CMake required to build the native library.

cmake_minimum_required(VERSION 3.6.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 \
                    -frtti -fexceptions -Wall -Wno-unused-const-variable \
                    -Wno-unused-variable -Wno-unused-private-field \
                    -Wno-constant-conversion -Wno-infinite-recursion \
                    -Wno-deprecated-declarations -Wno-reorder \
                    -Wno-sometimes-uninitialized -Wno-unused-function \
                    -Wno-tautological-constant-out-of-range-compare \
                    -Wno-implicit-function-declaration \
                    -DUSE_EQUALIZER=1 \
")

# Creates and names a library, sets it as either STATIC
# or SHARED, and provides the relative paths to its source code.
# You can define multiple libraries, and CMake builds them for you.
# Gradle automatically packages shared libraries with your APK.

include_directories(
        ./
)

file(GLOB SOURCE_FILES ./*.cpp)
file(GLOB HEADERS_FILES ./*.h)

set(FLUIDSYNTH_SOURCE_DIR src)
include(${FLUIDSYNTH_SOURCE_DIR}/CMakeLists.txt)

set(FLUIDSYNTH_SOURCE_FILES
        ${config_SOURCES}
        ${fluid_dbus_SOURCES}
        ${fluid_lash_SOURCES}
        ${fluid_ladspa_SOURCES}
        ${fluidsynth_SOURCES}
        ${public_HEADERS}
        ${public_main_HEADER}
        )

add_library( # Sets the name of the library.
        FluidSynthLib

        # Sets the library as a shared library.
        SHARED

        # Provides a relative path to your source file(s).
        ${SOURCE_FILES}
        ${HEADERS_FILES}
        ${FLUIDSYNTH_SOURCE_FILES}
        )

target_compile_options(FluidSynthLib PRIVATE ${FLUID_COMPILE_OPTIONS})

target_include_directories(FluidSynthLib PRIVATE
        ${FLUID_INCLUDE_DIRECTORIES}
        )

#target_link_directories(FluidSynthLib PRIVATE
#        ${FLUID_LINK_DIRECTORIES}
#)

#target_include_directories(FluidSynthLib INTERFACE
#        ${FLUIDSYNTH_INTERFACE_INCLUDE_DIRECTORIES}
#        )

# Searches for a specified prebuilt library and stores the path as a
# variable. Because CMake includes system libraries in the search path by
# default, you only need to specify the name of the public NDK library
# you want to add. CMake verifies that the library exists before
# completing its build.

find_library( # Sets the name of the path variable.
        log-lib

        # Specifies the name of the NDK library that
        # you want CMake to locate.
        log)

# Specifies libraries CMake should link to your target library. You
# can link multiple libraries, such as libraries you define in this
# build script, prebuilt third-party libraries, or system libraries.

target_link_libraries(FluidSynthLib
        PRIVATE
        # Links the target library to the log library
        # included in the NDK.
        android
        atomic
        ${FLUID_LINK_LIBRARIES}
        ${log-lib})