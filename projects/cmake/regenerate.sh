#!/bin/sh

set -e

SCRIPT_DIR=$(pwd)

if [ "$#" = 1 ] ; then
    BUILD_DIR="$1"
else
    BUILD_DIR=build
fi

echo "Making build directory '$BUILD_DIR'"
if [ ! -e "$BUILD_DIR" ] ; then
   mkdir "$BUILD_DIR"
fi

BUILD_DIR=$(cd "$BUILD_DIR"; pwd)

SRC_DIR="$BUILD_DIR"/../../../src
cd "$BUILD_DIR"/../../../src

echo "Generating src/generated_include_dirs.cmake"
######### Generate generated_include_dirs.cmake ############
list_include_dirs()
{
    find "$@" -type d |
        grep -v "svn" |
        sed "s|libs/\([^/]*\)/.*|libs/\1|" |
        sort |
        uniq |
        sed 's|^|    ${PROJECT_SOURCE_DIR}/|g'
}

(
    echo '
# rb-libs and rb-core see only these, so they build without the language
set(RB_CORE_INCLUDE_DIRS'
    list_include_dirs libs core
    echo ' ${Boost_INCLUDE_DIRS} )

set(RB_LANGUAGE_INCLUDE_DIRS'
    list_include_dirs revlanguage
    echo ')
'
) > "$SRC_DIR/generated_include_dirs.cmake"

######## Generate CMakeLists.txt for subdirs
generate_subdir_cmake()
{
    local subdir="$1"
    local libname=$2
    echo "Generating src/$subdir/CMakeLists.txt"
    (
        echo "set(${subdir}_FILES"
        ### Treat all files in each subdir as source files
        find "$subdir" | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' 
        echo ")
add_library($libname \${${subdir}_FILES})"
    ) > "$SRC_DIR/$subdir/CMakeLists.txt"
}

generate_subdir_cmake "libs" "rb-libs"

generate_subdir_cmake "core" "rb-core"

generate_subdir_cmake "revlanguage" "rb-parser"

# We will only USE this if we do subdir(help2yml)
generate_subdir_cmake "help2yml" "rb-help"

# We will only USE this if we do subdir(cmd)
generate_subdir_cmake "cmd" "rb-cmd-lib"
