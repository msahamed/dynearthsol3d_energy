#!/bin/bash
# Script to remove CMake-generated files from git tracking

echo "Removing CMake-generated files from git tracking without deleting the files..."

# Get a list of all CMake-generated files tracked by git
TRACKED_FILES=$(git ls-files | grep -iE 'CMakeFiles|CMakeCache.txt|cmake_install.cmake|Makefile|\.o\.d$|flags\.make|link\.txt|build\.make|DependInfo\.cmake|depend\.make|cmake_clean\.cmake|compiler_depend\.|\.cmake$')

# Remove files from git tracking without deleting them
if [ -n "$TRACKED_FILES" ]; then
    echo "$TRACKED_FILES" | xargs git rm --cached
    echo "Files have been removed from git tracking but kept in the filesystem."
    echo "These files will be ignored in future commits due to updates to .gitignore."
    echo ""
    echo "Number of files removed from tracking: $(echo "$TRACKED_FILES" | wc -l)"
    echo ""
    echo "To complete the process, commit these changes:"
    echo "git commit -m \"Remove CMake-generated files from git tracking\""
else
    echo "No CMake-generated files found in git tracking."
fi

echo "Done!"