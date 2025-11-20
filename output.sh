#!/bin/bash

# Define the destination folders
FOLDERS="src"

(
  # Generate directory structure
  echo "### FILE TREE ###"
  # tree $FOLDERS
  for folder in $FOLDERS; do
    tree "$folder"
    # tree -L 1 "$folder"
  done

  # Generate file contents (code files only)
  echo -e "\n### FILE CONTENTS ###"
  find $FOLDERS -type f \( -name "*.hpp" -o -name "*.cpp" -o -name "*.h" -o -name "*.py" \) \
    -exec sh -c 'echo -e "\n--- $1 ---\n"; cat "$1"' sh {} \;
  # CHANGED: Added -maxdepth 1 before -type f to limit to top-level files only
  # find $FOLDERS -maxdepth 1 -type f \( -name "*.hpp" -o -name "*.cpp" -o -name "*.h" -o -name "*.py" \) \
  #   -exec sh -c 'echo -e "\n--- $1 ---\n"; cat "$1"' sh {} \;
) > output.txt