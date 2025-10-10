#!/bin/bash

CMAKE_FILE="../CMakeLists.txt"
RELEASE_DIR="../release"

# ======================================================
# 1. Extract the version number

nibrary_VERSION_MAJOR=$(grep "set(nibrary_VERSION_MAJOR" "$CMAKE_FILE" | awk '{print $2}' | tr -d '()' | tr -d '\r')
nibrary_VERSION_MINOR=$(grep "set(nibrary_VERSION_MINOR" "$CMAKE_FILE" | awk '{print $2}' | tr -d '()' | tr -d '\r')
nibrary_VERSION_PATCH=$(grep "set(nibrary_VERSION_PATCH" "$CMAKE_FILE" | awk '{print $2}' | tr -d '()' | tr -d '\r')

VERSION="${nibrary_VERSION_MAJOR}.${nibrary_VERSION_MINOR}.${nibrary_VERSION_PATCH}"

echo "Creating release $VERSION"

# ======================================================
# 2. Create the release directory
RELEASE_FOLDER="$RELEASE_DIR/nibrary_v$VERSION"

# Check if the release directory already exists
if [ -d "$RELEASE_FOLDER" ]; then
  echo "Error: Release directory already exists: $RELEASE_FOLDER"
  exit 1
fi

mkdir -p "$RELEASE_FOLDER"

# Check if directory creation was successful
if [ $? -ne 0 ]; then
  echo "Error: Could not create release directory: $RELEASE_FOLDER"
  exit 1
fi

echo "Created release directory: $RELEASE_FOLDER"

# ======================================================
# 3. Copy files and directories to the release directory
ITEMS_TO_COPY=(
  "../cmake/"
  "../external/"
  "../scripts/"
  "../src/"
  "../build_linux.sh"
  "../build_macOS.sh"
  "../build_win.bat"
  "../CMakeLists.txt"
  "../LICENSE.md"
  "../README.md"
)

for ITEM in "${ITEMS_TO_COPY[@]}"; do
  # Handle directories and files differently
  if [ -d "$ITEM" ]; then
    echo "Copying directory: $ITEM to $RELEASE_FOLDER"
    cp -r "$ITEM" "$RELEASE_FOLDER"
  else
    echo "Copying file: $ITEM to $RELEASE_FOLDER"
    cp "$ITEM" "$RELEASE_FOLDER"
  fi

  # Check for errors after each copy operation
  if [ $? -ne 0 ]; then
    echo "Error copying $ITEM to $RELEASE_FOLDER"
    exit 1
  fi
done

echo "Copied items to release directory."



# ======================================================
# 4. Remove unnecessary directories from within submodules
PROXSUITE_TEST_FOLDER="${RELEASE_FOLDER}/external/proxsuite/test"
echo "Removing proxsuite test folder: $PROXSUITE_TEST_FOLDER"
rm -rf $PROXSUITE_TEST_FOLDER

PROXSUITE_EXAMPLE_FOLDER="${RELEASE_FOLDER}/external/proxsuite/examples"
echo "Removing proxsuite example folder: $PROXSUITE_EXAMPLE_FOLDER"
rm -rf $PROXSUITE_EXAMPLE_FOLDER

SIMDE_TEST_FOLDER="${RELEASE_FOLDER}/external/simde/test"
echo "Removing simde test folder: $SIMDE_TEST_FOLDER"
rm -rf $SIMDE_TEST_FOLDER

DCM2NIIX_GIT_FOLDER="${RELEASE_FOLDER}/external/dcm2niix/.git"
echo "Removing dcm2niix .git folder: $DCM2NIIX_GIT_FOLDER"
rm -rf $DCM2NIIX_GIT_FOLDER

GEOGRAM_GIT_FOLDER="${RELEASE_FOLDER}/external/geogram/.git"
echo "Removing geogram .git folder: $GEOGRAM_GIT_FOLDER"
rm -rf $GEOGRAM_GIT_FOLDER



# ======================================================
# 5. Create a .zip archive of the release folder
echo "Creating zip archive of release folder..."
ZIP_FILE="$RELEASE_DIR/nibrary_v$VERSION.zip"
cd "$RELEASE_FOLDER/.." || exit
zip -q -r "$ZIP_FILE" "nibrary_v$VERSION"

# Check for errors
if [ $? -ne 0 ]; then
  echo "Error creating zip archive: $ZIP_FILE"
  exit 1
fi

echo "Created zip archive: $ZIP_FILE"

# ======================================================


exit 0
