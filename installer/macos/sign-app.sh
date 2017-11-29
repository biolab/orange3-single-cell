#!/usr/bin/env bash

# Build app
rm -rf ~/dev/orange3-single-cell/dist/scOrange.app
./build-macos-app.sh ~/dev/orange3-single-cell/dist/scOrange.app
# Missing symlink Current messes with code signing
ln -s 3.6 ../../dist/scOrange.app/Contents/Frameworks/Python.framework/Versions/Current
# sign bundle
codesign -s "Developer ID" /Users/anze/dev/orange3-single-cell/dist/scOrange.app/Contents/Frameworks/Python.framework/Versions/3.6
codesign -s "Developer ID" /Users/anze/dev/orange3-single-cell/dist/scOrange.app/Contents/MacOS/pip
codesign -s "Developer ID" /Users/anze/dev/orange3-single-cell/dist/scOrange.app

# Create disk image
./create-dmg-installer.sh --app ../../dist/scOrange.app ../../dist/scOrange-1.0.0-alpha.dmg
# Sign disk image
codesign -s "Developer ID" /Users/anze/dev/orange3-single-cell/dist/scOrange-1.0.0-alpha.dmg
