#!/usr/bin/env bash

DIST="$( cd "$(dirname "$0")/../../dist" ; pwd -P )"
APP=$DIST/scOrange.app
CERT="Developer ID Application: Univerza v Ljubljani"
# Build app
rm -rf $APP
./build-macos-app.sh $APP
# Missing symlink Current messes with code signing
ln -s 3.6 $APP/Contents/Frameworks/Python.framework/Versions/Current
# sign bundle
codesign -s "$CERT" $APP/Contents/Frameworks/Python.framework/Versions/3.6
codesign -s "$CERT" $APP/Contents/MacOS/pip
codesign -s "$CERT" $APP

VERSION=$($APP/Contents/MacOS/python -c 'import pkg_resources; print(pkg_resources.get_distribution("Orange3-SingleCell").version)')
# Create disk image
./create-dmg-installer.sh --app $APP $DIST/scOrange-$VERSION.dmg
# Sign disk image
codesign -s "$CERT" $DIST/scOrange-$VERSION.dmg
