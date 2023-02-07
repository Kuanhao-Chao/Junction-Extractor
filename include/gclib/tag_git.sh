#!/bin/bash -e
git checkout master
ver=$(fgrep '#define GCLIB_VERSION ' GBase.h)
ver=${ver#*\"}
ver=${ver%%\"*}
#git fetch --tags
if [[ "$1" == "delete" || "$1" == "del" ]]; then
  echo "Deleting tag v$ver .."
  git tag -d v$ver
  git push origin :refs/tags/v$ver
  exit
fi
git tag -a "v$ver" -m "release $ver"
git push --tags
