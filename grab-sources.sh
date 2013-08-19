#! /bin/bash

set -e

function download_unpack_zip_if_not_there
{
  URL="$1"
  BASENAME="$(basename $1)"
  DIRNAME="${BASENAME%.zip}"
  UNVERSIONED_DIRNAME="${DIRNAME%-[0-9].*}"

  if ! test -f $BASENAME; then
    wget $URL || curl -O $URL
  fi
  if ! test -d $UNVERSIONED_DIRNAME; then
    unzip $BASENAME
    mv $DIRNAME $UNVERSIONED_DIRNAME
  fi
}

download_unpack_zip_if_not_there https://cims.nyu.edu/cmcl/fmm2dlib/fmmlib2d-1.2.zip
download_unpack_zip_if_not_there https://cims.nyu.edu/cmcl/fmm3dlib/fmmlib3d-1.2.zip
