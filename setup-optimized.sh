#! /bin/sh

if test "$1" = ""; then
  echo "usage: $0 <command> ..."
  echo "<command> will usually be 'install' or 'develop'."
  exit 1
fi

FOPT="-O3 -fopenmp" OPT="-O3 -fopenmp" EXTRA_LINK_ARGS="-fopenmp" python setup.py "$1" "$@"
