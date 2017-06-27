#! /bin/sh

if test "$2" = ""; then
  echo "usage: $0 <debug|opt> <command> ..."
  echo "<command> will usually be 'install' or 'develop'."
  exit 1
fi

style="$1"
shift

case "$style" in

  opt)
    export FOPT="-O3 -fopenmp"
    export OPT="-O3 -fopenmp"
    export EXTRA_LINK_ARGS="-fopenmp"
    ;;
  debug)
    export FOPT="-g"
    export OPT="-g"
    export EXTRA_LINK_ARGS="-g"
    ;;
  *)
    echo "unsupported style: $1"
    exit 1
    ;;
esac

python setup.py "$@"
