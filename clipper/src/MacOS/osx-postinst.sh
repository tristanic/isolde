#!/bin/bash
# Relativize dylib paths using @loader_path, for OSX 10.4+.
# If JHBUILD_PREFIX is not set the install prefix must be given as arg. 

set -eu -o pipefail

# set to 1 for more verbose output
VERBOSE=${VERBOSE:-0}
# set to 1 to print the output without actually changing anything
DRY_RUN=${DRY_RUN:-0}
# set to 1 to change install name of GCC dylibs
COPIED_COMPILER_LIBS=${COPIED_COMPILER_LIBS:-0}

if [ $# = 1 ]; then
  cd $1 || exit 1
elif [ "$JHBUILD_PREFIX" != "" ]; then
  cd $JHBUILD_PREFIX || exit 1
else
  echo "Usage: $0 install_prefix"
  echo "   or  ./cj run $0"
  exit 1
fi
prefix=`pwd`

# 2 args: install name and new prefix
relative_name()
{
  [ "$#" = 2 ] || exit 1
  case "$1" in
    ${prefix}*) # absolute path (vast majority of libs)
      echo $1 | sed "s|$prefix|$2|" | sed "s|/_jhbuild/.*$prefix||"
      ;;
    lib[a-zA-Z]*) # only filename (boost 1.56, LAPACK, BLT)
      echo "$2/lib/$1"
      ;;
    @rpath/libboost_*) # Boost 1.60
      echo "$2/lib/$(basename $1)"
      ;;
    /opt/X11/lib/*)
      echo $1 | sed "s|/opt/X11|/usr/X11|"
      ;;
    # gfortran from dmg has libs in /usr/local/gfortran/lib
    # homebrew has GCC libs in /usr/local/opt/gcc/lib/gcc/ and
    # /usr/local/lib/gcc/, MacPorts in /opt/local/lib/libgcc,
    # the last line is for Intel compiler libs
    /usr/local/gfortran/lib/lib* |\
    /usr/local/gfortran/lib/gcc/5/lib* |\
    /usr/local/opt/gcc/lib/gcc/* |\
    /usr/local/lib/gcc/* |\
    /opt/local/lib/libgcc/lib* |\
    */mac64/lib[^/]*.dylib | */libiomp*.dylib | */libirc*.dylib)
      if [ "$COPIED_COMPILER_LIBS" = 1 ]; then echo "$2/lib/$(basename $1)"; fi
      ;;
    /System/Library/Frameworks/Python.framework/Versions/2.7/Python)
      path="$2/lib/libpython2.7.dylib"
      if [ -e "$path" ]; then echo "$path"; fi
      ;;
    replace_with_path_to_ccp4/*)
      echo $1 | sed "s|replace_with_path_to_ccp4|$2|"
      ;;
    # distutils-built libs in MG use build paths of cmake-built dylibs
    /*/ccp4mg/lib[^/]*.dylib)
      echo "$2/lib/$(basename $1)"
      ;;
  esac
}

subst_loader_path()
{
  loader_dir="$(dirname "$2")"
  echo "${1/@loader_path/$loader_dir}"
}

warn_if_absent()
{
  [ "$#" = 2 ] || exit 1
  path=$(subst_loader_path "$1" "$2")
  if [ -e "$path" ]; then
    if [[ "$path" =~ ^/usr/local|^/opt|^/sw ]]; then
      echo ":::INFO::: local $1 referenced from $2"
    elif [[ "$path" =~ ^/System/Library/Frameworks/Python.framework ]]; then
      echo ":::INFO::: $1 referenced from $2"
    fi
  else
    if [[ "$path" =~ ^@ ]]; then
      echo ":::INFO::: $1 referenced from $2"
    else
      echo ":::WARNING::: missing $1 referenced from $2"
    fi
  fi
}

# change dependency paths in binaries
for b in bin/* libexec/* $(find lib \( -name \*.so -o -name \*.dylib \) ); do
 [ -h $b ] && continue  # skip symlinks, to make it faster
 if (file $b | grep -q 'Mach-O'); then
  [ $VERBOSE != 0 ] && /bin/echo $b || /bin/echo -n o
  [ $DRY_RUN = 0 ] && chmod u+w $b
  rel="@loader_path/$(dirname $b | sed -E 's=[^/]+=..=g')"
  lib_id=$(otool -D $b | sed '2!d;q')
  if [[ "$lib_id" =~ ^/ ]]; then
   new=$(relative_name $lib_id $rel)
   echo "$new"
   if [ -n "$new" ] && [ "$lib_id" != "$new" ]; then
    #if [ !-e ${new/replace_with_path_to_ccp4/.} ]; then
     [ $VERBOSE != 0 ] && /bin/echo "    \`-> $new" || /bin/echo -n ':'
     echo "Changing id to $new"
     [ $DRY_RUN = 0 ] && install_name_tool -id $new $b
     lib_id="$new"
    #else
     #echo "BZZZZZT no such file: $new"
    #fi
   fi
  fi

  for old in $(otool -L $b | tail +2 | awk '{print $1}'); do
   # -L displays also ID for libraries, skip them
   [ "$old" = "$lib_id" ] && continue
   new=$(relative_name $old $rel)
   # too long string may not fit - removing redundant ../lib may help
   [[ $b =~ ^lib/ && $new =~ ^@loader_path ]] && new=${new/\/..\/lib\//\/}
   if [ -n "$new" ] && [ -e "$(subst_loader_path "$new" "$b")" ]; then
    [ $VERBOSE != 0 ] && /bin/echo "    $old -> $new" || /bin/echo -n .
    [ $DRY_RUN = 0 ] && install_name_tool -change $old $new $b
    warn_if_absent $new $b
   else
    warn_if_absent $old $b
   fi
  done
 fi
done
echo
