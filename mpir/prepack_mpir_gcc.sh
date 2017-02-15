#!/bin/bash

die() {
    echo >&2 "$@"
    exit 1
}

if [ $# -ne 2 ]; then
    echo "Please specify an arch and a package revision number."
    exit
fi

name=mpir
version=2.2.1
arch=$1
revision=$2
archive=mpir-${version}

case $arch in
  Linux*)
    tar=tar
    ;;
  MacX)
    tar=gnutar

    # Disable extra resource fork ._* files.
    export COPYFILE_DISABLE
    export COPY_EXTENDED_ATTRIBUTES_DISABLE
    ;;
  *)
    die "Unsupported platform"
    ;;
esac  

url="http://www.mpir.org/mpir-${version}.tar.bz2"

[ -e ${archive}.tar.bz2 ] || curl --progress-bar --fail -O $url || die "Failed to download '$url'."

[ -d ${archive} ] || tar xvf ${archive}.tar.bz2 || die "Failed to extract ${archive}.tar.bz2"

cd ${archive}

./configure --enable-gmpcompat --disable-shared --enable-cxx --prefix=`pwd`/build || die "Failed to configure"

make -j12 || die "Build failed"
make check || die "Check failed"
make install  || die "Install failed"

cd ..

GCC_VERSION=$(gcc -dumpversion)
MAJOR_GCC_VERSION=$(echo $GCC_VERSION | cut -d'.' -f1)
MINOR_GCC_VERSION=$(echo $GCC_VERSION | cut -d'.' -f2)
PATCH_GCC_VERSION=$(echo $GCC_VERSION | cut -d'.' -f3)
compiler=gcc${MAJOR_GCC_VERSION}${MINOR_GCC_VERSION}

basename=${name}_${version}-${revision}

echo "Packaging 'base' ..."
for cfg in Debug Optimize
do
    echo "   ... '${cfg}' ..."
    mkdir -p lib/arch-${arch}-${cfg}

    cp -a ${archive}/build/lib/*.a lib/arch-${arch}-${cfg}
    
    filename=${basename}_dev_${arch}-${compiler}-$cfg
    $tar --numeric-owner --owner=0 --group=0 -cvjf ${filename}.tar.bz2 lib
    rm -fr lib
done

echo "Packaging 'dev' ..."
cp -rf ${archive}/build/include .
$tar --numeric-owner --owner=0 --group=0 -cvjf ${basename}_dev_${arch}.tar.bz2 include
rm -fr include


cat <<\EOF
========================================================================

Prepackeds are available in the current directory.
Add the following lines to PREPACKED in mpir/Package:

EOF

openssl sha1 *.tar.bz2 | sed -e 's/.tar.bz2//' -e 's/SHA1.//' -e 's/.= /:sha1-/'










