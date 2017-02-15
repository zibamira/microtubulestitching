#!/bin/bash

die() {
    echo >&2 "$@"
    exit 1
}

confirm() {
    echo "    Confirm with any key, or hit CTRL-C to break."
    read _
}

if [ $# -ne 2 ]; then
    echo "Please specify an arch and a package revision number."
    exit
fi

name=mpir
version=2.7.0
arch=$1
revision=$2
archive=mpir-${version}

url="http://www.mpir.org/mpir-${version}.tar.bz2"

[ -e ${archive}.tar.bz2 ] || curl --progress-bar --fail -O $url || die "Failed to download '$url'."

[ -d ${archive} ] || tar xvf ${archive}.tar.bz2 || die "Failed to extract ${archive}.tar.bz2"

scriptdir=`pwd`
if [ -d "$scriptdir/patches/$arch/$version" ]; then
    for p in $scriptdir/patches/$arch/$version/*.patch
    do
        echo "Applying patch $p ..."
        ( cd $archive ; patch --reject-file=$(basename $p).rej -p1 <"$p" ) || {
            echo "Failed to apply patch; continue anyway?"
            confirm
        }
    done
fi

cd ${archive}

case $arch in
    Win64VC12)
        platform="x64"
        env="%VS120COMNTOOLS%..\..\VC\bin\x86_amd64\vcvarsx86_amd64.bat"
        solution="build.vc12\\mpir.sln"
    ;;
    *)
        echo "Unsupported arch '$arch'."
        exit
    ;;
esac

for config in "Optimize" "Debug"; do
    if [ $config == "Optimize" ]; then
        cfg="Release"
    else
        cfg=$config
    fi
    
cat <<EOF > helper_mpir_win.cmd
call "$env"
devenv $solution /rebuild "$cfg|$platform" /project "lib_mpir_core2"
devenv $solution /rebuild "$cfg|$platform" /project "lib_mpir_cxx"
EOF

    cmd //c helper_mpir_win.cmd
    rm helper_mpir_win.cmd
done

cd ..

basename=${name}_${version}-${revision}

mkdir -p include
cp $archive/lib/$platform/Debug/*.h include
tar --numeric-owner --owner=0 --group=0 -cvjf ${basename}_dev_${arch}.tar.bz2 include
rm -fr include

mkdir -p lib/arch-${arch}-Debug
mkdir -p lib/arch-${arch}-Optimize
cp -a ${archive}/lib/$platform/Debug/*.lib lib/arch-${arch}-Debug
cp -a ${archive}/lib/$platform/Release/*.lib lib/arch-${arch}-Optimize
tar --numeric-owner --owner=0 --group=0 -cvjf ${basename}_dev_${arch}-Debug.tar.bz2 lib/arch-${arch}-Debug
tar --numeric-owner --owner=0 --group=0 -cvjf ${basename}_dev_${arch}-Optimize.tar.bz2 lib/arch-${arch}-Optimize
rm -fr lib

cat <<\EOF
========================================================================

Prepackeds are available in the current directory.
Add the following lines to PREPACKED in mpir/Package:

EOF

openssl sha1 *.tar.bz2 | sed -e 's/.tar.bz2//' -e 's/SHA1.//' -e 's/.= /:sha1-/'





