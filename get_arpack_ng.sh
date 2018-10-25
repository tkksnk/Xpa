# This script is based on https://gist.github.com/galanakis/4069435
# You need to clone the files from https://github.com/opencollab/arpack-ng in advance.
VERSION=3.1.2
BASENAME=arpack-ng #_$VERSION
#PACKAGE=$BASENAME.tar.gz
LIB=libarpack.a


# silently download the sources
#curl -sOL http://forge.scilab.org/index.php/p/arpack-ng/downloads/get/$PACKAGE

# unpack, delete the package and delete second.f and some empty files. Instead of second.f we will use second_NONE.f
#tar zxf $PACKAGE
#rm $PACKAGE
rm $BASENAME/UTIL/second.f
rm $BASENAME/SRC/dnaupe.f
rm $BASENAME/SRC/snaupe.f

PARENTDIR=`pwd`
cd $BASENAME

OUTLIB=$PARENTDIR/$LIB

# Different architectures
case `uname` in
	Darwin)
		echo Compiling on Mac OS X
		ifort -fast -nofor-main -staticlib -o $OUTLIB SRC/*.f UTIL/*.f
		;;
	*)
		echo compiling on `uname`
		ifort -c -fast -nofor-main SRC/*.f UTIL/*.f
		xiar rc $OUTLIB *.o
		;;
esac

cd $PARENTDIR
# deleting the source directory
#rm -rf $BASENAME
