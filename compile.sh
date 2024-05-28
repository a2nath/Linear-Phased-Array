#!/bin/sh

echo "\n"
echo "======================================================================================="
echo "                         Compile Script for building the model    "
echo "======================================================================================="
echo "\n"

CURRPWD=$(pwd) # current dir
OUT=$CURRPWD/model.out
LOG=$CURRPWD/last_compile.log

[ "$1" = "-release" ] && build="Release" || build="Debug"

rm $OUT 2> /dev/null

cppdir="$CURRPWD/src"
cpp="project.cpp"
include="-I$CURRPWD/inc"
#options="-ftime-report"
#options="-pass-exit-codes"
release_options="-Ofast"
debug_options="-g"


if [ "$1" = "-release" ]; then
	flags="$include $CURRPWD/$cppdir/$cpp -ljsoncpp -o $OUT $release_options -std=c++17 -Wall"
else
	flags="$include $CURRPWD/$cppdir/$cpp -ljsoncpp -o $OUT $debug_options -std=c++17 -Wall"
fi
echo "\tFlags:"
echo "-----------------------------------------------"
echo "\t$flags"
echo "\n"


echo "\tOutput:"
echo "-----------------------------------------------"
start=`date +%s`
g++ $flags 2>&1 | tee $LOG
end=`date +%s`


# g++ exit code instead of using g++ flag
gcc_exit_code=$(grep ": error:" $LOG)
if [ -n "$gcc_exit_code" ]; then
	echo "Error:"
	echo "-----------------------------------------------"
	echo $gcc_exit_code
fi



# give confirmation of compilation
if [ -f $OUT ]; then
	echo "\n\n"
	echo "======================================================================================="
	echo "       $build build compiled on $(date) and took $(( $end - $start )) seconds"
	echo "======================================================================================="
fi
