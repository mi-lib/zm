#!/bin/sh

make clean

TEST=`ls *test.c`
for testfile in $TEST
do
  testcase=`basename $testfile .c`
  echo "> Testing" $testcase
  make $testcase || exit 1
  ./$testcase || exit 1
done

make clean

echo "Done."

exit 0
