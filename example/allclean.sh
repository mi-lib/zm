#!/bin/sh

for sampledir in *
do
    if [ -d ${sampledir} ]; then
	cd ${sampledir} > /dev/null
	make clean
	cd - > /dev/null
    fi
done

exit 0
