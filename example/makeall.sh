#!/bin/sh

for sampledir in *
do
    if [ -d ${sampledir} ]; then
	cd ${sampledir} > /dev/null
	if [ -e makefile ]; then
	    make $1
	fi
	cd - > /dev/null
    fi
done

exit 0
