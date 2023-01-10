#!/bin/bash

while read TMP
do
	cd parts/$TMP/
	codeml ${TMP}.ctl
	cd ../../
done < ctls.txt
