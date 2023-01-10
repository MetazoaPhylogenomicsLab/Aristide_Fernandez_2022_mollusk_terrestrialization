#!/bin/bash

# this script moves files and modifies tmp.ctl files from mcmctree run with usedata=3 to then run codeml properly with WAG model

ls res_1stpart/*.ctl > a.txt
sed 's/\.ctl//' a.txt > b.txt
sed 's/res_1stpart\///' b.txt > ctls.txt
rm a.txt b.txt

mkdir parts

while read TMP
do
	mkdir parts/$TMP
	cp res_1stpart/${TMP}.ctl parts/$TMP
	cp res_1stpart/${TMP}.trees parts/$TMP
	cp res_1stpart/${TMP}.txt parts/$TMP
	sed 's/model = 0/model = 2/' parts/$TMP/${TMP}.ctl > parts/$TMP/tmp.ctl
	sed 's/aaRatefile =/aaRatefile = wag.dat/' parts/$TMP/tmp.ctl > parts/$TMP/tmp2.ctl
	mv parts/$TMP/tmp2.ctl parts/$TMP/${TMP}.ctl
	rm parts/$TMP/tmp.ctl
	cp wag.dat parts/$TMP
done <  ctls.txt 
