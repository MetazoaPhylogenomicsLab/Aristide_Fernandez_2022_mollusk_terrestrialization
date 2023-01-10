#!/bin/bash

while read TMP
do
	cat parts/${TMP}/rst2 >> out.BV
done < ctls.txt