#!/bin/sh
# 
# Make the .tar.gz with a version number.
#

tar -zcvf afnl--stable--$1.tar.gz LICENSE Makefile doc/*.tex \
doc/*.bbl doc/*.idx doc/lib.pdf src/*.f90 doc/CC*


