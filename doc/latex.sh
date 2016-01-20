#!/bin/bash

FILE=$1
TEX=${FILE}.tex
PDF=${FILE}.pdf

while [ 1 ]; do
	inotifywait -e modify $TEX
	pdflatex -interaction=nonstopmode -halt-on-error $TEX
	if [ "$?" != "0" ]; then
		echo -e "\a"
		continue
	fi
	kill -HUP $(pidof mupdf)
	kill -HUP $(pidof mu)
done
