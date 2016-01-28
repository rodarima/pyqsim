#!/bin/bash

FILE=$1
TEX=${FILE}.tex
PDF=${FILE}.pdf


while [ 1 ]; do
	# Files to look can change
	FILES=$(cat "$TEX" |\
		grep -o '\\input{.*}' |\
		sed -e 's/\\input{\(.*\)}/\1.tex/g' |\
		paste -sd " ")

	inotifywait -e modify $TEX $FILES

	pdflatex -interaction=nonstopmode -halt-on-error $TEX
	if [ "$?" != "0" ]; then
		echo -e "\a"
		continue
	fi
	pdflatex -interaction=nonstopmode -halt-on-error $TEX
	if [ "$?" != "0" ]; then
		echo -e "\a"
		continue
	fi
	kill -HUP $(pidof mupdf) $(pidof mu)
done
