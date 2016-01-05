#!/bin/bash

FILE=$1
TEX=${FILE}.tex
PDF=${FILE}.pdf

while [ 1 ]; do
	inotifywait -e modify $TEX
	pdflatex -interaction=nonstopmode -halt-on-error $TEX && kill -HUP $(pidof mupdf)
done
