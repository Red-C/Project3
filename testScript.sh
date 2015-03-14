#!/bin/bash

files=$(ls ./assignment3-tests-and-results/ | grep '[a-zA-Z]\.txt')

outputname=project3
outputfolder=results
g++ $1 $2 -o $outputname ./template-rt.cpp  

if [ ! -d "$outputfolder" ]; then
	mkdir $outputfolder
fi


for file in $files 
do
	./$outputname './assignment3-tests-and-results/'$file
done

mv *.ppm ./$outputfolder

