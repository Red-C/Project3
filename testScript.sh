#!/bin/bash

files=$(ls ./assignment3-tests-and-results/ | grep '[a-zA-Z]\.txt')

outputname=project3

g++ $1 $2 -o $outputname ./template-rt.cpp  

for file in $files 
do
	./$outputname './assignment3-tests-and-results/'$file
done

