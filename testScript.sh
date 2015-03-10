#!/bin/bash

files=$(ls ./assignment3-tests-and-results/ | grep '[a-zA-Z]\.txt')

for file in $files 
do
	./a.out './assignment3-tests-and-results/'$file
done

