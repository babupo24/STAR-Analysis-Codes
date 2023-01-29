#!/bin/bash

 for run in `sed -n  11,601p KevinsRunList.list` 
do
	star-submit-template -template ./add_trees.xml -entities run=$run 

done


