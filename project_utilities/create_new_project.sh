#!/bin/bash

projectName=$1

baseDir="/groups/itay_mayrose/nosnap/liorglic/Projects"
projectDir="$baseDir/$projectName"
mkdir $projectDir
mkdir "$projectDir/R"
mkdir "$projectDir/data"
mkdir "$projectDir/doc"
mkdir "$projectDir/figs"
mkdir "$projectDir/output"
mkdir "$projectDir/notebooks"
createDate=`date`
echo "Project $projectName - created $createDate" > "$projectDir/README"
