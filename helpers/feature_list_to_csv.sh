#!/bin/bash

# Generate a file listing CSV of the current folder as basis for the annotation file.
# TODO: change regular expression accordingly.

find . -type f -name "*_featureScores.csv" | while read filepath; do
    filename=$(basename "$filepath" "_featureScores.csv")
    fullpath=$(realpath "$filepath")
    echo "$filename,$fullpath"
done > feature_list_paths.csv
