#!/bin/bash

tsv_directory="$1"
output_file="$2"

cd "${tsv_directory}" || exit
for directory in *; do
    if [ -d "${directory}" ]; then
        echo "Appending ${directory}.tsv contents..."
        cat "${directory}/${directory}.tsv" >> "${output_file}"
    fi
done