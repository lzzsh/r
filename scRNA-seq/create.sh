#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <id> <description>"
    exit 1
fi

id=$1
description=$2

src="../../new_folder"
dest="./$(date +%Y-%m-%d)_${id}_${description// /-}"

cp -r "$src" "$dest"
echo "Files copied to $dest"