#!/bin/bash

TARGET_NAME="G.12_foldmason" 

/Users/dalsasso/foldmason/bin/foldmason easy-msa <path2folder_with_selected_files>/*pdb "$TARGET_NAME" ./tmp --report-mode 1

rm -r ./tmp
