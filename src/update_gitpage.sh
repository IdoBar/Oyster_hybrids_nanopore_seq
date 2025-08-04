#!/bin/bash 
cp *_analysis.html index.html
git add .
git commit -m "$1"
git push