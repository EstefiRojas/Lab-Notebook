#!/bin/bash

set -uex


for histone_name in $(ls ../data/epigenetic_data/histone_marks/); do
	sh run_histone.sh "$histone_name"
done
