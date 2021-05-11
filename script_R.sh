#!bin/bash

awk '{print $2}' MOVIE.xyz | grep -v '^ *$' | awk 'NR % 2 == 0' > R.dat
