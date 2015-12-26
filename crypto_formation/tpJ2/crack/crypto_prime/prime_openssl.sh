#!/bin/bash

# define start and ending points
AQUO=15114307486000000000
ADQUEM=15114307486330924284
for N in $(seq $AQUO $ADQUEM); do
  # use bc to convert hex to decimal
  openssl prime $N | awk '/is prime/ {print "ibase=16;"$1}' | bc
done
