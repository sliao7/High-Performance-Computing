#!/bin/bash
for n in 100 500 1000 5000 10000
do
for p in 1
do
for d in 1 5 9
do
	qsub ${n}A${p}A${d}A.pbs
done
done
done