#!/bin/bash 

START=$1
END=$2
STEP=$3

for ((i=$START;i<=$END;i+=$STEP))
do
	((j=$i+$STEP-1))
	echo $i $j
	printf -v x "%03d" $i
	printf -v y "%03d" $j
	sed -e s/AAA/$i/g -e s/BBB/$j/g -e s/CCC/$x/g -e s/DDD/$y/g < calc_avg_structure.config > temp.calc_avg_structure.config
	time ./calc_avg_structure.py temp.calc_avg_structure.config
done

