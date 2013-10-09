#!/bin/sh
# $1 = SCHEDULE
# $2 = BATCH SIZE

out_file="results_$1_$2_.stats"
aux_file="aux.tmp"
if [ -f $out_file ]
then rm $out_file
fi

count=0

echo "Recollection: Output file $out_file"

for file in $(ls $1_$2_*_.stats);
do
	echo "Recollection: Proccessing $file"	

	# First field is number of threads
	aux=$(echo $file | awk -F "_" '{printf "%d", $3}')

	# Next fields for every time slot
	./stats_mean.sh $file > $aux_file

	# Append number of threads to first column
	cat $aux_file | awk -v t=$aux '
	{
		printf "%d %d %.9f\n", t, $1, $2
	}' >> $out_file

	count=$((count+1))
done

if [ $count -eq 0 ] 
then
	echo "Recollection: No files to recollect!"
	exit 1
else
	echo "Recollection: $count files recollected"
fi

exit 0
