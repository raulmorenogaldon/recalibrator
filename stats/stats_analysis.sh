#!/bin/sh
# $1 = SCHEDULE
# $2 = BATCH SIZE

echo "Recollecting..."
if ./stats_recollect.sh $1 $2
then
	echo "Recollection done!"
	echo "Plotting..."
	./stats_plot.sh results_$1_$2_.stats
	echo "Successfull plot! Output: results_$1_$2_.stats"
else
	echo "Failed recollection!!!"
fi
