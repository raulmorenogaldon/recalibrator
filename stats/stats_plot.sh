#!/bin/sh

export DATAFILE="$1"

echo "Plotting: Launching GNUPlot..."

gnuplot <<\EOF

#file=`echo $DATAFILE`

set terminal postscript
set output "`echo $DATAFILE`.eps"
set autoscale
set title 'Total Execution Time Measurement'
set xlabel 'Number of threads'
set ylabel 'Execution time in s'

#Plot total times
plot 	"<awk '{if($2 == 1){print $1,$3}}' $DATAFILE" using 1:2 with linespoints title "D_SLOT_PH1_COLLECT_BAM" , \
	"<awk '{if($2 == 9){print $1,$3}}' $DATAFILE" using 1:2 with linespoints title "D_SLOT_PH2_RECALIBRATE"

#Plot Phase 1 times
set ylabel 'Execution time in ms'
set title 'PHASE 1 Time Measurement'
plot	"<awk '{if($2 == 2){print $1,$3}}' $DATAFILE" using 1:($2*1e3) with linespoints title "D_SLOT_PH1_READ_BATCH" , \
	"<awk '{if($2 == 3){print $1,$3}}' $DATAFILE" using 1:($2*1e3) with linespoints title "D_SLOT_PH1_COLLECT_BATCH" , \
	"<awk '{if($2 == 11){print $1,$3}}' $DATAFILE" using 1:($2*1e3) with linespoints title "D_SLOT_PH1_ITERATION"

#Plot alig recollect times
set title 'Alignment Recollection(ph1)/Recalibration(ph2) Time Measurement'
set ylabel 'Execution time in micros'
plot 	"<awk '{if($2 == 4){print $1,$3}}' $DATAFILE" using 1:($2*1e6) with linespoints title "D_SLOT_PH1_COLLECT_ALIG", \
	"<awk '{if($2 == 8){print $1,$3}}' $DATAFILE" using 1:($2*1e6) with linespoints title "D_SLOT_PH2_RECAL_ALIG"
 
#Plot Phase 2 times
set ylabel 'Execution time in ms'
set title 'PHASE 2 Time Measurement'
plot	"<awk '{if($2 == 5){print $1,$3}}' $DATAFILE" using 1:($2*1e3) with linespoints title "D_SLOT_PH2_READ_BATCH" , \
	"<awk '{if($2 == 6){print $1,$3}}' $DATAFILE" using 1:($2*1e3) with linespoints title "D_SLOT_PH2_PROCCESS_BATCH" , \
	"<awk '{if($2 == 7){print $1,$3}}' $DATAFILE" using 1:($2*1e3) with linespoints title "D_SLOT_PH2_WRITE_BATCH" , \
	"<awk '{if($2 == 12){print $1,$3}}' $DATAFILE" using 1:($2*1e3) with linespoints title "D_SLOT_PH2_ITERATION"

EOF
