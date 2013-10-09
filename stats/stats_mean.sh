#!/bin/sh
# $1 = FILE

#cat $1 | awk 'BEGIN {max_r = 0.0; max_p = 0.0; max_w = 0.0}
#	{sum_r += $1; sum_p += $2; sum_w += $3}
#	$1 > max_r { max_r = $1 }
#	$2 > max_p { max_p = $2 }
#	$3 > max_w { max_w = $3 }
#	END { printf "%.2f %.2f %.2f %.2f %.2f %.2f\n", 
#	sum_r/NR, sum_p/NR, sum_w/(NR-2), max_r, max_p, max_w}'

cat $1 | awk '
	BEGIN {
		maxindex = -1000
		minindex = 1000
	}

	{ 
		sums[$1] += $2
		num[$1]++

		if ($2 > max[$1])
			max[$1] = $2

		if ($1 > maxindex)
			maxindex = $1

		if ($1 < minindex)
			minindex = $1
	}

	END {
		#Print times
		for (x = minindex; x <= maxindex; x++)
		{
			printf "%d %.9f\n", x, sums[x] / num[x]
		}
	}'
