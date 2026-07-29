cat 11 | awk '{print $2}' | capital.sh 1 1 | awk '{print $1".bed.gz"}' > 22
cat 1 22
cat 1 22 | break_count.sh s_2 | sort
