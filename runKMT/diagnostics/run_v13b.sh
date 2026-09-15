#!/bin/bash
SP=$HOME/runKMT/$(basename $(dirname $(readlink -f $0)))
cd /home/ocs
for N in 8 12; do for F in BLG41 BLG01; do
  ( FIELD=$F NSYS=$N matlab -batch "run('$SP/v13chain.m')" > $SP/v13c${N}_$F.log 2>&1 ) &
done; done
wait
echo "V13 ALL FINISHED"
