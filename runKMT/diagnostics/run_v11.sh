#!/bin/bash
SP=$HOME/runKMT/$(basename $(dirname $(readlink -f $0)))
export SP
cd /home/ocs
for M in trim secz130; do for F in BLG41 BLG01; do
  ( FIELD=$F MODE=$M matlab -batch "run('$SP/v11chain.m')" > $SP/v11_${M}_$F.log 2>&1 ) &
done; done
wait
echo "V11 ALL FINISHED"
