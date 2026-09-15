#!/bin/bash
SP=$HOME/runKMT/$(basename $(dirname $(readlink -f $0)))
cd /home/ocs
for F in BLG41 BLG01; do
  ( FIELD=$F matlab -batch "run('$SP/v14chain.m')" > $SP/v14_$F.log 2>&1 ) &
done
wait
echo "V14 BOTH FINISHED"
