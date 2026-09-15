#!/bin/bash
SP=$HOME/runKMT/$(basename $(dirname $(readlink -f $0)))
cd /home/ocs
for F in BLG41 BLG01; do
  ( FIELD=$F matlab -batch "run('$SP/v12chain.m')" > $SP/v12_$F.log 2>&1 ) &
done
wait
echo "V12 BOTH FINISHED"
