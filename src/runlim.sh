#!/bin/bash
make limitSet &
make drawLim &
wait
SECONDS=0
build/limitSet od > /dev/null &
build/limitSet ee > /dev/null &
build/limitSet ed > /dev/null &
build/limitSet ec > /dev/null &
build/limitSet oc > /dev/null &
build/limitSet es > /dev/null &
build/limitSet ea > /dev/null &
build/limitSet oa > /dev/null &
build/limitSet oe > /dev/null &
build/limitSet eb > /dev/null &
build/limitSet ob > /dev/null &
build/limitSet os > /dev/null &
wait

build/drawLim

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED