#!/bin/bash
make draw &
make stack &
wait
SECONDS=0
    EXE='/build/draw'
    DATA='/output/analyse_out/'
    OUT='/output/draw_out/'
    VBSHOME='/mnt/SSD/VBS2/'

    rm -rf $OUT/*
    mkdir -p $VBSHOME$OUT
    i=1
    for entry in `ls $VBSHOME$DATA`; do
        echo [$i] $entry
        ($VBSHOME$EXE $VBSHOME$DATA$entry $VBSHOME$OUT$entry) &
        let i+=1
    done
    # echo "[$i] data1516"
    # ($VBSHOME$EXE $VBSHOME/output/analyse_out/Data1516.root $VBSHOME/output/draw_out/Data1516.root) &
    wait
    cp $VBSHOME/$OUT/Data.1516.root $VBSHOME/output/stack_out/
    $VBSHOME/build/stack
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
