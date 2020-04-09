#!/bin/bash
make binSizeEva &
make calResp &
make unfold &
wait

build/binSizeEva

    VBSHOME='/mnt/SSD/VBS2/'
    EXE='/build/calResp'
    DATA='/output/analyse_out/'
    OUT='/output/resp_out/'

    rm -rf $OUT/*
    mkdir -p $VBSHOME$OUT
    i=1
    for entry in `ls $VBSHOME$DATA`; do
        echo [$i] $entry
        ($VBSHOME$EXE $VBSHOME$DATA$entry $VBSHOME$OUT$entry) &
        let i+=1
    done
    wait

build/unfold
