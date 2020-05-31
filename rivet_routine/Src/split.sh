#!/bin/bash
c++ ./Src/splithepmc.cpp -o Build/splitHepmc
HEPMCHOME='/mnt/SSD/VBS/hepmc/'
EXE='/Build/splitHepmc'
DATA='/Data/'
OUT='/Split/'
rm -rf $HEPMCHOME/$OUT/*
for entry in `ls $HEPMCHOME$DATA`; do
    $HEPMCHOME$EXE $HEPMCHOME$DATA$entry $HEPMCHOME$OUT$entry &
done
wait
for entry in `ls $HEPMCHOME$OUT`; do
    gzip $HEPMCHOME$OUT$entry &
done
wait

    