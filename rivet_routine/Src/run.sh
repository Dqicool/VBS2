#!/bin/bash
HEPMCHOME='/mnt/SSD/rivet_routine/'
DATA=/Data/
OUT=/Plots/

rm $HEPMCHOME/Build/Rivetqidong.so
rm $HEPMCHOME/$OUT/*.yoda
rm $HEPMCHOME/$OUT/*.root
rm $HEPMCHOME/log/*.log

docker run  --rm  -u `id -u $USER`:`id -g`  -v `pwd`:`pwd` -w `pwd`  hepstore/rivet rivet-buildplugin $HEPMCHOME/Build/Rivetqidong.so $HEPMCHOME/Src/qidong.cc

EXE="docker run  --rm  -u `id -u $USER`:`id -g` -v `pwd`:`pwd` -w `pwd` hepstore/rivet rivet -q --analysis-path $HEPMCHOME/Build -a qidong"

mkdir -p $HEPMCHOME$OUT
rm -rf $HEPMCHOME$OUT/*.root
rm -rf $HEPMCHOME$OUT/*.yoda
for entry in `ls $HEPMCHOME$DATA`; do
    $EXE $HEPMCHOME$DATA$entry -H $HEPMCHOME$OUT$entry.yoda > $HEPMCHOME/log/$entry.log &
done
wait

# cd $HEPMCHOME$OUT/
# #echo "docker run -it --rm -u `id -u $USER`:`id -g`  -v `pwd`:`pwd` -w `pwd`  hepstore/rivet yodamerge -o $HEPMCHOME$OUT/result.yoda ./*"
#docker run --rm -u `id -u $USER`:`id -g`  -v `pwd`:`pwd` -w `pwd`  hepstore/rivet yodamerge --add `ls $HEPMCHOME$OUT` -o $HEPMCHOME/result.yoda
for entry in `ls $HEPMCHOME$OUT`; do
    yoda2root $HEPMCHOME$OUT$entry $HEPMCHOME$OUT/root_$entry.root 
done
