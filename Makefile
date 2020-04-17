PWD = $(shell pwd)
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = -L/usr/lib/root -lGenVector -lEve -lEG -lGeom -lGed -lRGL -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMathMore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic
ROOTINC   = -I/mnt/SSD/VBS2/

CXX = c++

UNFOLD_COMP_FLAG = -g -Wall -rdynamic -fPIC -fPIE -pthread -std=c++17 -m64 -DHAVE_TSVDUNFOLD=1 -DMAKEBUILD 
UNFOLD_COMP_INCL = -I/mnt/SSD/VBS2/RooUnfold/src/ -I/mnt/SSD/VBS2/RooUnfold -I/mnt/SSD/VBS2

UNFOLD_LINK_FLAG =  -g -m64 -rdynamic
UNFOLD_LINK_LIBS =  -L/mnt/SSD/VBS2/RooUnfold/ -lRooUnfold_static  -L/usr/lib/root -lGenVector -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMathMore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic -lUnfold -lRooFit -lRooFitCore -lThread -lMinuit -lFoam -lHtml

ana: src/analyse.cpp libs/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/analyse.cpp -o build/analyse $(ROOTLIBS) $(ROOTINC)

ananew: src/ana_new.cpp libs/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/ana_new.cpp -o build/ana_new $(ROOTLIBS) $(ROOTINC)

draw: src/draw.cpp libs/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/draw.cpp  -o build/draw $(ROOTLIBS) $(ROOTINC)

stack: src/stack.cpp libs/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/stack.cpp  -o build/stack $(ROOTLIBS) $(ROOTINC)

cutflow: src/cutflow.cpp libs/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/cutflow.cpp  -o build/cutflow $(ROOTLIBS) $(ROOTINC)

stackcutflow: src/stackcutflow.cpp libs/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/stackcutflow.cpp -o build/stackcutflow $(ROOTLIBS) $(ROOTINC)

binSizeEva: src/binSizeEva.cpp libs/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/binSizeEva.cpp -o build/binSizeEva $(ROOTLIBS) $(ROOTINC)

calResp: src/calResp.cpp libs/genAna.h
	$(CXX) $(UNFOLD_COMP_FLAG) -c src/calResp.cpp -o build/calResp.o $(UNFOLD_COMP_INCL)
	$(CXX) $(UNFOLD_LINK_FLAG) -g  build/calResp.o -o build/calResp $(UNFOLD_LINK_LIBS)
	rm -rf build/calResp.o

unfold: src/unfold.cpp libs/genAna.h
	$(CXX) $(UNFOLD_COMP_FLAG) -g -c src/unfold.cpp -o build/unfold.o $(UNFOLD_COMP_INCL)
	$(CXX) $(UNFOLD_LINK_FLAG) -g build/unfold.o -o build/unfold $(UNFOLD_LINK_LIBS)
	rm -rf build/unfold.o

SMvali: src/SMvali.cpp libs/genAna.h
	$(CXX) $(ROOTFLAGS) -g src/SMvali.cpp -o build/SMvali $(ROOTLIBS) $(ROOTINC)

NPvSM: src/NPvSM.cpp 
	$(CXX) $(ROOTFLAGS) -g src/NPvSM.cpp -o build/NPvSM $(ROOTLIBS) $(ROOTINC)

limitSet: src/limitSet.cpp
	$(CXX) $(ROOTFLAGS) -g src/limitSet.cpp -o build/limitSet $(ROOTLIBS) $(ROOTINC)

test: src/test.cpp
	$(CXX) $(ROOTFLAGS) -g src/test.cpp -o build/test $(ROOTLIBS) $(ROOTINC)