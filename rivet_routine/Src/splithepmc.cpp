#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#define FILEHEAD1 "HepMC::Version 2.06.09"
#define FILEHEAD2 "HepMC::IO_GenEvent-START_EVENT_LISTING"
#define FILEEND   "HepMC::IO_GenEvent-END_EVENT_LISTING"
// #define DEBUG
#ifndef DEBUG
int main(int argc, char** argv){
    std::ifstream infile(argv[1]);
    std::string line;
    std::string out_file_name(argv[1]);
    std::ofstream outfile((std::string)argv[2] + ".1");
    size_t nline = 0, n_evt = 0, n_file = 1;
    while (std::getline(infile, line))
    {  
        nline++;
        if(line.data()[0] == 'E'){
            n_evt++;
            if ((n_evt - 1) % 10000 == 0 && n_evt > 1){
                n_file++;
                outfile<<FILEEND<<std::endl;
                std::stringstream tmp;
                tmp<<argv[2]<<"."<<n_file;
                outfile.close();
                outfile.open(tmp.str());
                outfile<<FILEHEAD1<<std::endl;
                outfile<<FILEHEAD2<<std::endl;
            }
        }
        outfile<<line<<std::endl;
    }
}
#else
int main(int argc, char** argv){
    std::ifstream infile("./hepmc/Data/SM_EW_1E4.hepmc");
    std::string line;
    std::string out_file_name("./hepmc/Data/SM_EW_1E4.hepmc");
    std::ofstream outfile((std::string)"./hepmc/Data/SM_EW_1E4.hepmc" + ".1");
    size_t nline = 0, n_evt = 0, n_file = 1;
    while (std::getline(infile, line))
    {  
        nline++;
        if(line.data()[0] == 'E'){
            n_evt++;
            if (n_evt % 1000 == 0){
                n_file++;
                outfile<<FILEEND<<std::endl;
                std::stringstream tmp;
                tmp<<"./hepmc/Data/SM_EW_1E4.hepmc."<<n_file;
                outfile.close();
                outfile.open(tmp.str());
                outfile<<FILEHEAD1<<std::endl;
                outfile<<FILEHEAD2<<std::endl;
            }
        }
        outfile<<line<<std::endl;
    }
}
#endif

