#include "libs/genAna.h"
#include <THStack.h>

#define XSPERBIN 14./139.

using namespace std;

std::vector<double> mulScalar(std::vector<double> mat, double x){
    std::vector<double> ret{};
    for(double i = 0; i < mat.size(); i++){
        ret.push_back(mat[i]*x);
    }
    return ret;
}

std::vector<double> pluScalar(std::vector<double> mat, double x){
    std::vector<double> ret{};
    for(double i = 0; i < mat.size(); i++){
        ret.push_back(mat[i] + x);
    }
    return ret;
}

void binSizeEva(const char* fn_sig,      const char* hn_sig, 
                const char* fn_sigtrue,  const char* hn_sigtrue,
                const char* fn_bkg,      const char* hn_bkg,
                const char* fn_dat,      const char* hn_dat, 

                const char* fn_CW_NP,    const char* hn_CW_NP,
                const char* fn_CW_NP2,   const char* hn_CW_NP2,
                const char* fn_CWtil_NP, const char* hn_CWtil_NP,
                const char* fn_CWtil_NP2,const char* hn_CWtil_NP2,

                const char* fn_CHW_NP,    const char* hn_CHW_NP,
                const char* fn_CHW_NP2,   const char* hn_CHW_NP2,
                const char* fn_CHWtil_NP, const char* hn_CHWtil_NP,
                const char* fn_CHWtil_NP2,const char* hn_CHWtil_NP2,

                const char* fn_out,   
                bool sym_cent){

    TH1::SetDefaultSumw2();
    TFile* f_sig = TFile::Open(fn_sig, "read");
    TFile* f_sigtrue = TFile::Open(fn_sigtrue, "read");
    TFile* f_bkg = TFile::Open(fn_bkg, "read");

    TFile* f_CW_NP = TFile::Open(fn_CW_NP, "read");
    TFile* f_CW_NP2 = TFile::Open(fn_CW_NP2, "read");
    TFile* f_CWtil_NP = TFile::Open(fn_CWtil_NP, "read");
    TFile* f_CWtil_NP2 = TFile::Open(fn_CWtil_NP2, "read");

    TFile* f_CHW_NP = TFile::Open(fn_CHW_NP, "read");
    TFile* f_CHW_NP2 = TFile::Open(fn_CHW_NP2, "read");
    TFile* f_CHWtil_NP = TFile::Open(fn_CHWtil_NP, "read");
    TFile* f_CHWtil_NP2 = TFile::Open(fn_CHWtil_NP2, "read");


    string sn_sig = (string)hn_sig + "sig";
    string sn_bkg = (string)hn_bkg + "bkg";

    string sn_CW_NP = (string)hn_sig + "_CW_NP";
    string sn_CW_NP2 = (string)hn_sig + "_CW_NP2";
    string sn_CWtil_NP = (string)hn_sig + "_CWtil_NP";
    string sn_CWtil_NP2 = (string)hn_sig + "_CWtil_NP2";

    string sn_CHW_NP = (string)hn_sig + "_CHW_NP";
    string sn_CHW_NP2 = (string)hn_sig + "_CHW_NP2";
    string sn_CHWtil_NP = (string)hn_sig + "_CHWtil_NP";
    string sn_CHWtil_NP2 = (string)hn_sig + "_CHWtil_NP2";


    TH1D * h_sig = (TH1D*)f_sig->Get(hn_sig);
    TH1D * h_sigtrue = (TH1D*)f_sigtrue->Get(hn_sigtrue);
    TH1D * h_bkg = (TH1D*)f_bkg->Get(hn_bkg);

    TH1D * h_CW_NP = (TH1D*)f_CW_NP->Get(hn_CW_NP);
    TH1D * h_CW_NP2 = (TH1D*)f_CW_NP2->Get(hn_CW_NP2);
    TH1D * h_CWtil_NP = (TH1D*)f_CWtil_NP->Get(hn_CWtil_NP);
    TH1D * h_CWtil_NP2 = (TH1D*)f_CWtil_NP2->Get(hn_CWtil_NP2);

    TH1D * h_CHW_NP = (TH1D*)f_CHW_NP->Get(hn_CHW_NP);
    TH1D * h_CHW_NP2 = (TH1D*)f_CHW_NP2->Get(hn_CHW_NP2);
    TH1D * h_CHWtil_NP = (TH1D*)f_CHWtil_NP->Get(hn_CHWtil_NP);
    TH1D * h_CHWtil_NP2 = (TH1D*)f_CHWtil_NP2->Get(hn_CHWtil_NP2);
    

    auto nbin = h_sig->GetNbinsX();
    auto wbin = h_sig->GetBinWidth(1);
    auto xmin = h_sig->GetBinLowEdge(1);

    std::vector<double> sepa;
    if(sym_cent)    sepa = {(double)nbin / 2};
    else            sepa = {0};
    for(int i=sepa[0]+1; i<nbin; i++){
        int dam = sepa[sepa.size() - 1];
        double damn = h_sig->Integral(dam,i);
        if(damn > XSPERBIN) {sepa.push_back(i);}
    }
    if (*(sepa.end()) != nbin) sepa.push_back(nbin);
    double tmp1 = h_sig->Integral(sepa[sepa.size()-2], nbin);
    double tmp2 = (XSPERBIN * 3. / 5.);
    if(tmp1 < tmp2)
        {sepa.erase(sepa.end() - 2);}
    std::vector<double> sepaL=sepa;
    if(sym_cent){
        std::reverse(sepaL.begin(), sepaL.end());
        sepaL = pluScalar(mulScalar(sepaL, -1), nbin);
        sepa.insert(sepa.begin(), sepaL.begin(), sepaL.end()-1);
    }
    sepa = mulScalar(sepa, wbin);
    sepa = pluScalar(sepa, xmin);


    auto h_rebin_sig = h_sig->Rebin(sepa.size()-1, &sn_sig[0], &sepa[0]);
    auto h_rebin_bkg = h_bkg->Rebin(sepa.size()-1, &sn_bkg[0], &sepa[0]);
    auto h_rebin_sigtrue = h_sigtrue->Rebin(sepa.size()-1, hn_sigtrue, &sepa[0]);

    auto h_rebin_CW_NP   = h_CW_NP->Rebin(sepa.size()-1, &sn_CW_NP[0], &sepa[0]);
    auto h_rebin_CW_NP2   = h_CW_NP2->Rebin(sepa.size()-1, &sn_CW_NP2[0], &sepa[0]);
    auto h_rebin_CWtil_NP   = h_CWtil_NP->Rebin(sepa.size()-1, &sn_CWtil_NP[0], &sepa[0]);
    auto h_rebin_CWtil_NP2   = h_CWtil_NP2->Rebin(sepa.size()-1, &sn_CWtil_NP2[0], &sepa[0]);
    
    auto h_rebin_CHW_NP   = h_CHW_NP->Rebin(sepa.size()-1, &sn_CHW_NP[0], &sepa[0]);
    auto h_rebin_CHW_NP2   = h_CHW_NP2->Rebin(sepa.size()-1, &sn_CHW_NP2[0], &sepa[0]);
    auto h_rebin_CHWtil_NP   = h_CHWtil_NP->Rebin(sepa.size()-1, &sn_CHWtil_NP[0], &sepa[0]);
    auto h_rebin_CHWtil_NP2   = h_CHWtil_NP2->Rebin(sepa.size()-1, &sn_CHWtil_NP2[0], &sepa[0]);


    string sn_dat = (string)hn_dat + "dat";
    TFile* f_dat = TFile::Open(fn_dat, "read");
    TH1D * h_dat = (TH1D*)f_dat->Get(hn_dat);
    auto h_rebin_dat = h_dat->Rebin(sepa.size()-1, &sn_dat[0], &sepa[0]);
    h_rebin_dat->Add(h_rebin_bkg, -1);
    TFile* outfile = TFile::Open(fn_out, "recreate");
    
    h_rebin_sig->Write();
    h_rebin_bkg->Write();
    h_rebin_sigtrue->Write();
    h_rebin_dat->Write();
    
    h_rebin_CW_NP->Write();
    h_rebin_CW_NP2->Write();
    h_rebin_CWtil_NP->Write();
    h_rebin_CWtil_NP2->Write();

    h_rebin_CHW_NP->Write();
    h_rebin_CHW_NP2->Write();
    h_rebin_CHWtil_NP->Write();
    h_rebin_CHWtil_NP2->Write();
    outfile->Close();
}

int main(){

    binSizeEva( "output/stack_out/sigs/jjDelPhi.root",  "jj_dphi_cut_h",
                "output/stack_out/sigs/jj_true_delphi.root",            "jj_truthBorn_dphi_cut_h",
                "output/stack_out/bkgs/jjDelPhi.root",                  "jj_dphi_cut_h",
                "output/stack_out/Data.1516.root",                      "h_jj_dphi_pass_cut",

                "output/rivet_out/root_CW_NP.hepmc.gz.yoda.root",       "qidong/jjDphi",
                "output/rivet_out/root_CW_NP2.hepmc.gz.yoda.root",      "qidong/jjDphi",
                "output/rivet_out/root_CWtil_NP.hepmc.gz.yoda.root",    "qidong/jjDphi",
                "output/rivet_out/root_CWtil_NP2.hepmc.gz.yoda.root",   "qidong/jjDphi",

                "output/rivet_out/root_CHW_NP.hepmc.gz.yoda.root",      "qidong/jjDphi",
                "output/rivet_out/root_CHW_NP2.hepmc.gz.yoda.root",     "qidong/jjDphi",
                "output/rivet_out/root_CHWTIL_NP.hepmc.gz.yoda.root",   "qidong/jjDphi",
                "output/rivet_out/root_CHWTIL_NP2.hepmc.gz.yoda.root",  "qidong/jjDphi",

                "output/rebin_out/dPhiJJ.root", 
                1);
    
    binSizeEva("output/stack_out/sigs/at1.root",        "Angle_theta1_cut_h",
                "output/stack_out/sigs/at1_true.root",  "Angle_truthBorn_theta1_cut_h",
                "output/stack_out/bkgs/at1.root",       "Angle_theta1_cut_h",
                "output/stack_out/Data.1516.root",      "h_Angle_theta1_pass_cut",
                
                "output/rivet_out/root_CW_NP.hepmc.gz.yoda.root",       "qidong/at1",
                "output/rivet_out/root_CW_NP2.hepmc.gz.yoda.root",      "qidong/at1",
                "output/rivet_out/root_CWtil_NP.hepmc.gz.yoda.root",    "qidong/at1",
                "output/rivet_out/root_CWtil_NP2.hepmc.gz.yoda.root",   "qidong/at1",

                "output/rivet_out/root_CHW_NP.hepmc.gz.yoda.root",      "qidong/at1",
                "output/rivet_out/root_CHW_NP2.hepmc.gz.yoda.root",     "qidong/at1",
                "output/rivet_out/root_CHWTIL_NP.hepmc.gz.yoda.root",   "qidong/at1",
                "output/rivet_out/root_CHWTIL_NP2.hepmc.gz.yoda.root",  "qidong/at1",

                "output/rebin_out/at1.root", 
                0);
    
    binSizeEva( "output/stack_out/sigs/at2.root",                        "Angle_theta2_cut_h", 
                "output/stack_out/sigs/at2_true.root",                   "Angle_truthBorn_theta2_cut_h", 
                "output/stack_out/bkgs/at2.root",                        "Angle_theta2_cut_h", 
                "output/stack_out/Data.1516.root",                       "h_Angle_theta2_pass_cut", 
               
                "output/rivet_out/root_CW_NP.hepmc.gz.yoda.root",       "qidong/at2",
                "output/rivet_out/root_CW_NP2.hepmc.gz.yoda.root",      "qidong/at2",
                "output/rivet_out/root_CWtil_NP.hepmc.gz.yoda.root",    "qidong/at2",
                "output/rivet_out/root_CWtil_NP2.hepmc.gz.yoda.root",   "qidong/at2",

                "output/rivet_out/root_CHW_NP.hepmc.gz.yoda.root",      "qidong/at2",
                "output/rivet_out/root_CHW_NP2.hepmc.gz.yoda.root",     "qidong/at2",
                "output/rivet_out/root_CHWTIL_NP.hepmc.gz.yoda.root",   "qidong/at2",
                "output/rivet_out/root_CHWTIL_NP2.hepmc.gz.yoda.root",  "qidong/at2",
               
                "output/rebin_out/at2.root",
                0);

    binSizeEva( "output/stack_out/sigs/ats.root",                        "Angle_theta_star_cut_h", 
                "output/stack_out/sigs/ats_true.root",                   "Angle_truthBorn_theta_star_cut_h", 
                "output/stack_out/bkgs/ats.root",                        "Angle_theta_star_cut_h", 
                "output/stack_out/Data.1516.root",                       "h_Angle_theta_star_pass_cut", 
               
                "output/rivet_out/root_CW_NP.hepmc.gz.yoda.root",       "qidong/ats",
                "output/rivet_out/root_CW_NP2.hepmc.gz.yoda.root",      "qidong/ats",
                "output/rivet_out/root_CWtil_NP.hepmc.gz.yoda.root",    "qidong/ats",
                "output/rivet_out/root_CWtil_NP2.hepmc.gz.yoda.root",   "qidong/ats",

                "output/rivet_out/root_CHW_NP.hepmc.gz.yoda.root",      "qidong/ats",
                "output/rivet_out/root_CHW_NP2.hepmc.gz.yoda.root",     "qidong/ats",
                "output/rivet_out/root_CHWTIL_NP.hepmc.gz.yoda.root",   "qidong/ats",
                "output/rivet_out/root_CHWTIL_NP2.hepmc.gz.yoda.root",  "qidong/ats",
               
                "output/rebin_out/ats.root",
                0);
    binSizeEva( "output/stack_out/sigs/aps.root",                        "Angle_phi_star_cut_h", 
                "output/stack_out/sigs/aps_true.root",                   "Angle_truthBorn_phi_star_cut_h", 
                "output/stack_out/bkgs/aps.root",                        "Angle_phi_star_cut_h", 
                "output/stack_out/Data.1516.root",                       "h_Angle_phi_star_pass_cut", 
               
                "output/rivet_out/root_CW_NP.hepmc.gz.yoda.root",       "qidong/aps",
                "output/rivet_out/root_CW_NP2.hepmc.gz.yoda.root",      "qidong/aps",
                "output/rivet_out/root_CWtil_NP.hepmc.gz.yoda.root",    "qidong/aps",
                "output/rivet_out/root_CWtil_NP2.hepmc.gz.yoda.root",   "qidong/aps",

                "output/rivet_out/root_CHW_NP.hepmc.gz.yoda.root",      "qidong/aps",
                "output/rivet_out/root_CHW_NP2.hepmc.gz.yoda.root",     "qidong/aps",
                "output/rivet_out/root_CHWTIL_NP.hepmc.gz.yoda.root",   "qidong/aps",
                "output/rivet_out/root_CHWTIL_NP2.hepmc.gz.yoda.root",  "qidong/aps",
               
                "output/rebin_out/aps.root",
                1);
    binSizeEva( "output/stack_out/sigs/ap0.root",                        "Angle_phi0_cut_h", 
                "output/stack_out/sigs/ap0_true.root",                   "Angle_truthBorn_phi0_cut_h", 
                "output/stack_out/bkgs/ap0.root",                        "Angle_phi0_cut_h", 
                "output/stack_out/Data.1516.root",                       "h_Angle_phi0_pass_cut", 
               
                "output/rivet_out/root_CW_NP.hepmc.gz.yoda.root",       "qidong/ap0",
                "output/rivet_out/root_CW_NP2.hepmc.gz.yoda.root",      "qidong/ap0",
                "output/rivet_out/root_CWtil_NP.hepmc.gz.yoda.root",    "qidong/ap0",
                "output/rivet_out/root_CWtil_NP2.hepmc.gz.yoda.root",   "qidong/ap0",

                "output/rivet_out/root_CHW_NP.hepmc.gz.yoda.root",      "qidong/ap0",
                "output/rivet_out/root_CHW_NP2.hepmc.gz.yoda.root",     "qidong/ap0",
                "output/rivet_out/root_CHWTIL_NP.hepmc.gz.yoda.root",   "qidong/ap0",
                "output/rivet_out/root_CHWTIL_NP2.hepmc.gz.yoda.root",  "qidong/ap0",
               
                "output/rebin_out/ap0.root",
                1);
    
    binSizeEva( "output/stack_out/sigs/ap1.root",                        "Angle_phi1_cut_h", 
                "output/stack_out/sigs/ap1_true.root",                   "Angle_truthBorn_phi1_cut_h", 
                "output/stack_out/bkgs/ap1.root",                        "Angle_phi1_cut_h", 
                "output/stack_out/Data.1516.root",                       "h_Angle_phi1_pass_cut", 
               
                "output/rivet_out/root_CW_NP.hepmc.gz.yoda.root",       "qidong/ap1",
                "output/rivet_out/root_CW_NP2.hepmc.gz.yoda.root",      "qidong/ap1",
                "output/rivet_out/root_CWtil_NP.hepmc.gz.yoda.root",    "qidong/ap1",
                "output/rivet_out/root_CWtil_NP2.hepmc.gz.yoda.root",   "qidong/ap1",

                "output/rivet_out/root_CHW_NP.hepmc.gz.yoda.root",      "qidong/ap1",
                "output/rivet_out/root_CHW_NP2.hepmc.gz.yoda.root",     "qidong/ap1",
                "output/rivet_out/root_CHWTIL_NP.hepmc.gz.yoda.root",   "qidong/ap1",
                "output/rivet_out/root_CHWTIL_NP2.hepmc.gz.yoda.root",  "qidong/ap1",
               
                "output/rebin_out/ap1.root",
                1);
    /*binSizeEva( "output/stack_out/sigs/m4l.root",                        "llll_m_cut_h", 
                "output/stack_out/sigs/llll_true_m.root",                   "llll_truthBorn_m_cut_h", 
                "output/stack_out/bkgs/m4l.root",                        "llll_m_cut_h", 
                "output/stack_out/Data.1516.root",                       "llll_m_cut_h", 
               
                "output/rivet_out/root_CW_NP.hepmc.gz.yoda.root",       "qidong/m4l",
                "output/rivet_out/root_CW_NP2.hepmc.gz.yoda.root",      "qidong/m4l",
                "output/rivet_out/root_CWtil_NP.hepmc.gz.yoda.root",    "qidong/m4l",
                "output/rivet_out/root_CWtil_NP2.hepmc.gz.yoda.root",   "qidong/m4l",

                "output/rivet_out/root_CHW_NP.hepmc.gz.yoda.root",      "qidong/m4l",
                "output/rivet_out/root_CHW_NP2.hepmc.gz.yoda.root",     "qidong/m4l",
                "output/rivet_out/root_CHWTIL_NP.hepmc.gz.yoda.root",   "qidong/m4l",
                "output/rivet_out/root_CHWTIL_NP2.hepmc.gz.yoda.root",  "qidong/m4l",
               
                "output/rebin_out/m4l.root",
                0);
    binSizeEva( "output/stack_out/sigs/mjj.root",                        "jj_m_cut_h", 
                "output/stack_out/sigs/jj_true_m.root",                   "jj_truthBorn_m_cut_h", 
                "output/stack_out/bkgs/mjj.root",                        "jj_m_cut_h", 
                "output/stack_out/Data.1516.root",                       "jj_m_cut_h", 
               
                "output/rivet_out/root_CW_NP.hepmc.gz.yoda.root",       "qidong/mjj",
                "output/rivet_out/root_CW_NP2.hepmc.gz.yoda.root",      "qidong/mjj",
                "output/rivet_out/root_CWtil_NP.hepmc.gz.yoda.root",    "qidong/mjj",
                "output/rivet_out/root_CWtil_NP2.hepmc.gz.yoda.root",   "qidong/mjj",

                "output/rivet_out/root_CHW_NP.hepmc.gz.yoda.root",      "qidong/mjj",
                "output/rivet_out/root_CHW_NP2.hepmc.gz.yoda.root",     "qidong/mjj",
                "output/rivet_out/root_CHWTIL_NP.hepmc.gz.yoda.root",   "qidong/mjj",
                "output/rivet_out/root_CHWTIL_NP2.hepmc.gz.yoda.root",  "qidong/mjj",
               
                "output/rebin_out/mjj.root",
                0);*/
    return 0;
}