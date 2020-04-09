#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <RooUnfoldResponse.h>
#include <cstdlib>
#include <TLegend.h>
#include <TROOT.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include "vector"
#include <TMath.h>
#include <iostream>
#include <string.h>
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include<TROOT.h>
#include<TStyle.h>
#include<TLatex.h>
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#endif
//#define debug
#define PRWLUMI (36.21 + 58.45 + 44.31)

using namespace std;
void calResp(const char* infile, const char* out_resp,  const char* dist, const char* dist_true, TH1D* templ, bool recre)
{
    //load data file
        TH1::SetDefaultSumw2();
        TFile* f = TFile::Open(infile,"read");
        TTree* t = (TTree*)f->Get("SM4L_Nominal");
        RooUnfoldResponse resp(templ, templ, dist, infile);
        ROOT::RDataFrame df("SM4L_Nominal" ,infile);

    
    //declare variable pointer in tree
        float x_measured=0;
        float x_true = 0;
        float weight = 0;
        float weight_true = 0;
        bool pass_cut = 0;
        bool pass_true_cut = 0;
    
    if (df.HasColumn("fid_weight")){
        //setting pointers
            t->SetBranchAddress("pass_cut", &pass_cut);
            t->SetBranchAddress("pass_cut_truthBorn", &pass_true_cut);
            t->SetBranchAddress(dist,&x_measured);
            t->SetBranchAddress(dist_true,&x_true);
            t->SetBranchAddress("weight", &weight);
            t->SetBranchAddress("fid_weight", &weight_true);
            auto n_entry = t->GetEntries();
        
        // loop fills
            for( Long64_t i = 0; i<(n_entry); i++)
            {
                t->GetEntry(i);
                if(pass_cut && pass_true_cut){
                    resp.Fill(x_measured, x_true, weight / PRWLUMI);
                }
            }
        //save representation response
            const char* file_opt;
            if (recre)  file_opt = "recreate";
            else        file_opt = "update";
            TFile* out = TFile::Open(out_resp, file_opt);
            string disp_save_name = "disp_resp_" + (string)(dist);
            resp.SetName(&disp_save_name[0]);
            resp.Write();
        
        //loop over fake and miss
            for( Long64_t i = 0; i<(n_entry); i++)
            {
                t->GetEntry(i);
                if((pass_cut==0) && pass_true_cut){
                    resp.Miss(x_true, weight_true / PRWLUMI);
                }
                else if(pass_cut && (pass_true_cut==0)){
                    resp.Fake(x_measured, weight / PRWLUMI);
                }
            }
        //save respond
            string use_save_name = "use_resp_" + (string)(dist);
            resp.SetName(&use_save_name[0]);
            resp.Write();
            out->Close();
            f->Close();
    }
}
#ifndef debug
int main(int argc, char** argv){
    TFile* indelphijj  = TFile::Open("output/rebin_out/dPhiJJ.root","read");
    TFile* inat1       = TFile::Open("output/rebin_out/at1.root");
    TFile* inat2       = TFile::Open("output/rebin_out/at2.root");
    TFile* inap0       = TFile::Open("output/rebin_out/ap0.root");
    auto Templdelphijj = (TH1D*)indelphijj->Get("jj_dphi_cut_hsig");
    auto Templat1      = (TH1D*)inat1->Get("Angle_theta1_cut_hsig");
    auto Templat2      = (TH1D*)inat2->Get("Angle_theta2_cut_hsig");
    auto Templap0      = (TH1D*)inap0->Get("Angle_phi0_cut_hsig");
    const char* infile=argv[1];
    const char* out_resp = argv[2];
    calResp(infile, out_resp, "jj_dphi", "jj_truthBorn_dphi", Templdelphijj, 1);
    calResp(infile, out_resp, "Angle_theta1", "Angle_truthBorn_theta1", Templat1, 0);
    calResp(infile, out_resp, "Angle_theta2", "Angle_truthBorn_theta2", Templat2, 0);
    calResp(infile, out_resp, "Angle_phi0", "Angle_truthBorn_phi0", Templap0, 0);
}
#else
int main(int argc, char** argv){
    TFile* indelphijj = TFile::Open("output/rebin_out/dPhiJJ.root","read");
    auto Templdelphijj = (TH1D*)indelphijj->Get("jj_dphi_cut_hsig");
    const char* infile=argv[1];
    const char* out_resp = argv[2];
    calResp("output/analyse_out/364250.Sherpa_222_NNPDF30NNLO_llll.root", "debug.root", "jj_dphi", "jj_truthBorn_dphi", Templdelphijj, 1);
}
#endif