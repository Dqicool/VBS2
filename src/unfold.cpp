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
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <Math/GenVector/PtEtaPhiM4Dfwd.h>
#include <Math/Vector4Dfwd.h>
#include<Math/Vector3Dfwd.h>
#include<Math/GenVector/Boost.h>
#include <TChain.h>
#include<TROOT.h>
#include<TStyle.h>
#include<TLatex.h>

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/EulerAngles.h"
#include "Math/GenVector/AxisAngle.h"
#include "Math/GenVector/Quaternion.h"
#include "Math/GenVector/LorentzRotation.h"
#include "Math/GenVector/Boost.h"
#include "Math/GenVector/Transform3D.h"
#include "Math/GenVector/Plane3D.h"
#include "Math/GenVector/VectorUtil.h"

#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>
#include <RooUnfoldErrors.h>
#include <cstdlib>
#include <TLegend.h>
#include<THashList.h>
#include"libs/ErrorProp.h"

//#define UNFOLDSELF
#define UNFLOLDDATA
using namespace std;

TH1* devideBinW(TH1* h){
    uint nbins = h->GetNbinsX();
    h->SetBinContent(0,0);
    h->SetBinContent(nbins+1,0);
    h->SetBinError(0,0);
    h->SetBinError(nbins+1,0);
    for(uint i = 1;i <= nbins;i++){
        double ibinw = h->GetBinWidth(i);
        double ibinc = h->GetBinContent(i) / ibinw;
        double ibine = h->GetBinError(i) / ibinw;
        h->SetBinContent(i, ibinc);
        h->SetBinError(i, ibine);
    }
    return (TH1*)h->Clone();
}

void distUnfold(const char* dist, const char* labeltex, TH1D* h_data, const char* outf, double ymin, double ymax)
{   

    TH1::SetDefaultSumw2();
    gStyle->SetPaintTextFormat("1.2f");
    //load file name into a vector
        std::vector<string> filevec{};
        filevec.push_back("output/resp_out/364364.Sherpa_222_NNPDF30NNLO_lllljj_EW6_noHiggs.root");
        filevec.push_back("output/resp_out/346228.PowhegPy8EG_NNPDF30_AZNLOCTEQ6L1_VBFH125_ZZ4lep_notau.root");

        filevec.push_back("output/resp_out/345060.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_ZZ4l.root");
        filevec.push_back("output/resp_out/345706.Sherpa_222_NNPDF30NNLO_ggllll_130M4l.root");
        //filevec.push_back("output/resp_out/345708.Sherpa_222_NNPDF30NNLO_ggllllNoHiggs_0M4l130.root");

        filevec.push_back("output/resp_out/364250.Sherpa_222_NNPDF30NNLO_llll.root");
        filevec.push_back("output/resp_out/346340.PowhegPy8EG_A14NNPDF23_NNPDF30ME_ttH125_ZZ4l_allhad.root");
        filevec.push_back("output/resp_out/346341.PowhegPy8EG_A14NNPDF23_NNPDF30ME_ttH125_ZZ4l_semilep.root");
        filevec.push_back("output/resp_out/346342.PowhegPy8EG_A14NNPDF23_NNPDF30ME_ttH125_ZZ4l_dilep.root");

    //getting calculated response and stack
        std::vector<RooUnfoldResponse> resp_use_vec{};
        std::vector<RooUnfoldResponse> resp_disp_vec{};
        for (uint i = 0; i< filevec.size(); i++)
        {   
            string disp_save_name = "disp_resp_" + (string)(dist);
            string use_save_name = "use_resp_" + (string)(dist);
            TFile* f = TFile::Open(&(filevec[i])[0]);
            resp_use_vec.push_back(*((RooUnfoldResponse *)f->Get(&use_save_name[0])));
            resp_disp_vec.push_back(*((RooUnfoldResponse *)f->Get(&disp_save_name[0])));
        }
        RooUnfoldResponse resp_disp = resp_disp_vec[0];
        RooUnfoldResponse resp_use = resp_use_vec[0];
        for(uint i = 1; i < resp_disp_vec.size(); i++)
        {
            resp_disp.Add(resp_disp_vec[i]);
            resp_use.Add(resp_use_vec[i]);
        }

    //unfold
        auto h_meas = resp_use.Hmeasured();
        auto h_true = resp_use.Htruth();
        auto h_fake = resp_use.Hfakes();

        //RooUnfoldBayes unfoldself(&resp_use, h_meas, 2);
        RooUnfoldBayes unfolddata(&resp_use, h_data, 1);
        auto cov = unfolddata.Ereco(RooUnfold::kCovariance);
        TH1D* h_unfold= (TH1D*) unfolddata.Hreco(RooUnfold::kCovariance);
        unfolddata.PrintTable (std::cout);

    //save 
        string sn_reco = (string)"h_" + (string)dist + (string)"_reco";
        string sn_true = (string)"h_" + (string)dist + (string)"_true";

        TFile* out = TFile::Open(outf, "recreate");
        
        auto h_cov = new TH2D(cov);
        h_unfold->Clone(sn_reco.data())->Write();
        h_true->Clone(sn_true.data())->Write();
        h_cov->Clone("h_Cov")->Write();
        string pic_name = "plots/unfold/" + (string)(const char*)dist;
        auto unfold_name = pic_name + "_unfold.png";
        auto resp_name = pic_name + "_resp.png";
        auto resph_name = pic_name + "_resph.png";
        out->Close();
    //plot
        auto* c1 = new TCanvas("c1","",1200,1200);
                //h_unfold->SetAxisRange(0,25,"Y");
                h_unfold->SetTitle("");
                h_unfold->SetLineColorAlpha(kBlack,0);
                h_unfold->SetMarkerColor(kBlack);
                h_unfold->SetFillColorAlpha(kBlack, 0.2);
                h_unfold->SetMarkerStyle(kFullCircle);
                h_unfold->SetMarkerSize(2.0);
                h_unfold->SetLineWidth(2.0);
                h_unfold->GetXaxis()->SetTitle(labeltex);
                h_unfold->GetYaxis()->SetTitle("Events");
                h_unfold->SetAxisRange(ymin,ymax,"Y");
                h_unfold->SetStats(0);
                auto h_diff_unf = devideBinW(h_unfold);
                h_diff_unf->Draw("E1");
                h_diff_unf->Draw("E2 same");

                h_data->SetLineColorAlpha(kBlack,0);
                h_data->SetMarkerColor(kBlue);
                h_data->SetFillColorAlpha(kBlue, 0.1);
                h_data->SetMarkerStyle(kFullSquare);
                h_data->SetMarkerSize(2.0);
                h_data->SetLineWidth(2.0);
                h_data->SetStats(0);
                auto h_diff_data = devideBinW(h_data);
                h_diff_data->Draw("same E1");
                h_diff_data->Draw("E2 same");

                h_meas->SetLineColor(kBlue);
                h_meas->SetLineWidth(2);
                h_meas->Draw("SAME HIST");
                h_meas->SetStats(0);
                auto h_diff_meas = devideBinW(h_meas);
                h_diff_meas->SetFillStyle(3004);
                h_diff_meas->SetFillColorAlpha(kBlue,1);
                
                h_diff_meas->Draw("SAME E2");

                h_true->SetLineColor(kBlack);
                h_true->SetLineWidth(2);
                h_true->Draw("SAME HIST");
                h_true->SetStats(0);
                auto h_diff_true = devideBinW(h_true);
                h_diff_true->SetFillStyle(3004);
                h_diff_true->SetFillColorAlpha(kBlack,1);
                h_diff_true->Draw("SAME E2");
                

                TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
                legend->AddEntry(h_diff_unf,"Unfold","lepf");
                legend->AddEntry(h_diff_true,"Truth Pred","lepf");
                legend->AddEntry(h_diff_data, "Data","lepf");
                legend->AddEntry(h_diff_meas, "Detec Pred","lepf");
                legend->Draw();
            c1->SaveAs((&unfold_name[0]));

        auto* c2 = new TCanvas("c2","",2000,2000);
            c2->SetBottomMargin(0.2);
            c2->SetLeftMargin(0.2);
            auto m_resp = resp_disp.Mresponse();
            auto tmp = new TH2D(m_resp);
            for(int i=0; i < h_meas->GetNbinsX(); i++){
                std::stringstream buffer;
                buffer <<std::setprecision(3)<< h_meas->GetBinLowEdge(i+1);
                tmp->GetXaxis()->ChangeLabel(i+1,-1,-1,-1,-1,-1, &(buffer.str())[0]);
                tmp->GetYaxis()->ChangeLabel(i+1,-1,-1,-1,-1,-1, &(buffer.str())[0]);
            }
            tmp->SetTitle("Colums sum to unity, indecates how truth smeared");
            tmp->GetXaxis()->SetTitle("True");
            tmp->GetYaxis()->SetTitle("Meas");
            tmp->SetStats(0);
            tmp->SetMaximum(1);
            tmp->SetMinimum(0);
            tmp->Draw("COLZ");
            tmp->Draw("TEXT SAME");
            c2->SaveAs(&resp_name[0]);

        auto* c3 = new TCanvas("c3","",2000,2000);
            auto h_resp = resp_disp.Hresponse();
            h_resp->SetTitle("Colums sum to unity, indecates how measurments distributed in truth");
            h_resp->GetXaxis()->SetTitle("Meas");
            h_resp->GetYaxis()->SetTitle("True");
            h_resp->SetStats(0);
            h_resp->Draw("COLZ");
            h_resp->Draw("TEXT SAME");
            auto nx = h_resp->GetNbinsX();
            auto ny = h_resp->GetNbinsY();
            for(int i = 1; i <= nx; i++){
                auto linetotal = h_resp->Integral(i, i, 1, ny);
                for(int j = 1; j <= ny; j++){
                    h_resp->SetBinContent(i,j, h_resp->GetBinContent(i,j)/linetotal);
                }
            }
            h_resp->SetMaximum(1);
            h_resp->SetMinimum(0);
            //h_resp->Draw("COLZ");
            //h_resp->Draw("TEXT SAME");
            c3->SaveAs(&resph_name[0]);

}


int main(){
    TFile* inphijj = TFile::Open("output/rebin_out/dPhiJJ.root","read");
    auto datadelphijj  = (TH1D*)inphijj->Get("h_jj_dphi_pass_cutdat");
    distUnfold("jj_dphi", "#Delta#phi_{jj} [rad]", datadelphijj, "output/unfold_out/dphijj.root", 0, 0.6);

    TFile* inat1 = TFile::Open("output/rebin_out/at1.root","read");
    auto dataat1 = (TH1D*)inat1->Get("h_Angle_theta1_pass_cutdat");
    distUnfold("Angle_theta1", "#theta 1 [rad]", dataat1,"output/unfold_out/at1.root", 0, 4);

    TFile* inat2 = TFile::Open("output/rebin_out/at2.root","read");
    auto dataat2 = (TH1D*)inat2->Get("h_Angle_theta2_pass_cutdat");
    distUnfold("Angle_theta2", "#theta 2 [rad]", dataat2, "output/unfold_out/at2.root", 0, 4);

    TFile* inphi0 = TFile::Open("output/rebin_out/ap0.root","read");
    auto dataphi0 = (TH1D*)inphi0->Get("h_Angle_phi0_pass_cutdat");
    distUnfold("Angle_phi0", "#phi 0", dataphi0, "output/unfold_out/ap0.root", 0, 0.6);
}

