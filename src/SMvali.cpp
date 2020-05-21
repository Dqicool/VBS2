#include "libs/genAna.h"
#include "libs/ErrorProp.h"
#include<THStack.h>
#include<TLegend.h>
#include<cmath>
using namespace std;
TH1* getPlot(const char* file, const char* plot, double eft_c, Color_t col){
    TH1::SetDefaultSumw2();
    TFile* infile = TFile::Open(file, "read");
    TH1D * h = (TH1D*)infile->Get(plot);
    h->SetBins(h->GetNbinsX(), -M_PI, M_PI);
    h->SetFillColor(col);
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    auto h_rebin = h->Rebin(10, plot);
    h_rebin->Scale(eft_c / 5.0);
    return h_rebin;
}

TH1* getSMPlot(const char* file, const char* plot, Color_t col){
    TH1::SetDefaultSumw2();
    TFile* infile = TFile::Open(file, "read");
    TH1D * h = (TH1D*)infile->Get(plot);
    h->SetBins(h->GetNbinsX(), -M_PI, M_PI);
    h->SetFillColor(col);
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    auto h_rebin = h->Rebin(10, plot);
    return h_rebin;
}

std::vector<double> validate2z4l(){
    auto h_sm_2z_dpjj = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_SM_2z.hepmc.gz.yoda.root", "qidong/jjDphi", 5, kBlue);//5 for not scale SM
    auto h_sm_4l_dpjj = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_SM_4l.hepmc.gz.yoda.root", "qidong/jjDphi", 5, kRed);//5 for not scale SM
    h_sm_4l_dpjj = h_sm_4l_dpjj->Rebin(10,"");
    h_sm_2z_dpjj = h_sm_2z_dpjj->Rebin(10,"");
    h_sm_2z_dpjj->Scale(10./M_PI);
    h_sm_4l_dpjj->Scale(10./M_PI);
    TCanvas c1("c1","",2000,2000);
        c1.SetMargin(0.15,0.05,0.1,0.05);
        h_sm_2z_dpjj->SetMarkerColor(kBlue);
        h_sm_2z_dpjj->SetLineColor(kBlue);
        h_sm_2z_dpjj->SetFillColor(kBlue);
        h_sm_2z_dpjj->SetFillStyle(3017);

        h_sm_4l_dpjj->SetMarkerColor(kRed);
        h_sm_4l_dpjj->SetLineColor(kRed);
        h_sm_4l_dpjj->SetFillColor(kRed);
        h_sm_4l_dpjj->SetFillStyle(3017);

        h_sm_4l_dpjj->SetStats(0);
        h_sm_2z_dpjj->SetStats(0);

        h_sm_2z_dpjj->GetXaxis()->SetTitle("#Delta#phi_{jj}");
        h_sm_2z_dpjj->GetYaxis()->SetTitle("#frac{d#sigma}{d#Delta#phi_{jj}}");
        h_sm_2z_dpjj->Draw("E2");
        h_sm_4l_dpjj->Draw("E2,same");
        h_sm_2z_dpjj->Draw("E1 same");
        h_sm_4l_dpjj->Draw("E1 same");
        TLegend legend(0.6, 0.75, 0.89, 0.94);
                legend.AddEntry(h_sm_2z_dpjj,"SM EW ZZjj4l on-shell","lepf");
                legend.AddEntry(h_sm_4l_dpjj,"SM EW jj4l inclusive","lepf");
                legend.SetLineColorAlpha(kBlack,0);
                legend.SetTextFont(22);
                legend.Draw();
        auto cut_save_name = "plots/rivet/sm_validation.png";
    c1.SaveAs(&cut_save_name[0]);
    double corr,  err2z, err4l, err;
    double int4l = h_sm_4l_dpjj->IntegralAndError(1, 50, err4l);
    double int2z = h_sm_2z_dpjj->IntegralAndError(1, 50, err2z);
    corr = int4l / int2z;
    err  = ErrADiviB(corr, int4l, int2z, err4l, err2z, 0);
    return {corr, err};
}

void validatetheta(){
    
}

TH1* HistOverHist(TH1* up, TH1* down){
    for(int i=1; i <= up->GetNbinsX(); i++){
        double cont = up->GetBinContent(i)/down->GetBinContent(i);
        double err = ErrADiviB(cont, up->GetBinContent(i), down->GetBinContent(i), up->GetBinError(i), down->GetBinError(i), 0);
        up->SetBinContent(i, cont);
        up->SetBinError(i, err);
    }
    auto ret = (TH1*)(up->Clone());
    return ret;
}

int main(){
    auto corr = validate2z4l();
    // std::cout<<corr[0]<<"+-"<<corr[1]<<std::endl;
    // //get_histos from file
    //     //SM
    //     auto h_SM_dpj = getSMPlot("output/stack_out/sigs/jj_true_delphi.root", "jj_truthBorn_dphi_cut_h", kRed);
    //     auto h_SM_ats = getSMPlot("output/stack_out/sigs/ats_true.root", "Angle_truthBorn_theta_star_cut_h", kRed);
    //     auto h_SM_aps = getSMPlot("output/stack_out/sigs/aps_true.root", "Angle_truthBorn_phi_star_cut_h", kRed);
    //     auto h_SM_ap0 = getSMPlot("output/stack_out/sigs/ap0_true.root", "Angle_truthBorn_phi0_cut_h", kRed);
    //     auto h_SM_ap1 = getSMPlot("output/stack_out/sigs/ap1_true.root", "Angle_truthBorn_phi1_cut_h", kRed);
    //     auto h_SM_at1 = getSMPlot("output/stack_out/sigs/at1_true.root", "Angle_truthBorn_theta1_cut_h", kRed);
    //     auto h_SM_at2 = getSMPlot("output/stack_out/sigs/at2_true.root", "Angle_truthBorn_theta2_cut_h", kRed);
    //     //CW_NP
    //     auto h_CW_NP_dpj = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP.hepmc.gz.yoda.root", "qidong/jjDphi", 0.1, kBlue);
    //     auto h_CW_NP_ats = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP.hepmc.gz.yoda.root", "qidong/ats", 0.1, kBlue);
    //     auto h_CW_NP_aps = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP.hepmc.gz.yoda.root", "qidong/aps", 0.1, kBlue);
    //     auto h_CW_NP_ap0 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP.hepmc.gz.yoda.root", "qidong/ap0", 0.1, kBlue);
    //     auto h_CW_NP_ap1 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP.hepmc.gz.yoda.root", "qidong/ap1", 0.1, kBlue);
    //     auto h_CW_NP_at1 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP.hepmc.gz.yoda.root", "qidong/at1", 0.1, kBlue);
    //     auto h_CW_NP_at2 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP.hepmc.gz.yoda.root", "qidong/at2", 0.1, kBlue);
    //     //CW_NP2
    //     auto h_CW_NP2_dpj = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP2.hepmc.gz.yoda.root", "qidong/jjDphi", 0.1, kMagenta);
    //     auto h_CW_NP2_ats = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP2.hepmc.gz.yoda.root", "qidong/ats", 0.1, kMagenta);
    //     auto h_CW_NP2_aps = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP2.hepmc.gz.yoda.root", "qidong/aps", 0.1, kMagenta);
    //     auto h_CW_NP2_ap0 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP2.hepmc.gz.yoda.root", "qidong/ap0", 0.1, kMagenta);
    //     auto h_CW_NP2_ap1 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP2.hepmc.gz.yoda.root", "qidong/ap1", 0.1, kMagenta);
    //     auto h_CW_NP2_at1 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP2.hepmc.gz.yoda.root", "qidong/at1", 0.1, kMagenta);
    //     auto h_CW_NP2_at2 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CW_NP2.hepmc.gz.yoda.root", "qidong/at2", 0.1, kMagenta);
    //     //CWTIL_NP
    //     auto h_CWTIL_NP_dpj = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP.hepmc.gz.yoda.root", "qidong/jjDphi", 0.1, kBlue);
    //     auto h_CWTIL_NP_ats = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP.hepmc.gz.yoda.root", "qidong/ats", 0.1, kBlue);
    //     auto h_CWTIL_NP_aps = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP.hepmc.gz.yoda.root", "qidong/aps", 0.1, kBlue);
    //     auto h_CWTIL_NP_ap0 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP.hepmc.gz.yoda.root", "qidong/ap0", 0.1, kBlue);
    //     auto h_CWTIL_NP_ap1 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP.hepmc.gz.yoda.root", "qidong/ap1", 0.1, kBlue);
    //     auto h_CWTIL_NP_at1 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP.hepmc.gz.yoda.root", "qidong/at1", 0.1, kBlue);
    //     auto h_CWTIL_NP_at2 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP.hepmc.gz.yoda.root", "qidong/at2", 0.1, kBlue);
    //     //CWTIL_NP2
    //     auto h_CWTIL_NP2_dpj = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP2.hepmc.gz.yoda.root", "qidong/jjDphi", 0.1, kMagenta);
    //     auto h_CWTIL_NP2_ats = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP2.hepmc.gz.yoda.root", "qidong/ats", 0.1, kMagenta);
    //     auto h_CWTIL_NP2_aps = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP2.hepmc.gz.yoda.root", "qidong/aps", 0.1, kMagenta);
    //     auto h_CWTIL_NP2_ap0 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP2.hepmc.gz.yoda.root", "qidong/ap0", 0.1, kMagenta);
    //     auto h_CWTIL_NP2_ap1 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP2.hepmc.gz.yoda.root", "qidong/ap1", 0.1, kMagenta);
    //     auto h_CWTIL_NP2_at1 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP2.hepmc.gz.yoda.root", "qidong/at1", 0.1, kMagenta);
    //     auto h_CWTIL_NP2_at2 = getPlot("/mnt/MS/Projects/Qichen/rivet_routine/Plots/root_CWtil_NP2.hepmc.gz.yoda.root", "qidong/at2", 0.1, kMagenta);
    
    // //Draw and save
    //     TH1* h_tmp;
    //     //CW
    //         //dphijj
    //         TCanvas c00("c00","",2000,2000);
    //         h_CW_NP_dpj->SetAxisRange(-0.001,0.01, "Y");
    //         h_CW_NP_dpj->SetStats(0);
    //         h_CW_NP2_dpj->SetStats(0);
    //         (h_CW_NP_dpj)->Draw("E2"); 
    //         (h_CW_NP2_dpj)->Draw("E2 SAME");
    //         c00.SaveAs("plots/rivet/CW_NP_NP2_dpj.png");
            
    //         TCanvas c01("c01","",2000,2000);
    //         h_tmp = (TH1*)(h_SM_dpj->Clone());
    //         auto h_diff_cw = (TH1*)(h_CW_NP_dpj->Clone());
    //         h_diff_cw->Add(h_CW_NP2_dpj);
    //         h_diff_cw->Add(h_SM_dpj);
    //         h_diff_cw = HistOverHist(h_diff_cw, h_SM_dpj);
    //         h_diff_cw->SetAxisRange(1,1.2,"Y");
    //         h_diff_cw->Draw("E2");
    //         c01.SaveAs("plots/rivet/CW_SM_diff_dpj.png");
    //         //ap0
    //         TCanvas c10("c10","",2000,2000);
    //         h_CW_NP_ap0->SetStats(0);
    //         h_CW_NP2_ap0->SetStats(0);
    //         h_SM_ap0->SetStats(0);
    //         h_SM_ap0->SetAxisRange(-0.01,0.01,"Y");
    //         h_SM_ap0->Draw("E2");
    //         (h_CW_NP_ap0)->Draw("E2 SAME"); 
    //         (h_CW_NP2_ap0)->Draw("E2 SAME");
    //         c10.SaveAs("plots/rivet/CW_NP_NP2_ap0.png");
    //         //ap1
    //         TCanvas c3("c3","",2000,2000);
    //         h_CW_NP_ap1->SetStats(0);
    //         h_CW_NP2_ap1->SetStats(0);
    //         h_SM_ap1->SetStats(0);
    //         h_SM_ap1->SetAxisRange(-0.01,0.01,"Y");
    //         h_SM_ap1->Draw("E2");
    //         (h_CW_NP_ap1)->Draw("E2 SAME"); 
    //         (h_CW_NP2_ap1)->Draw("E2 SAME");
    //         c3.SaveAs("plots/rivet/CW_NP_NP2_ap1.png");
    //         //aps
    //         TCanvas c4("c4","",2000,2000);
    //         h_CW_NP_aps->SetStats(0);
    //         h_CW_NP2_aps->SetStats(0);
    //         h_SM_aps->SetStats(0);
    //         h_SM_aps->SetAxisRange(-0.01,0.01,"Y");
    //         h_SM_aps->Draw("E2");
    //         (h_CW_NP_aps)->Draw("E2 SAME"); 
    //         (h_CW_NP2_aps)->Draw("E2 SAME");
    //         c4.SaveAs("plots/rivet/CW_NP_NP2_aps.png");
    //         //ats
    //         TCanvas c5("c5","",2000,2000);
    //         h_CW_NP_ats->SetStats(0);
    //         h_CW_NP2_ats->SetStats(0);
    //         h_SM_ats->SetStats(0);
    //         h_SM_ats->SetAxisRange(0, M_PI);
    //         h_SM_ats->SetAxisRange(-0.1,0.1,"Y");
    //         h_SM_ats->Draw("E2");
    //         (h_CW_NP_ats)->Draw("E2 SAME"); 
    //         (h_CW_NP2_ats)->Draw("E2 SAME");
    //         c5.SaveAs("plots/rivet/CW_NP_NP2_ats.png");
    //         //at1
    //         TCanvas c6("c6","",2000,2000);
    //         h_CW_NP_at1->SetStats(0);
    //         h_CW_NP2_at1->SetStats(0);
    //         h_SM_at1->SetStats(0);
    //         h_SM_at1->SetAxisRange(0, M_PI);
    //         h_SM_at1->SetAxisRange(-0.1,0.1,"Y");
    //         h_SM_at1->Draw("E2");
    //         (h_CW_NP_at1)->Draw("E2 SAME"); 
    //         (h_CW_NP2_at1)->Draw("E2 SAME");
    //         c6.SaveAs("plots/rivet/CW_NP_NP2_at1.png");
    //         //at2
    //         TCanvas c7("c7","",2000,2000);
    //         h_CW_NP_at2->SetStats(0);
    //         h_CW_NP2_at2->SetStats(0);
    //         h_SM_at2->SetStats(0);
    //         h_SM_at2->SetAxisRange(0, M_PI);
    //         h_SM_at2->SetAxisRange(-0.1,0.1,"Y");
    //         h_SM_at2->Draw("E2");
    //         (h_CW_NP_at2)->Draw("E2 SAME"); 
    //         (h_CW_NP2_at2)->Draw("E2 SAME");
    //         c7.SaveAs("plots/rivet/CW_NP_NP2_at2.png");


    //     //CWTIL
    //         TCanvas c0("c0","",2000,2000);
    //         h_CWTIL_NP_dpj->SetAxisRange(-0.004,0.008, "Y");
    //         h_CWTIL_NP_dpj->SetStats(0);
    //         h_CWTIL_NP2_dpj->SetStats(0);
    //         (h_CWTIL_NP_dpj)->Draw("E2"); 
    //         (h_CWTIL_NP2_dpj)->Draw("E2 SAME");
    //         c0.SaveAs("plots/rivet/CWTIL_NP_NP2_dpj.png");

    //         TCanvas c1("c11","",2000,2000);
    //         h_tmp = (TH1*)(h_SM_dpj->Clone());
    //         auto h_diff_cwtil = (TH1*)(h_CWTIL_NP_dpj->Clone());
    //         h_diff_cwtil->Add(h_CWTIL_NP2_dpj);
    //         h_diff_cwtil->Add(h_SM_dpj);
    //         h_diff_cwtil = HistOverHist(h_diff_cw, h_SM_dpj);
    //         h_diff_cwtil->SetAxisRange(0.9,1.2,"Y");
    //         h_diff_cwtil->Draw("E2");
    //         c1.SaveAs("plots/rivet/CWTIL_SM_diff_dpj.png");


}

