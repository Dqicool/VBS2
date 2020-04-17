#include"libs/genAna.h"
#include<THStack.h>
#include<TLegend.h>
#include<TPaveText.h>
#define PRWLUMI (36.21 + 58.45 + 44.31)

void testat1(){
    TFile* insig = TFile::Open("output/stack_out/sigs/at1.root");
    TFile* inbkg = TFile::Open("output/stack_out/bkgs/at1.root");
    TFile* indat = TFile::Open("output/stack_out/Data.1516.root");
    TFile* intmp = TFile::Open("output/rebin_out/at1.root");

    THStack* s_sig = (THStack*)insig->Get("Angle_theta1_cut_stack");
    THStack* s_bkg = (THStack*)inbkg->Get("Angle_theta1_cut_stack");
    TH1D* h_dat = (TH1D*)indat->Get("h_Angle_theta1_pass_cut");
    TH1D* h_tmp = (TH1D*)intmp->Get("Angle_theta1_cut_hsig");
    int nbin = h_tmp->GetNbinsX();
    std::vector<double>sepa{};
    for(int i = 1; i<=nbin; i++){
        double aaa = h_tmp->GetBinLowEdge(i);
        sepa.push_back(aaa);
    }
    sepa.push_back(h_tmp->GetBinLowEdge(nbin) + h_tmp->GetBinWidth(nbin));
    std::string name = s_sig->begin().Next()->GetName();

    THStack s_rebin;
    //bkg
        TIter iter1(s_bkg->GetHists());
        while(1){
            auto h = (TH1*)(iter1.Next());
            if(h==nullptr) break;
            h = h->Rebin(nbin, name.data(), sepa.data());
            h->Scale(PRWLUMI);
            s_rebin.Add(h);
        }
    //sig
        TIter iter2(s_sig->GetHists());
        while(1){
            auto h = (TH1*)(iter2.Next());
            if(h==nullptr) break;
            h = h->Rebin(nbin, name.data(), sepa.data());
            h->Scale(PRWLUMI);
            s_rebin.Add(h);
        }
    
    h_dat = (TH1D*)h_dat->Rebin(nbin, name.data(), sepa.data());
    h_dat->Scale(PRWLUMI);
    TCanvas c1("c1","",1200,1200);
        c1.SetMargin(0.15,0.05,0.1,0.05);
        h_dat->SetStats(0);
        h_dat->SetAxisRange(0,30,"Y");
        h_dat->SetMarkerStyle(kFullCircle);
        h_dat->SetMarkerSize(2.);
        h_dat->SetLineWidth(4.);
        h_dat->SetLineColor(kBlack);
        h_dat->SetMarkerColor(kBlack);
        h_dat->GetXaxis()->SetTitle("#theta_{1} [ rad ]");
        h_dat->GetXaxis()->SetTitleFont(22);
        //h_dat->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta_{1}} [ fb / rad ]");
        h_dat->GetYaxis()->SetTitle("Events / bin");
        h_dat->GetYaxis()->SetTitleFont(22);
        h_dat->GetYaxis()->SetTitleOffset(1.5);
        h_dat->Draw("E1");
        s_rebin.Draw("HIST SAME");
        h_dat->Draw("E1 SAME");

        //legends
        TLegend legend(0.64,0.74,0.94,0.94);
        legend.AddEntry(h_dat, "Data", "lep");
        
        legend.AddEntry((s_sig->GetHists())->At(1), "EW ZZjj", "f");
        legend.AddEntry((s_sig->GetHists())->At(3), "QCD gg ZZjj", "f");
        legend.AddEntry((s_sig->GetHists())->At(7), "QCD qq ZZjj", "f");
        legend.AddEntry((s_bkg->GetHists())->At(1), "Other", "f");
        legend.SetTextFont(22);
        legend.SetLineColorAlpha(kBlack,0);
        legend.Draw();
        //text
        TLatex tex;
        tex.SetTextAlign(11);
        tex.SetTextFont(22);
        tex.SetTextSize(0.04);
        tex.DrawLatex(0.9,28,"#sqrt{s} = 13 TeV");
        tex.DrawLatex(0.9,26, "139 fb^{-1}");
    c1.SaveAs("debug.png");
}

void testdpj(){
    TFile* insig = TFile::Open("output/stack_out/sigs/jjDelPhi.root");
    TFile* inbkg = TFile::Open("output/stack_out/bkgs/jjDelPhi.root");
    TFile* indat = TFile::Open("output/stack_out/Data.1516.root");
    TFile* intmp = TFile::Open("output/rebin_out/dPhiJJ.root");

    THStack* s_sig = (THStack*)insig->Get("jj_dphi_cut_stack");
    TH1* h_sig = (TH1*)insig->Get("jj_dphi_cut_h");
    THStack* s_bkg = (THStack*)inbkg->Get("jj_dphi_cut_stack");
    TH1* h_bkg = (TH1*)inbkg->Get("jj_dphi_cut_h");
    TH1D* h_dat = (TH1D*)indat->Get("h_jj_dphi_pass_cut");
    TH1D* h_tmp = (TH1D*)intmp->Get("jj_dphi_cut_hsig");
    int nbin = h_tmp->GetNbinsX();
    std::vector<double>sepa{};
    for(int i = 1; i<=nbin; i++){
        double aaa = h_tmp->GetBinLowEdge(i);
        sepa.push_back(aaa);
    }
    sepa.push_back(h_tmp->GetBinLowEdge(nbin) + h_tmp->GetBinWidth(nbin));
    std::string name = s_sig->begin().Next()->GetName();
    THStack s_rebin;
    //bkg
        TIter iter1(s_bkg->GetHists());
        while(1){
            auto h = (TH1*)(iter1.Next());
            if(h==nullptr) break;
            h = h->Rebin(nbin, name.data(), sepa.data());
            h->Scale(PRWLUMI);
            s_rebin.Add(h);
        }
    //sig
        TIter iter2(s_sig->GetHists());
        while(1){
            auto h = (TH1*)(iter2.Next());
            if(h==nullptr) break;
            h = h->Rebin(nbin, name.data(), sepa.data());
            h->Scale(PRWLUMI);
            s_rebin.Add(h);
        }
        h_dat = (TH1D*)h_dat->Rebin(nbin, name.data(), sepa.data());
        h_dat->Scale(PRWLUMI);


        double err;
        double inte;
        h_sig = (TH1D*)h_sig->Rebin(nbin, name.data(), sepa.data());
        h_sig->Scale(PRWLUMI);
        inte = h_sig->IntegralAndError(1,h_sig->GetNbinsX(),err);
        std::cout<<"whole:\t"<<inte<<"+-"<<err<<std::endl;


    TCanvas c("c","",1200,1200);
        c.SetMargin(0.15,0.05,0.1,0.05);
        h_dat->SetStats(0);
        h_dat->SetAxisRange(0,30,"Y");
        h_dat->SetMarkerStyle(kFullCircle);
        h_dat->SetMarkerSize(2.);
        h_dat->SetLineWidth(4.);
        h_dat->SetLineColor(kBlack);
        h_dat->SetMarkerColor(kBlack);
        h_dat->GetXaxis()->SetTitle("#Delta#phi_{jj} [ rad ]");
        h_dat->GetXaxis()->SetTitleFont(22);
        //h_dat->GetYaxis()->SetTitle("#frac{d#sigma}{d#Delta#phi_{jj}} [ fb / rad ]");
        h_dat->GetYaxis()->SetTitle("Events / bin");
        h_dat->GetYaxis()->SetTitleFont(22);
        h_dat->Draw("E1");
        s_rebin.Draw("HIST SAME");
        h_dat->Draw("E1 SAME");

        //legends
        TLegend legend(0.64,0.74,0.94,0.94);
        legend.AddEntry(h_dat, "Data", "lep");
        legend.AddEntry((s_sig->GetHists())->At(1), "EW ZZjj", "f");
        legend.AddEntry((s_sig->GetHists())->At(3), "QCD gg ZZjj", "f");
        legend.AddEntry((s_sig->GetHists())->At(7), "QCD qq ZZjj", "f");
        legend.AddEntry((s_bkg->GetHists())->At(1), "Other", "f");
        legend.SetTextFont(22);
        legend.SetLineColorAlpha(kBlack,0);
        legend.Draw();
        //text
        TLatex texx;
        texx.SetTextAlign(11);
        texx.SetTextFont(22);
        texx.SetTextSize(0.04);
        texx.DrawLatex(-1.342,28,"#sqrt{s} = 13 TeV");
        texx.DrawLatex(-1.342,26, "139 fb^{-1}");
    c.SaveAs("dpjj.png");
}


int main(){
    testat1();
    testdpj();
}