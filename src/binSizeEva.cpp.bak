#include "genAna.h"
#include <THStack.h>

#define NEVTSPERBIN 10

using namespace std;

std::vector<double> mulScalar(std::vector<double> mat, double x){
    std::vector<double> ret{};
    for(int i = 0; i < mat.size(); i++){
        ret.push_back(mat[i]*x);
    }
    return ret;
}

std::vector<double> pluScalar(std::vector<double> mat, double x){
    std::vector<double> ret{};
    for(int i = 0; i < mat.size(); i++){
        ret.push_back(mat[i] + x);
    }
    return ret;
}

void m4lBinSizeEva()
{
    TH1::SetDefaultSumw2();
    TFile* infile1 = TFile::Open("output/draw_out/Data1516.root", "read");
    TH1D * h_llll_m_data_cut = (TH1D*)infile1->Get("h_llll_m_pass_cut");

    TFile* infile2 = TFile::Open("output/stack_out/sig_m4l.root", "read");
    TFile* infile3 = TFile::Open("output/stack_out/bak_m4l.root", "read");
    TH1D * h_llll_m_bak_cut = (TH1D*)infile3->Get("llll_m_cut_h");
    TH1D * h_llll_m_cut = (TH1D*)infile2->Get("llll_m_cut_h");
    //cout<<(h_llll_m_cut->Integral())<<endl;
    auto nbin = h_llll_m_cut->GetNbinsX();
    auto wbin = h_llll_m_cut->GetBinWidth(1);
    std::vector<double> sepa{0, 110, 140};
    for(int i=140; i<nbin; ){
        double store = 0;
        while (store < NEVTSPERBIN && i < nbin){
            store += h_llll_m_cut->GetBinContent(i);
            i++;
        }
        sepa.push_back(i);
    }
    if(h_llll_m_cut->Integral((int)sepa[sepa.size()-2], nbin) < 5) 
        {sepa.erase(sepa.end() - 2);}
    sepa = mulScalar(sepa, wbin);
    auto h_rebin_m4l_cut = h_llll_m_cut->Rebin(sepa.size()-1, "h_rebin_m4l_cut",&sepa[0]);
    auto h_rebin_llll_m_data_cut = h_llll_m_data_cut->Rebin(sepa.size()-1, "h_rebin_m4l_data_cut",&sepa[0]);
    auto h_rebin_llll_m_bak_cut = h_llll_m_bak_cut->Rebin(sepa.size()-1, "h_rebin_m4l_bak_cut",&sepa[0]);
    cout<<h_rebin_llll_m_data_cut->Integral()<<endl;
    h_rebin_llll_m_data_cut->Add(h_rebin_llll_m_bak_cut, -1);
    cout<<h_rebin_llll_m_data_cut->Integral()<<endl;
    h_rebin_m4l_cut->SetAxisRange(0,20,"Y");
    h_rebin_m4l_cut->SetStats(0);
    h_rebin_m4l_cut->Draw("");
    h_llll_m_cut->Scale(10.0);
    h_llll_m_cut->SetFillColor(kRed);
    h_llll_m_cut->SetMarkerColor(kRed);
    h_llll_m_cut->SetLineColor(kRed);
    h_llll_m_cut->SetStats(0);
    h_llll_m_cut->Draw("Histo same");
    
    //save in files
    TFile* outfile = TFile::Open("output/histo_out/dem_unfold_m4l.root", "recreate");
    h_rebin_m4l_cut->Write();
    h_rebin_llll_m_data_cut->Write();
    
    outfile->Close();
}

void mjjBinSizeEva()
{
    TH1::SetDefaultSumw2();
    TFile* infile1 = TFile::Open("output/stack_out/sig_mjj.root", "read");
    TH1D * h_jj_m_cut = (TH1D*)infile1->Get("jj_m_cut_h");

    TFile* infile2 = TFile::Open("output/draw_out/Data1516.root", "read");
    TH1D * h_mjj_data_cut = (TH1D*)infile2->Get("h_jj_m_pass_cut");

    TFile* infile3 = TFile::Open("output/stack_out/bak_mjj.root", "read");
    TH1D * h_mjj_bak_cut = (TH1D*)infile3->Get("jj_m_cut_h");
    //cout<<(h_jj_m_cut->Integral())<<endl;
    auto nbin = h_jj_m_cut->GetNbinsX();
    auto wbin = h_jj_m_cut->GetBinWidth(1);
    auto xmin = h_jj_m_cut->GetBinLowEdge(1);
    std::vector<double> sepa{0};
    for(int i=0; i<nbin; ){
        double store = 0;
        while (store < NEVTSPERBIN && i < nbin){
            store += h_jj_m_cut->GetBinContent(i);
            i++;
        }
        sepa.push_back(i);
    }
    if(h_jj_m_cut->Integral((int)sepa[sepa.size()-2], nbin) < 5) 
        {sepa.erase(sepa.end() - 2);}
    sepa = mulScalar(sepa, wbin);
    sepa = pluScalar(sepa,xmin);
    auto h_rebin_mjj_cut = h_jj_m_cut->Rebin(sepa.size()-1, "h_rebin_mjj_cut",&sepa[0]);
    auto h_rebin_mjj_data_cut =  h_mjj_data_cut->Rebin(sepa.size()-1, "h_rebin_mjj_data_cut",&sepa[0]);
    auto h_rebin_mjj_bak_cut =  h_mjj_bak_cut->Rebin(sepa.size()-1, "h_rebin_mjj_bak_cut",&sepa[0]);
    h_rebin_mjj_data_cut->Add(h_rebin_mjj_bak_cut, -1);
    //cout<<h_rebin_mjj_cut->Integral()<<endl;
    h_rebin_mjj_cut->SetStats(0);
    h_rebin_mjj_cut->SetAxisRange(0,20,"Y");
    h_rebin_mjj_cut->Draw("");
    h_jj_m_cut->Scale(10.0);
    h_jj_m_cut->SetFillColor(kRed);
    h_jj_m_cut->SetMarkerColor(kRed);
    h_jj_m_cut->SetLineColor(kRed);
    h_jj_m_cut->SetStats(0);
    h_jj_m_cut->Draw("Histo same");
    
    //save in files
    TFile* outfile = TFile::Open("output/histo_out/dem_unfold_mjj.root", "recreate");
    h_rebin_mjj_cut->Write();
    h_rebin_mjj_data_cut->Write();
    outfile->Close();
}

void delphijjBinSizeEva(){
    TH1::SetDefaultSumw2();
    TFile* infile = TFile::Open("output/stack_out/sig_jjDelPhi.root", "read");
    TH1D * h_jj_delta_phi_cut = (TH1D*)infile->Get("jj_delta_phi_cut_h");

    TFile* infile2 = TFile::Open("output/draw_out/Data1516.root", "read");
    TH1D * h_jj_delta_phi_data_cut = (TH1D*)infile2->Get("h_jj_delta_phi_pass_cut");

    TFile* infile3 = TFile::Open("output/stack_out/bak_jjDelPhi.root", "read");
    TH1D * h_jj_delta_phi_bak_cut = (TH1D*)infile3->Get("jj_delta_phi_cut_h");
    //cout<<(h_jj_delta_phi_cut->Integral())<<endl;
    auto nbin = h_jj_delta_phi_cut->GetNbinsX();
    auto wbin = h_jj_delta_phi_cut->GetBinWidth(1);
    auto xmin = h_jj_delta_phi_cut->GetBinLowEdge(1);

    std::vector<double> sepaR{(double)nbin/2};
    std::vector<double> sepaL{};
    for(int i=nbin/2 +1; i<nbin; ){
        double store = 0;
        while (store < NEVTSPERBIN && i < nbin){
            store += h_jj_delta_phi_cut->GetBinContent(i);
            i++;
        }
        sepaR.push_back(i);
    }
    if(h_jj_delta_phi_cut->Integral((int)sepaR[sepaR.size()-2], nbin) < 5) 
        {sepaR.erase(sepaR.end() - 2);}
    sepaL = sepaR;
    std::reverse(sepaL.begin(), sepaL.end());
    sepaL = pluScalar(mulScalar(sepaL, -1), nbin);

    std::vector<double> sepa{};
    sepa.insert(sepa.end(), sepaL.begin(), sepaL.end()-1);
    sepa.insert(sepa.end(), sepaR.begin(), sepaR.end());
    sepa = mulScalar(sepa, wbin);
    sepa = pluScalar(sepa,xmin);
    auto h_rebin_delphijj_cut = h_jj_delta_phi_cut->Rebin(sepa.size()-1, "h_rebin_delphijj_cut",&sepa[0]);
    auto h_rebin_jj_delta_phi_data_cut = h_jj_delta_phi_data_cut->Rebin(sepa.size()-1, "h_rebin_delphijj_data_cut",&sepa[0]);
    auto h_rebin_jj_delta_phi_bak_cut = h_jj_delta_phi_bak_cut->Rebin(sepa.size()-1, "h_rebin_delphijj_bak_cut",&sepa[0]);
    h_rebin_jj_delta_phi_data_cut->Add(h_rebin_jj_delta_phi_bak_cut, -1);
    //cout<<h_rebin_delphijj_cut->Integral()<<endl;
    h_rebin_delphijj_cut->SetAxisRange(0,20,"Y");
    h_rebin_delphijj_cut->SetStats(0);
    h_rebin_delphijj_cut->Draw("");
    h_jj_delta_phi_cut->Scale(5.0);
    h_jj_delta_phi_cut->SetFillColor(kRed);
    h_jj_delta_phi_cut->SetMarkerColor(kRed);
    h_jj_delta_phi_cut->SetLineColor(kRed);
    h_jj_delta_phi_cut->SetStats(0);
    h_jj_delta_phi_cut->Draw("Histo same");
    
    //save in files
    TFile* outfile = TFile::Open("output/histo_out/dem_unfold_delphijj.root", "recreate");
    h_rebin_jj_delta_phi_data_cut->Write();
    h_rebin_delphijj_cut->Write();
    outfile->Close();
}

int main(){
    m4lBinSizeEva();
    mjjBinSizeEva();
    delphijjBinSizeEva();
}