#include"libs/genAna.h"
#include"libs/ErrorProp.h"
#include<cmath>

enum MCType { SM, NPLinear, NPQuad };
TH1* getPlot(const char* file, const char* plot, double eft_C, Color_t col, MCType types){
    TH1::SetDefaultSumw2();
    
    TFile* infile = TFile::Open(file, "read");
    TH1D * h = (TH1D*)infile->Get(plot);
    switch(types)
    {
        case SM      :                                          break;
        case NPLinear: h->Scale(eft_C / 5.0);                   break;
        case NPQuad  : h->Scale((eft_C / 5.0)*(eft_C / 5.0));   break;
    }
    h->SetFillColor(col);
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    return h;
}

TH1* HistOverHist(TH1* up, TH1* down){
    for(int i=1; i <= up->GetNbinsX(); i++){
        if (down->GetBinContent(i) != 0){
            double cont = up->GetBinContent(i) / down->GetBinContent(i);
            double err = ErrADiviB(cont, up->GetBinContent(i), down->GetBinContent(i), up->GetBinError(i), down->GetBinError(i), 0);
            up->SetBinContent(i, cont);
            up->SetBinError(i, err);
        }
        else {up->SetBinContent(i, 0);}
    }
    auto ret = (TH1*)(up->Clone());
    return ret;
}

void drawRatio(const char* fSM, const char* fNP, const char* fNP2, const char* hSM, const char* hNP, const char* hNP2, double eftC, const char* fout, double yext, double xmin, double xmax){
    TH1* h_SM  = getPlot(fSM,   hSM,    eftC, kRed,    SM);
    TH1* h_NP  = getPlot(fNP,   hNP,    eftC, kBlue,   NPQuad);
    TH1* h_NP2 = getPlot(fNP2,  hNP2,   eftC, kGreen,  NPLinear);
    //calculate just Ratio
    auto h_ratio_NP = HistOverHist(h_NP, h_SM);
    auto h_ratio_NP2 = HistOverHist(h_NP2, h_SM);

    TCanvas c0("c0","",2000,2000);
    h_ratio_NP2->SetLineStyle(10);
    h_ratio_NP2->SetLineWidth(10);
    h_ratio_NP2->SetFillColorAlpha(kGreen, 0);
    h_ratio_NP2->SetAxisRange(-yext,yext, "Y");
    h_ratio_NP2->SetAxisRange(xmin, xmax, "X");
    h_ratio_NP2->SetStats(0);

    h_ratio_NP->SetLineWidth(10);
    h_ratio_NP->SetLineStyle(2);
    h_ratio_NP->SetFillColorAlpha(kBlue, 0);
    h_ratio_NP->SetAxisRange(-yext,yext, "Y");
    //h_ratio_NP->SetAxisRange(xmin, xmax, "X");
    //std::cout<< h_ratio_NP->GetBinLowEdge(1) <<std::endl;
    h_ratio_NP->SetStats(0);
    //h_SM_rebin->Draw();
    h_ratio_NP->Draw("HIST");
    h_ratio_NP2->Draw("SAME HIST");
    auto h_tmp1 = (TH1*)h_ratio_NP->Clone();
    auto h_tmp2 = (TH1*)h_ratio_NP2->Clone();
    h_tmp1->SetFillColorAlpha(kBlue, 0.3);
    h_tmp2->SetFillColorAlpha(kGreen, 0.3);
    h_tmp1->Draw("SAME E2");
    h_tmp2->Draw("SAME E2");
    c0.SaveAs(fout);
}

int main(){
    //dphijj 
        //CW 
        drawRatio("output/rebin_out/dPhiJJ.root", "output/rebin_out/dPhiJJ.root",    "output/rebin_out/dPhiJJ.root",    "jj_truthBorn_dphi_cut_h", "jj_dphi_cut_h_CW_NP", "jj_dphi_cut_h_CW_NP2",           1, "plots/rivet/Dphijj_CW_ratio.png",     0.5, -M_PI, M_PI);
        //CWtil
        drawRatio("output/rebin_out/dPhiJJ.root", "output/rebin_out/dPhiJJ.root",    "output/rebin_out/dPhiJJ.root",    "jj_truthBorn_dphi_cut_h", "jj_dphi_cut_h_CWtil_NP", "jj_dphi_cut_h_CWtil_NP2",     1, "plots/rivet/Dphijj_CWtil_ratio.png",  0.5, -M_PI, M_PI);
        //CHW 
        drawRatio("output/rebin_out/dPhiJJ.root", "output/rebin_out/dPhiJJ.root",    "output/rebin_out/dPhiJJ.root",    "jj_truthBorn_dphi_cut_h", "jj_dphi_cut_h_CHW_NP", "jj_dphi_cut_h_CHW_NP2",         1, "plots/rivet/Dphijj_CHW_ratio.png",    0.5, -M_PI, M_PI);
        //CHWTIL
        drawRatio("output/rebin_out/dPhiJJ.root", "output/rebin_out/dPhiJJ.root",    "output/rebin_out/dPhiJJ.root",    "jj_truthBorn_dphi_cut_h", "jj_dphi_cut_h_CHWtil_NP", "jj_dphi_cut_h_CHWtil_NP2",   1, "plots/rivet/Dphijj_CHWtil_ratio.png", 0.5, -M_PI, M_PI);

    //theta_star
        drawRatio("output/rebin_out/ats.root", "output/rebin_out/ats.root",    "output/rebin_out/ats.root",    "Angle_truthBorn_theta_star_cut_h", "Angle_theta_star_cut_h_CW_NP", "Angle_theta_star_cut_h_CW_NP2",     1, "plots/rivet/ats_CW_ratio.png",     0.2, -M_PI, M_PI);
        //CWtil
        drawRatio("output/rebin_out/ats.root", "output/rebin_out/ats.root",    "output/rebin_out/ats.root",    "Angle_truthBorn_theta_star_cut_h", "Angle_theta_star_cut_h_CWtil_NP", "Angle_theta_star_cut_h_CWtil_NP2",     1, "plots/rivet/ats_CWtil_ratio.png",  0.2, -M_PI, M_PI);
        //CHW 
        drawRatio("output/rebin_out/ats.root", "output/rebin_out/ats.root",    "output/rebin_out/ats.root",    "Angle_truthBorn_theta_star_cut_h", "Angle_theta_star_cut_h_CHW_NP", "Angle_theta_star_cut_h_CHW_NP2",     1, "plots/rivet/ats_CHW_ratio.png",    0.2, -M_PI, M_PI);
        //CHWTIL
        drawRatio("output/rebin_out/ats.root", "output/rebin_out/ats.root",    "output/rebin_out/ats.root",    "Angle_truthBorn_theta_star_cut_h", "Angle_theta_star_cut_h_CHWtil_NP", "Angle_theta_star_cut_h_CHWtil_NP2",     1, "plots/rivet/ats_CHWtil_ratio.png", 0.2, -M_PI, M_PI);

    //theta1
        drawRatio("output/rebin_out/at1.root", "output/rebin_out/at1.root",    "output/rebin_out/at1.root",    "Angle_truthBorn_theta1_cut_h", "Angle_theta1_cut_h_CW_NP", "Angle_theta1_cut_h_CW_NP2",     1, "plots/rivet/at1_CW_ratio.png",     0.2, -M_PI, M_PI);
        //CWtil
        drawRatio("output/rebin_out/at1.root", "output/rebin_out/at1.root",    "output/rebin_out/at1.root",    "Angle_truthBorn_theta1_cut_h", "Angle_theta1_cut_h_CWtil_NP", "Angle_theta1_cut_h_CWtil_NP2",     1, "plots/rivet/at1_CWtil_ratio.png",  0.2, -M_PI, M_PI);
        //CHW 
        drawRatio("output/rebin_out/at1.root", "output/rebin_out/at1.root",    "output/rebin_out/at1.root",    "Angle_truthBorn_theta1_cut_h", "Angle_theta1_cut_h_CHW_NP", "Angle_theta1_cut_h_CHW_NP2",     1, "plots/rivet/at1_CHW_ratio.png",    0.2, -M_PI, M_PI);
        //CHWTIL
        drawRatio("output/rebin_out/at1.root", "output/rebin_out/at1.root",    "output/rebin_out/at1.root",    "Angle_truthBorn_theta1_cut_h", "Angle_theta1_cut_h_CHWtil_NP", "Angle_theta1_cut_h_CHWtil_NP2",     1, "plots/rivet/at1_CHWtil_ratio.png", 0.2, -M_PI, M_PI);

    //theta2
        drawRatio("output/rebin_out/at2.root", "output/rebin_out/at2.root",    "output/rebin_out/at2.root",    "Angle_truthBorn_theta2_cut_h", "Angle_theta2_cut_h_CW_NP", "Angle_theta2_cut_h_CW_NP2",     1, "plots/rivet/at2_CW_ratio.png",     0.2, -M_PI, M_PI);
        //CWtil
        drawRatio("output/rebin_out/at2.root", "output/rebin_out/at2.root",    "output/rebin_out/at2.root",    "Angle_truthBorn_theta2_cut_h", "Angle_theta2_cut_h_CWtil_NP", "Angle_theta2_cut_h_CWtil_NP2",     1, "plots/rivet/at2_CWtil_ratio.png",  0.2, -M_PI, M_PI);
        //CHW 
        drawRatio("output/rebin_out/at2.root", "output/rebin_out/at2.root",    "output/rebin_out/at2.root",    "Angle_truthBorn_theta2_cut_h", "Angle_theta2_cut_h_CHW_NP", "Angle_theta2_cut_h_CHW_NP2",     1, "plots/rivet/at2_CHW_ratio.png",    0.2, -M_PI, M_PI);
        //CHWTIL
        drawRatio("output/rebin_out/at2.root", "output/rebin_out/at2.root",    "output/rebin_out/at2.root",    "Angle_truthBorn_theta2_cut_h", "Angle_theta2_cut_h_CHWtil_NP", "Angle_theta2_cut_h_CHWtil_NP2",     1, "plots/rivet/at2_CHWtil_ratio.png", 0.2, -M_PI, M_PI);
    //phi_star
        //CW
        drawRatio("output/rebin_out/aps.root", "output/rebin_out/aps.root",    "output/rebin_out/aps.root",    "Angle_truthBorn_phi_star_cut_h", "Angle_phi_star_cut_h_CW_NP", "Angle_phi_star_cut_h_CW_NP2",     1, "plots/rivet/aps_CW_ratio.png",     0.2, -M_PI, M_PI);
        //CWtil
        drawRatio("output/rebin_out/aps.root", "output/rebin_out/aps.root",    "output/rebin_out/aps.root",    "Angle_truthBorn_phi_star_cut_h", "Angle_phi_star_cut_h_CWtil_NP", "Angle_phi_star_cut_h_CWtil_NP2",     1, "plots/rivet/aps_CWtil_ratio.png",  0.2, -M_PI, M_PI);
        //CHW 
        drawRatio("output/rebin_out/aps.root", "output/rebin_out/aps.root",    "output/rebin_out/aps.root",    "Angle_truthBorn_phi_star_cut_h", "Angle_phi_star_cut_h_CHW_NP", "Angle_phi_star_cut_h_CHW_NP2",     1, "plots/rivet/aps_CHW_ratio.png",    0.2, -M_PI, M_PI);
        //CHWTIL
        drawRatio("output/rebin_out/aps.root", "output/rebin_out/aps.root",    "output/rebin_out/aps.root",    "Angle_truthBorn_phi_star_cut_h", "Angle_phi_star_cut_h_CHWtil_NP", "Angle_phi_star_cut_h_CHWtil_NP2",     1, "plots/rivet/aps_CHWtil_ratio.png", 0.2, -M_PI, M_PI);

    //theta1
        drawRatio("output/rebin_out/ap0.root", "output/rebin_out/ap0.root",    "output/rebin_out/ap0.root",    "Angle_truthBorn_phi0_cut_h", "Angle_phi0_cut_h_CW_NP", "Angle_phi0_cut_h_CW_NP2",     1, "plots/rivet/ap0_CW_ratio.png",     0.2, -M_PI, M_PI);
        //CWtil
        drawRatio("output/rebin_out/ap0.root", "output/rebin_out/ap0.root",    "output/rebin_out/ap0.root",    "Angle_truthBorn_phi0_cut_h", "Angle_phi0_cut_h_CWtil_NP", "Angle_phi0_cut_h_CWtil_NP2",     1, "plots/rivet/ap0_CWtil_ratio.png",  0.2, -M_PI, M_PI);
        //CHW 
        drawRatio("output/rebin_out/ap0.root", "output/rebin_out/ap0.root",    "output/rebin_out/ap0.root",    "Angle_truthBorn_phi0_cut_h", "Angle_phi0_cut_h_CHW_NP", "Angle_phi0_cut_h_CHW_NP2",     1, "plots/rivet/ap0_CHW_ratio.png",    0.2, -M_PI, M_PI);
        //CHWTIL
        drawRatio("output/rebin_out/ap0.root", "output/rebin_out/ap0.root",    "output/rebin_out/ap0.root",    "Angle_truthBorn_phi0_cut_h", "Angle_phi0_cut_h_CHWtil_NP", "Angle_phi0_cut_h_CHWtil_NP2",     1, "plots/rivet/ap0_CHWtil_ratio.png", 0.2, -M_PI, M_PI);

    //phi2
        drawRatio("output/rebin_out/ap1.root", "output/rebin_out/ap1.root",    "output/rebin_out/ap1.root",    "Angle_truthBorn_phi1_cut_h", "Angle_phi1_cut_h_CW_NP", "Angle_phi1_cut_h_CW_NP2",     1, "plots/rivet/ap1_CW_ratio.png",     0.2, -M_PI, M_PI);
        //CWtil
        drawRatio("output/rebin_out/ap1.root", "output/rebin_out/ap1.root",    "output/rebin_out/ap1.root",    "Angle_truthBorn_phi1_cut_h", "Angle_phi1_cut_h_CWtil_NP", "Angle_phi1_cut_h_CWtil_NP2",     1, "plots/rivet/ap1_CWtil_ratio.png",  0.2, -M_PI, M_PI);
        //CHW 
        drawRatio("output/rebin_out/ap1.root", "output/rebin_out/ap1.root",    "output/rebin_out/ap1.root",    "Angle_truthBorn_phi1_cut_h", "Angle_phi1_cut_h_CHW_NP", "Angle_phi1_cut_h_CHW_NP2",     1, "plots/rivet/ap1_CHW_ratio.png",    0.2, -M_PI, M_PI);
        //CHWTIL
        drawRatio("output/rebin_out/ap1.root", "output/rebin_out/ap1.root",    "output/rebin_out/ap1.root",    "Angle_truthBorn_phi1_cut_h", "Angle_phi1_cut_h_CHWtil_NP", "Angle_phi1_cut_h_CHWtil_NP2",     1, "plots/rivet/ap1_CHWtil_ratio.png", 0.2, -M_PI, M_PI);
    //mjj
        drawRatio("output/rebin_out/ap1.root", "output/rebin_out/ap1.root",    "output/rebin_out/ap1.root",    "Angle_truthBorn_phi1_cut_h", "Angle_phi1_cut_h_CW_NP", "Angle_phi1_cut_h_CW_NP2",     1, "plots/rivet/ap1_CW_ratio.png",     0.2, -M_PI, M_PI);
        //CWtil
        drawRatio("output/rebin_out/ap1.root", "output/rebin_out/ap1.root",    "output/rebin_out/ap1.root",    "Angle_truthBorn_phi1_cut_h", "Angle_phi1_cut_h_CWtil_NP", "Angle_phi1_cut_h_CWtil_NP2",     1, "plots/rivet/ap1_CWtil_ratio.png",  0.2, -M_PI, M_PI);
        //CHW 
        drawRatio("output/rebin_out/ap1.root", "output/rebin_out/ap1.root",    "output/rebin_out/ap1.root",    "Angle_truthBorn_phi1_cut_h", "Angle_phi1_cut_h_CHW_NP", "Angle_phi1_cut_h_CHW_NP2",     1, "plots/rivet/ap1_CHW_ratio.png",    0.2, -M_PI, M_PI);
        //CHWTIL
        drawRatio("output/rebin_out/ap1.root", "output/rebin_out/ap1.root",    "output/rebin_out/ap1.root",    "Angle_truthBorn_phi1_cut_h", "Angle_phi1_cut_h_CHWtil_NP", "Angle_phi1_cut_h_CHWtil_NP2",     1, "plots/rivet/ap1_CHWtil_ratio.png", 0.2, -M_PI, M_PI);

    //m4l
}