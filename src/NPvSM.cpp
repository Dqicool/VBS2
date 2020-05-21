#include"libs/genAna.h"
#include"libs/ErrorProp.h"
#include<cmath>
#include<TLegend.h>
#include<TLine.h>
#include"libs/ErrorProp.h"

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

void drawRatio(const char* fSM, const char* fNP, const char* fNP2, const char* hSM, const char* hNP, const char* hNP2, double eftC, const char* fout, double yext, double xmin, double xmax, const char* xlabel = "", const char* ylabel = "", const char* wilson_coe = ""){
    TH1* h_SM  = getPlot(fSM,   hSM,    eftC, kRed,    SM);
    TH1* h_NP  = getPlot(fNP,   hNP,    eftC, kBlue,   NPQuad);
    TH1* h_NP2 = getPlot(fNP2,  hNP2,   eftC, kGreen,  NPLinear);
    if((std::string)hSM == (std::string)"jj_truthBorn_dphi_cut_h"){
        auto nbin = h_NP->GetNbinsX();
        
        double errxsNP, errxsNP2, errxsSM;
        double xsEFTNP = h_NP->IntegralAndError(1, nbin, errxsNP);
        double xsEFTNP2 =h_NP2->IntegralAndError(1, nbin, errxsNP2);
        double xsEFT = xsEFTNP + xsEFTNP2;
        double xsSM  = h_SM->IntegralAndError(1, nbin, errxsSM);
        double errxsEFT = ErrAPlusB(errxsNP, errxsNP2,0);
        double errratio = ErrADiviB(xsEFT/xsSM, xsEFT, xsSM, errxsEFT, errxsSM, 0);
        std::cout<<"-------------------"<<std::endl;
        std::cout<<"for "<<hNP<<std::endl;
        std::cout<<"Lin xs\t"<<xsEFTNP2<<"+-"<<errxsNP2<<std::endl;
        std::cout<<"Qua xs\t"<<xsEFTNP<<"+-"<<errxsNP<<std::endl;
        std::cout<<"SM  xs\t"<<xsSM<<"+-"<<errxsSM<<std::endl;
        std::cout<<"EFT/SM\t"<<xsEFT/xsSM<<"+-"<<errratio<<std::endl;
        
        double errplusNP, errminusNP, errplusNP2, errminusNP2, errplusSM, errminusSM; 
        double xsminus=h_NP->IntegralAndError(1,nbin/2, errminusNP) + h_NP2->IntegralAndError(1,nbin/2, errminusNP2)+h_SM->IntegralAndError(1,nbin/2, errminusSM);
        double xsplus=h_NP->IntegralAndError(nbin/2+1, nbin, errplusNP) + h_NP2->IntegralAndError(nbin/2+1, nbin, errplusNP2)+h_SM->IntegralAndError(nbin/2+1, nbin, errplusSM);
        double xsplusSM = h_SM->IntegralAndError(nbin/2+1, nbin, errplusSM);
        double xsminusSM = h_SM->IntegralAndError(1,nbin/2, errminusSM);
        double asySM    = (xsplusSM - xsminusSM) / (xsplusSM + xsminusSM);
        double asy      = (xsplus - xsminus) / (xsplus + xsminus);
        double errplus = ErrAPlusB(ErrAPlusB(errplusNP, errplusNP2,0), errplusSM,0);
        double errminus = ErrAPlusB(ErrAPlusB(errminusNP, errminusNP2,0), errminusSM,0);

        double errup    = ErrAPlusB(errplus,errminus,0);
        double errdown  = ErrAMinusB(errplus,errminus,0);
        double errasy   = ErrADiviB(asy, (xsplus - xsminus), (xsplus + xsminus), errup, errdown, 0);

        double errupSM    = ErrAPlusB(errplusSM,errminusSM,0);
        double errdownSM  = ErrAMinusB(errplusSM,errminusSM,0);
        double errasySM   = ErrADiviB(asySM, (xsplusSM - xsminusSM), (xsplusSM + xsminusSM), errupSM, errdownSM, 0);
        std::cout<<"Asymm SM\t"<<asySM<<"+-"<<errasySM<<std::endl;
        std::cout<<"Asymm\t"<<asy<<"+-"<<errasy<<std::endl;
        std::cout<<"-------------------"<<std::endl;

    }
    //calculate just Ratio
    auto h_ratio_NP = HistOverHist(h_NP, h_SM);
    auto h_ratio_NP2 = HistOverHist(h_NP2, h_SM);

    TCanvas c0("c0","",2000,1500);
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
    h_ratio_NP->GetXaxis()->SetTitle(xlabel);
    h_ratio_NP->GetXaxis()->SetTitleFont(22);
    h_ratio_NP->GetYaxis()->SetTitle(ylabel);
    h_ratio_NP->GetYaxis()->SetTitleFont(22);
    h_ratio_NP->Draw("HIST");
    h_ratio_NP2->Draw("SAME HIST");
    auto h_ratio_tot = (TH1*)(h_ratio_NP->Clone());
    h_ratio_tot->Add(h_ratio_NP2);
    h_ratio_tot->SetLineStyle(7);
    h_ratio_tot->SetLineColor(kRed);
    h_ratio_tot->Draw("SAME HIST");
    auto h_tmp1 = (TH1*)h_ratio_NP->Clone();
    auto h_tmp2 = (TH1*)h_ratio_NP2->Clone();
    auto h_tmp3 = (TH1*)h_ratio_tot->Clone();
    h_tmp1->SetFillColorAlpha(kBlue, 0.3);
    h_tmp2->SetFillColorAlpha(kGreen, 0.3);
    h_tmp3->SetFillColorAlpha(kRed, 0.3);
    h_tmp1->Draw("SAME E2");
    h_tmp2->Draw("SAME E2");
    //h_tmp3->Draw("SAME E2");
    TLegend leg(0.7,0.7,0.89,0.89);
    leg.AddEntry(h_tmp1, "EFT Quad", "lpf");
    leg.AddEntry(h_tmp2, "EFT Linear", "lpf");
    leg.AddEntry(h_tmp3, "EFT total", "lp");
    leg.SetTextFont(22);
    leg.SetLineColorAlpha(kBlack,0);
    leg.Draw();
    TLine horiline(xmin,0,xmax,0);
    horiline.SetLineWidth(2);
    horiline.Draw();
    TLatex texx;
    texx.SetTextAlign(11);
    texx.SetTextFont(22);
    texx.SetTextSize(0.05);
    texx.DrawLatex((xmax - xmin)*0.2 + xmin, yext*0.8, wilson_coe);
    texx.DrawLatex((xmax - xmin)*0.2 + xmin, yext*0.65, "#Lambda = 1000");
    c0.SaveAs(fout);
    
}

int main(){
    //dphijj 
        //CW 
        drawRatio("output/rebin_out/dPhiJJ.root", "output/rebin_out/dPhiJJ.root",    "output/rebin_out/dPhiJJ.root",    "jj_truthBorn_dphi_cut_h", "jj_dphi_cut_h_CW_NP", "jj_dphi_cut_h_CW_NP2",           1, "plots/rivet/Dphijj_CW_ratio.png",     0.5, -M_PI, M_PI, "#Delta#phi_{jj} [ rad ]", "#frac{d#sigma^{EFT}}{d#sigma^{SM}}", "C_{W} = 1");
        //CWtil
        drawRatio("output/rebin_out/dPhiJJ.root", "output/rebin_out/dPhiJJ.root",    "output/rebin_out/dPhiJJ.root",    "jj_truthBorn_dphi_cut_h", "jj_dphi_cut_h_CWtil_NP", "jj_dphi_cut_h_CWtil_NP2",     1, "plots/rivet/Dphijj_CWtil_ratio.png",  0.5, -M_PI, M_PI, "#Delta#phi_{jj} [ rad ]", "#frac{d#sigma^{EFT}}{d#sigma^{SM}}", "C_{#tilde{W}} = 1");
        //CHW 
        drawRatio("output/rebin_out/dPhiJJ.root", "output/rebin_out/dPhiJJ.root",    "output/rebin_out/dPhiJJ.root",    "jj_truthBorn_dphi_cut_h", "jj_dphi_cut_h_CHW_NP", "jj_dphi_cut_h_CHW_NP2",         1, "plots/rivet/Dphijj_CHW_ratio.png",    0.5, -M_PI, M_PI, "#Delta#phi_{jj} [ rad ]", "#frac{d#sigma^{EFT}}{d#sigma^{SM}}", "C_{HW} = 1");
        //CHWTIL
        drawRatio("output/rebin_out/dPhiJJ.root", "output/rebin_out/dPhiJJ.root",    "output/rebin_out/dPhiJJ.root",    "jj_truthBorn_dphi_cut_h", "jj_dphi_cut_h_CHWtil_NP", "jj_dphi_cut_h_CHWtil_NP2",   1, "plots/rivet/Dphijj_CHWtil_ratio.png", 0.5, -M_PI, M_PI, "#Delta#phi_{jj} [ rad ]", "#frac{d#sigma^{EFT}}{d#sigma^{SM}}", "C_{H#tilde{W}} = 1");
        //CHWB
        drawRatio("output/rebin_out/dPhiJJ.root", "output/rebin_out/dPhiJJ.root",    "output/rebin_out/dPhiJJ.root",    "jj_truthBorn_dphi_cut_h", "jj_dphi_cut_h_CHWB_NP", "jj_dphi_cut_h_CHWB_NP2",   1, "plots/rivet/Dphijj_CHWB_ratio.png", 0.5, -M_PI, M_PI, "#Delta#phi_{jj} [ rad ]", "#frac{d#sigma^{EFT}}{d#sigma^{SM}}", "C_{HWB} = 1");
        //CHWBtil
        drawRatio("output/rebin_out/dPhiJJ.root", "output/rebin_out/dPhiJJ.root",    "output/rebin_out/dPhiJJ.root",    "jj_truthBorn_dphi_cut_h", "jj_dphi_cut_h_CHWBtil_NP", "jj_dphi_cut_h_CHWBtil_NP2",   1, "plots/rivet/Dphijj_CHWBtil_ratio.png", 0.5, -M_PI, M_PI, "#Delta#phi_{jj} [ rad ]", "#frac{d#sigma^{EFT}}{d#sigma^{SM}}", "C_{H#tilde{W}B} = 1");
    //theta_star
        drawRatio("output/rebin_out/ats.root", "output/rebin_out/ats.root",    "output/rebin_out/ats.root",    "Angle_truthBorn_theta_star_cut_h", "Angle_theta_star_cut_h_CW_NP", "Angle_theta_star_cut_h_CW_NP2",     1, "plots/rivet/ats_CW_ratio.png",     0.2, 0, M_PI);
        //CWtil
        drawRatio("output/rebin_out/ats.root", "output/rebin_out/ats.root",    "output/rebin_out/ats.root",    "Angle_truthBorn_theta_star_cut_h", "Angle_theta_star_cut_h_CWtil_NP", "Angle_theta_star_cut_h_CWtil_NP2",     1, "plots/rivet/ats_CWtil_ratio.png",  0.2, 0, M_PI);
        //CHW 
        drawRatio("output/rebin_out/ats.root", "output/rebin_out/ats.root",    "output/rebin_out/ats.root",    "Angle_truthBorn_theta_star_cut_h", "Angle_theta_star_cut_h_CHW_NP", "Angle_theta_star_cut_h_CHW_NP2",     1, "plots/rivet/ats_CHW_ratio.png",    0.2, 0, M_PI);
        //CHWTIL
        drawRatio("output/rebin_out/ats.root", "output/rebin_out/ats.root",    "output/rebin_out/ats.root",    "Angle_truthBorn_theta_star_cut_h", "Angle_theta_star_cut_h_CHWtil_NP", "Angle_theta_star_cut_h_CHWtil_NP2",     1, "plots/rivet/ats_CHWtil_ratio.png", 0.2, 0, M_PI);

    //theta1
        drawRatio("output/rebin_out/at1.root", "output/rebin_out/at1.root",    "output/rebin_out/at1.root",    "Angle_truthBorn_theta1_cut_h", "Angle_theta1_cut_h_CW_NP", "Angle_theta1_cut_h_CW_NP2",     1, "plots/rivet/at1_CW_ratio.png",     0.3, 0, M_PI, "#theta_{1} [ rad ]", "#frac{d#sigma^{EFT}}{d#sigma^{SM}}", "C_{W} = 1");
        //CWtil
        drawRatio("output/rebin_out/at1.root", "output/rebin_out/at1.root",    "output/rebin_out/at1.root",    "Angle_truthBorn_theta1_cut_h", "Angle_theta1_cut_h_CWtil_NP", "Angle_theta1_cut_h_CWtil_NP2",     1, "plots/rivet/at1_CWtil_ratio.png",  0.3, 0, M_PI, "#theta_{1} [ rad ]", "#frac{d#sigma^{EFT}}{d#sigma^{SM}}", "C_{#tilde{W}} = 1");
        //CHW 
        drawRatio("output/rebin_out/at1.root", "output/rebin_out/at1.root",    "output/rebin_out/at1.root",    "Angle_truthBorn_theta1_cut_h", "Angle_theta1_cut_h_CHW_NP", "Angle_theta1_cut_h_CHW_NP2",     1, "plots/rivet/at1_CHW_ratio.png",    0.3, 0, M_PI, "#theta_{1} [ rad ]", "#frac{d#sigma^{EFT}}{d#sigma^{SM}}", "C_{HW} = 1");
        //CHWTIL
        drawRatio("output/rebin_out/at1.root", "output/rebin_out/at1.root",    "output/rebin_out/at1.root",    "Angle_truthBorn_theta1_cut_h", "Angle_theta1_cut_h_CHWtil_NP", "Angle_theta1_cut_h_CHWtil_NP2",     1, "plots/rivet/at1_CHWtil_ratio.png", 0.3, 0, M_PI, "#theta_{1} [ rad ]", "#frac{d#sigma^{EFT}}{d#sigma^{SM}}", "C_{H#tilde{W}} = 1");
        //CHWB 
        drawRatio("output/rebin_out/at1.root", "output/rebin_out/at1.root",    "output/rebin_out/at1.root",    "Angle_truthBorn_theta1_cut_h", "Angle_theta1_cut_h_CHWB_NP", "Angle_theta1_cut_h_CHWB_NP2",     1, "plots/rivet/at1_CHWB_ratio.png",    0.3, 0, M_PI, "#theta_{1} [ rad ]", "#frac{d#sigma^{EFT}}{d#sigma^{SM}}", "C_{HWB} = 1");
        //CHWBTIL
        drawRatio("output/rebin_out/at1.root", "output/rebin_out/at1.root",    "output/rebin_out/at1.root",    "Angle_truthBorn_theta1_cut_h", "Angle_theta1_cut_h_CHWBtil_NP", "Angle_theta1_cut_h_CHWBtil_NP2",     1, "plots/rivet/at1_CHWBtil_ratio.png", 0.3, 0, M_PI, "#theta_{1} [ rad ]", "#frac{d#sigma^{EFT}}{d#sigma^{SM}}", "C_{H#tilde{W}B} = 1");

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