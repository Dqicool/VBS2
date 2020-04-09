#include<libs/genAna.h>
#define EPS 1e-9


double getChi2(TH2D* hcov, TH1D* hpred, TH1D* hdata){
    std::vector<double> stdv_pred, stdv_data, stdv_cov;
    for (int i = 1; i <= hpred->GetNbinsX(); i++){
        stdv_pred.push_back((double)(hpred->GetBinContent(i)));
        stdv_data.push_back((double)(hdata->GetBinContent(i)));
        
    }
    for (int i = 1; i <= hcov->GetNbinsX(); i++){
        for (int j = 1; j <= hcov->GetNbinsY(); j++){
            stdv_cov.push_back((double)hcov->GetBinContent(i,j));
        }
    }
    TMatrixD m_cov;
    TMatrixD v_pred, v_data;
    v_pred.Use(stdv_pred.size(),1, stdv_pred.data());//row vec
    v_data.Use(stdv_data.size(),1 , stdv_data.data());
    m_cov.Use(hcov->GetNbinsX(), hcov->GetNbinsY(), stdv_cov.data());

    // TCanvas c1("c1","c1",2000,2000);
    // hcov->Draw("COLZ");
    // hcov->Draw("TEXT SAME");
    // c1.SaveAs("debug.png");
    // TCanvas c2("c2","c2",2000,2000);
    // TH2D* hcov_reco = new TH2D(m_cov);
    // hcov_reco->Draw("COLZ");
    // hcov_reco->Draw("TEXT SAME");
    // c2.SaveAs("debug2.png");

    auto tmp = (v_data - v_pred).T() * m_cov.Invert();
    auto m_chi2 = tmp * (v_data - v_pred);

    return m_chi2[0][0];
}

double getMin(TH2D* h){
    double ret = 1e20;
    for(int i = 1; i <= h->GetNbinsX();i++){
        for(int j = 1; j <= h->GetNbinsY(); j++){
            if(h->GetBinContent(i,j) < ret) ret = h->GetBinContent(i,j);
        }
    }
    return ret;
}

void limCWCWtil(){
    TFile* inunf = TFile::Open("output/unfold_out/dphijj.root", "read");
    TFile* inNP  = TFile::Open("output/rebin_out/dPhiJJ.root", "read");

    TH2D* h_cov = (TH2D*)inunf->Get("h_Cov");
    TH1D* h_SM  = (TH1D*)inunf->Get("h_jj_dphi_true");
    TH1D* h_data = (TH1D*)inunf->Get("h_jj_dphi_reco");

    TH1D* h_CWtil_NP = (TH1D*)inNP->Get("jj_dphi_cut_h_CWtil_NP");
    TH1D* h_CWtil_NP2 = (TH1D*)inNP->Get("jj_dphi_cut_h_CWtil_NP2");

    TH1D* h_CW_NP = (TH1D*)inNP->Get("jj_dphi_cut_h_CW_NP");
    TH1D* h_CW_NP2 = (TH1D*)inNP->Get("jj_dphi_cut_h_CW_NP2");


    TH2D* h_chi2 = new TH2D("chi2_hist", "",76,-0.01,1.51, 76, -0.01, 1.51);
    for (double i = 0.00; i<=1.5+EPS; i+=0.02){
        for (double j = 0.00; j<=1.5+EPS; j+=0.02){
            auto h_loop_CWtil_NP = (TH1D*)(h_CWtil_NP->Clone());
            auto h_loop_CWtil_NP2 = (TH1D*)(h_CWtil_NP2->Clone());

            auto h_loop_CW_NP = (TH1D*)(h_CW_NP->Clone());
            auto h_loop_CW_NP2 = (TH1D*)(h_CW_NP2->Clone());

            auto h_loop_SM = (TH1D*)(h_SM->Clone());
            auto h_loop_data = (TH1D*)(h_data->Clone());

            h_loop_CWtil_NP->Scale(std::pow(i/5.,2));
            h_loop_CWtil_NP2->Scale(i/5.);
            h_loop_CW_NP->Scale(std::pow(j/5.,2));
            h_loop_CW_NP2->Scale(j/5.);

            h_loop_SM->Add(h_loop_CWtil_NP);
            h_loop_SM->Add(h_loop_CWtil_NP2);
            h_loop_SM->Add(h_loop_CW_NP);
            h_loop_SM->Add(h_loop_CW_NP2);
            h_chi2->Fill(i,j, getChi2(h_cov, h_loop_SM, h_data));
        }
    }

    TH2D* h_lim  = (TH2D*)h_chi2->Clone();
    double mini = getMin(h_chi2);
    for(int i = 1; i <= h_chi2->GetNbinsX();i++){
        for(int j = 1; j <= h_chi2->GetNbinsY(); j++){
            h_chi2->SetBinContent(i, j, h_chi2->GetBinContent(i,j) - mini);
            h_lim->SetBinContent(i,j, (double)(h_chi2->GetBinContent(i,j) < 2.705543454));
        }
    }
    
    TH2D* h_bond  = (TH2D*)h_lim->Clone();
    for(int i = 1; i <= h_lim->GetNbinsX();i++){
        for(int j = 1; j <= h_lim->GetNbinsY(); j++){
            h_bond->SetBinContent(i,j,0);
            if(h_lim->GetBinContent(i,j) == 1){
                if(i > 1 && j > 1 && h_lim->Integral(i-1,i+1,j-1,j+1) < 9){ h_bond->SetBinContent(i,j,1); }
                else if(((i == 1 && j != 1) || (i!=1 && j == 1)) && h_lim->Integral(i-1,i+1,j-1,j+1) < 6){ h_bond->SetBinContent(i,j,1); }
                else if(i==1 && j==1 && h_lim->Integral(i-1,i+1,j-1,j+1) < 4){ h_bond->SetBinContent(i,j,1); }
            }
        }
    }

    TCanvas c0("c0","c0",2000,2000);
    gStyle->SetPalette();
    h_chi2->SetStats(0);
    h_chi2->SetContour(64);
    h_chi2->Draw("COLZ0");
    c0.SaveAs("chi2cw.png");

    TCanvas c00("c00","c00",2000,2000);
    h_lim->SetStats(0);
    h_lim->Draw("COLZ");
    c00.SaveAs("limitcw.png");

    TCanvas c000("c000","c000",2000,2000);
    h_lim->SetStats(0);
    gStyle->SetPalette(kCherry);
    h_bond->SetAxisRange(0,2,"Z");
    h_bond->Draw("COLZ");
    c000.SaveAs("bondcw.png");

    TFile* out = TFile::Open("output/limit_out/CW.root", "recreate");
    h_chi2->Write();
    h_lim->Write();
    h_bond->Write();
    out->Close();
}

void limCHWCHWtil(){
    TFile* inunf = TFile::Open("output/unfold_out/dphijj.root", "read");
    TFile* inNP  = TFile::Open("output/rebin_out/dPhiJJ.root", "read");

    TH2D* h_cov = (TH2D*)inunf->Get("h_Cov");
    TH1D* h_SM  = (TH1D*)inunf->Get("h_jj_dphi_true");
    TH1D* h_data = (TH1D*)inunf->Get("h_jj_dphi_reco");

    TH1D* h_CHWtil_NP = (TH1D*)inNP->Get("jj_dphi_cut_h_CHWtil_NP");
    TH1D* h_CHWtil_NP2 = (TH1D*)inNP->Get("jj_dphi_cut_h_CHWtil_NP2");

    TH1D* h_CHW_NP = (TH1D*)inNP->Get("jj_dphi_cut_h_CHW_NP");
    TH1D* h_CHW_NP2 = (TH1D*)inNP->Get("jj_dphi_cut_h_CHW_NP2");


    TH2D* h_chi2 = new TH2D("chi2_hist", "",76,-0.01,1.51, 51, -0.05, 5.05);
    
    for (double i = 0.00; i<=1.5+EPS; i+=0.02){
        for (double j = 0.00; j<=5+EPS; j+=0.1){
            auto h_loop_CHWtil_NP = (TH1D*)(h_CHWtil_NP->Clone());
            auto h_loop_CHWtil_NP2 = (TH1D*)(h_CHWtil_NP2->Clone());

            auto h_loop_CHW_NP = (TH1D*)(h_CHW_NP->Clone());
            auto h_loop_CHW_NP2 = (TH1D*)(h_CHW_NP2->Clone());

            auto h_loop_SM = (TH1D*)(h_SM->Clone());
            auto h_loop_data = (TH1D*)(h_data->Clone());

            h_loop_CHWtil_NP->Scale(std::pow(i/5.,2));
            h_loop_CHWtil_NP2->Scale(i/5.);
            h_loop_CHW_NP->Scale(std::pow(j/5.,2));
            h_loop_CHW_NP2->Scale(j/5.);

            h_loop_SM->Add(h_loop_CHWtil_NP);
            h_loop_SM->Add(h_loop_CHWtil_NP2);
            h_loop_SM->Add(h_loop_CHW_NP);
            h_loop_SM->Add(h_loop_CHW_NP2);

            h_chi2->Fill(i,j, getChi2(h_cov, h_loop_SM, h_data));
        }
    }
    TH2D* h_lim  = (TH2D*)h_chi2->Clone();
    double mini = getMin(h_chi2);
    for(int i = 1; i <= h_chi2->GetNbinsX();i++){
        for(int j = 1; j <= h_chi2->GetNbinsY(); j++){
            h_chi2->SetBinContent(i, j, h_chi2->GetBinContent(i,j) - mini);
            h_lim->SetBinContent(i,j, (double)(h_chi2->GetBinContent(i,j) < 2.705543454));
        }
    }
    TH2D* h_bond  = (TH2D*)h_lim->Clone();
    for(int i = 1; i <= h_lim->GetNbinsX();i++){
        for(int j = 1; j <= h_lim->GetNbinsY(); j++){
            h_bond->SetBinContent(i,j,0);
            if(h_lim->GetBinContent(i,j) == 1){
                if(i > 1 && j > 1 && h_lim->Integral(i-1,i+1,j-1,j+1) < 9){ h_bond->SetBinContent(i,j,1); }
                else if(((i == 1 && j != 1) || (i!=1 && j == 1)) && h_lim->Integral(i-1,i+1,j-1,j+1) < 6){ h_bond->SetBinContent(i,j,1); }
                else if(i==1 && j==1 && h_lim->Integral(i-1,i+1,j-1,j+1) < 4){ h_bond->SetBinContent(i,j,1); }
            }
        }
    }
    TCanvas c0("c0","c0",2000,2000);
    gStyle->SetPalette();
    h_chi2->SetStats(0);
    h_chi2->SetContour(64);
    h_chi2->Draw("COLZ0");
    c0.SaveAs("chi2chw.png");

    TCanvas c00("c00","c00",2000,2000);
    h_lim->SetStats(0);
    h_lim->Draw("COLZ");
    c00.SaveAs("limitchw.png");

    TCanvas c000("c000","c000",2000,2000);
    h_lim->SetStats(0);
    gStyle->SetPalette(kCherry);
    h_bond->SetAxisRange(0,2,"Z");
    h_bond->Draw("COLZ");
    c000.SaveAs("bondchw.png");

    TFile* out = TFile::Open("output/limit_out/CHW.root", "recreate");
    h_chi2->Write();
    h_bond->Write();
    h_lim->Write();
    out->Close();

}

int main(){
    limCHWCHWtil();
    limCWCWtil();
}