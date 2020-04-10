#include<libs/genAna.h>
#define EPS 1e-9

enum MCType { SM, NPLinear, NPQuad };
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

void     limSet(const char* fn_unf, const char* fn_np, 
                const char* hn_cov, const char* hn_sm, const char* hn_reco,
                const char* hn_c1_np, const char* hn_c1_np2, double c1min, int c1n, double c1max,
                const char* hn_c2_np, const char* hn_c2_np2, double c2min, int c2n, double c2max,
                uint ndf, double CL1, double CL2,
                const char* out_chi2, const char* out_lim, const char* out_root){

    double lim1 = ROOT::MathMore::chisquared_quantile(CL1, ndf);
    double lim2 = ROOT::MathMore::chisquared_quantile(CL2, ndf);
    double c1int = (c1max - c1min) / c1n;
    double c2int = (c2max - c2min) / c2n;
    TFile* inunf = TFile::Open(fn_unf, "read");
    TFile* inNP  = TFile::Open(fn_np, "read");

    TH2D* h_cov = (TH2D*)inunf->Get(hn_cov);
    TH1D* h_SM  = (TH1D*)inunf->Get(hn_sm);
    TH1D* h_data = (TH1D*)inunf->Get(hn_reco);

    TH1D* h_CWtil_NP = (TH1D*)inNP->Get(hn_c1_np);
    TH1D* h_CWtil_NP2 = (TH1D*)inNP->Get(hn_c1_np2);

    TH1D* h_CW_NP = (TH1D*)inNP->Get(hn_c2_np);
    TH1D* h_CW_NP2 = (TH1D*)inNP->Get(hn_c2_np2);


    TH2D* h_chi2 = new TH2D("chi2_hist", "",c1n+1,c1min - c1int/2,c1max + c1int/2, c2n+1, c2min - c2int/2, c2max + c2int/2);
    for (double i = c1min; i<=c1max+EPS; i+=c1int){
        std::system("clear");
        std::cout<<(i - c1min)/(c1max - c1min)*100<<"%"<<std::endl;
        for (double j = c2min; j<=c2max+EPS; j+=c2int){
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
            if (h_chi2->GetBinContent(i,j) < lim1)
                h_lim->SetBinContent(i,j, 1);
            else if(h_chi2->GetBinContent(i,j) < lim2)
                h_lim->SetBinContent(i,j, 2);
            else 
                h_lim->SetBinContent(i,j, 0);
        }
    }

    TCanvas c0("c0","c0",2000,2000);
    h_chi2->SetStats(0);
    h_chi2->SetContour(64);
    h_chi2->Draw("COLZ0");
    c0.SaveAs(out_chi2);

    TCanvas c1("c1","c1",2000,2000);
    h_lim->SetStats(0);
    h_lim->SetAxisRange(0,2,"Z");
    h_lim->SetContour(3);
    h_lim->Draw("COL");
    c1.SaveAs(out_lim);

    TFile* out = TFile::Open(out_root,"recreate");
    h_lim->Write();
    h_chi2->Write();
}


int main(){
    //CWt CHWt phi
    limSet("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
           "h_Cov",                         "h_jj_dphi_true",               "h_jj_dphi_reco", 
           "jj_dphi_cut_h_CWtil_NP",        "jj_dphi_cut_h_CWtil_NP2",      -2.0, 100, 2.0,
           "jj_dphi_cut_h_CHWtil_NP",           "jj_dphi_cut_h_CHWtil_NP2",         -2.0, 100, 2.0,
           2, 0.9, 0.95,
           "plots/limit/CWt_CHWt_phi_chi2.png", "plots/limit/CWt_CHWt_phi_lim.png", "output/limit_out/CWt_CHWt_phi.root");
    //CW CHW theta1
    limSet("output/unfold_out/at1.root", "output/rebin_out/at1.root",
           "h_Cov",                         "h_Angle_theta1_true",               "h_Angle_theta1_reco", 
           "Angle_theta1_cut_h_CW_NP",        "Angle_theta1_cut_h_CW_NP2",      -2.0, 100, 2.0,
           "Angle_theta1_cut_h_CHW_NP",           "Angle_theta1_cut_h_CHW_NP2",         -5.0, 100, 5.0,
           2, 0.9, 0.95,
           "plots/limit/CW_CHW_at1_chi2.png", "plots/limit/CW_CHW_at1_lim.png", "output/limit_out/CWt_CHWt_at1.root");
    //
}