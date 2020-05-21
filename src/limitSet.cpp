#include<libs/genAna.h>
#include<TLegend.h>
#include<TLine.h>
#include<libs/ErrorProp.h>
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

double getMin1d(TH1D* h){
    double ret = 1e20;
    for(int i = 1; i <= h->GetNbinsX();i++){
            if(h->GetBinContent(i) < ret) ret = h->GetBinContent(i);
    }
    return ret;
}

void     limSet(const char* fn_unf, const char* fn_np, 
                const char* hn_cov, const char* hn_sm, const char* hn_reco,
                const char* hn_c1_np, const char* hn_c1_np2, double c1min, int c1n, double c1max,
                const char* hn_c2_np, const char* hn_c2_np2, double c2min, int c2n, double c2max,
                uint ndf, double CL1, double CL2,
                const char* out_chi2, const char* out_lim, const char* out_root,
                const char* xtitle = "", const char* ytitle = "", const char* ztitle="",
                const char* htitle=""){

    double lim1 = ROOT::MathMore::chisquared_quantile(CL1, ndf);
    double lim2 = ROOT::MathMore::chisquared_quantile(CL2, ndf);
    std::cout<<lim1<<"\t"<<lim2<<std::endl;
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
    TH2D* h_asy  = new TH2D("asy_hist",  "",c1n+1,c1min - c1int/2,c1max + c1int/2, c2n+1, c2min - c2int/2, c2max + c2int/2);
    for (double i = c1min; i<=c1max+EPS; i+=c1int){
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
            if ((std::string)hn_sm == (std::string)"h_jj_dphi_true"){
                double errplus, errminus;
                double xsplus   = h_loop_SM->IntegralAndError(h_loop_SM->GetNbinsX() / 2 + 1, h_loop_SM->GetNbinsX(), errplus);
                double xsminus  = h_loop_SM->IntegralAndError(1, h_loop_SM->GetNbinsX() / 2, errminus);
                double asy      = (xsplus - xsminus) / (xsplus + xsminus);
                double errup    = ErrAPlusB(xsplus,xsminus,0);
                double errdown  = ErrAMinusB(xsplus,xsminus,0);
                double errasy   = ErrADiviB(asy, (xsplus - xsminus), (xsplus + xsminus), errup, errdown, 0);
                h_asy->Fill(i, j, asy);
            }
            

            h_loop_SM->Delete();
            h_loop_data->Delete();
            h_loop_CWtil_NP->Delete();
            h_loop_CWtil_NP2->Delete();
            h_loop_CW_NP->Delete();
            h_loop_CW_NP2->Delete();
        }
    }

    TH2D* h_lim1  = (TH2D*)h_chi2->Clone();
    TH2D* h_lim2  = (TH2D*)h_chi2->Clone();
    double mini = getMin(h_chi2);
    for(int i = 1; i <= h_chi2->GetNbinsX();i++){
        for(int j = 1; j <= h_chi2->GetNbinsY(); j++){
            h_chi2->SetBinContent(i, j, h_chi2->GetBinContent(i,j) - mini);
            h_lim1->SetBinContent(i,j, 0);
            h_lim2->SetBinContent(i,j, 0);
            if (h_chi2->GetBinContent(i,j) < lim1)
                h_lim1->SetBinContent(i,j, 0.35);
            if(h_chi2->GetBinContent(i,j) < lim2)
                h_lim2->SetBinContent(i,j, -0.35);
        }
    }
    
    TCanvas c0("c0","c0",2000,2000);
    c0.SetMargin(0.15,0.15,0.15,0.15);
    h_chi2->SetStats(0);
    h_chi2->SetContour(64);

    h_chi2->SetTitle(htitle);
    h_chi2->GetXaxis()->SetTitle(xtitle);
    h_chi2->GetYaxis()->SetTitle(ytitle);
    h_chi2->GetZaxis()->SetTitle(ztitle);
    h_chi2->Draw("COLZ0");
    c0.SaveAs(out_chi2);
    
    if ((std::string)hn_sm == (std::string)"h_jj_dphi_true"){
        TCanvas c1("c1","c1",2000,2000);
        c1.SetMargin(0.15,0.15,0.15,0.15);
        h_asy->SetStats(0);
        h_asy->SetContour(64);

        h_asy->SetTitle(htitle);
        h_asy->GetXaxis()->SetTitle(xtitle);
        h_asy->GetYaxis()->SetTitle(ytitle);
        h_asy->GetZaxis()->SetTitle(ztitle);
        h_asy->Draw("COLZ0");
        c1.SaveAs(out_lim);
    }

    TFile* out = TFile::Open(out_root,"recreate");
    h_lim1->Clone("h_lim1")->Write();
    h_lim2->Clone("h_lim2")->Write();
    h_chi2->Clone("h_chi2")->Write();
    if ((std::string)hn_sm == (std::string)"h_jj_dphi_true"){
        h_asy->Write();
    }
    out->Close();
}

void     limSet1d(const char* fn_unf, const char* fn_np, 
                const char* hn_cov, const char* hn_sm, const char* hn_reco,
                const char* hn_c_np, const char* hn_c_np2, double cmin, int cn, double cmax,
                uint ndf, double CL1, double CL2){

    double lim1 = ROOT::MathMore::chisquared_quantile(CL1, ndf);
    
    double lim2 = ROOT::MathMore::chisquared_quantile(CL2, ndf);
    //std::cout<<lim1<<"\t"<<lim2<<std::endl;
    double cint = (cmax - cmin) / cn;
    TFile* inunf = TFile::Open(fn_unf, "read");
    TFile* inNP  = TFile::Open(fn_np, "read");

    TH2D* h_cov = (TH2D*)inunf->Get(hn_cov);
    TH1D* h_SM  = (TH1D*)inunf->Get(hn_sm);
    TH1D* h_data = (TH1D*)inunf->Get(hn_reco);

    TH1D* h_NP = (TH1D*)inNP->Get(hn_c_np);
    TH1D* h_NP2 = (TH1D*)inNP->Get(hn_c_np2);

    TH1D* h_chi2 = new TH1D("chi2_hist", "",cn+1,cmin - cint/2,cmax + cint/2);
    TH1D* h_asy  = new TH1D("asy_hist",  "",cn+1,cmin - cint/2,cmax + cint/2);
    for (double i = cmin; i<=cmax+EPS; i+=cint){
        //std::cout<<(i - cmin)/(cmax - cmin)*100<<"%"<<std::endl;
        auto h_loop_NP = (TH1D*)(h_NP->Clone());
        auto h_loop_NP2 = (TH1D*)(h_NP2->Clone());

        auto h_loop_SM = (TH1D*)(h_SM->Clone());
        auto h_loop_data = (TH1D*)(h_data->Clone());

        h_loop_NP->Scale(std::pow(i/5.,2));
        h_loop_NP2->Scale(i/5.);

        h_loop_SM->Add(h_loop_NP);
        h_loop_SM->Add(h_loop_NP2);
        h_loop_SM->Add(h_loop_NP);
        h_loop_SM->Add(h_loop_NP2);
        h_chi2->Fill(i, getChi2(h_cov, h_loop_SM, h_data));

        h_loop_SM->Delete();
        h_loop_data->Delete();
        h_loop_NP->Delete();
        h_loop_NP2->Delete();
    }

    TH1D* h_lim1  = (TH1D*)h_chi2->Clone();
    TH1D* h_lim2  = (TH1D*)h_chi2->Clone();
    double mini = getMin1d(h_chi2);
    for(int i = 1; i <= h_chi2->GetNbinsX();i++){
            h_chi2->SetBinContent(i, h_chi2->GetBinContent(i) - mini);
            h_lim1->SetBinContent(i, 0);
            h_lim2->SetBinContent(i, 0);
            if (h_chi2->GetBinContent(i) < lim1)
                h_lim1->SetBinContent(i, 0.35);
            if(h_chi2->GetBinContent(i) < lim2)
                h_lim2->SetBinContent(i, -0.35);
    }
    std::cout<<hn_c_np<<std::endl;
    for(int i = 2; i <= h_chi2->GetNbinsX();i++){
        if(h_lim1->GetBinContent(i)-h_lim1->GetBinContent(i-1) != 0){
            std::cout<<"lim1\t"<<cmin+(cmax - cmin)/cn*i<<std::endl;
        }
    }

     for(int i = 2; i <= h_chi2->GetNbinsX();i++){
        if(h_lim2->GetBinContent(i)-h_lim2->GetBinContent(i-1) != 0){
            std::cout<<"lim2\t"<<cmin+(cmax - cmin)/cn*i<<std::endl;
        }
    }
    TCanvas c1("deby", "deby", 2000,2000);
        h_chi2->Draw();
    c1.SaveAs(((std::string)hn_reco + ".png").data());


}


int main(int argc, const char** argv){
    //CWt CHWBt phi
    int reso = 1000;
    int reso1d = 10000;
    std::cout<<argv[0]<<std::endl;
    std::cout<<argv[1]<<std::endl;
    //observed limits
        if ((std::string)(argv[1]) == (std::string)"os"){
            limSet("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
            "h_Cov",                         "h_jj_dphi_true",               "h_jj_dphi_reco", 
            "jj_dphi_cut_h_CWtil_NP",        "jj_dphi_cut_h_CWtil_NP2",      -2.0, reso, 2.0,
            "jj_dphi_cut_h_CHWBtil_NP",           "jj_dphi_cut_h_CHWBtil_NP2",         -15.0, reso, 15.0,
            2, 0.68, 0.95,
            "plots/limit/CWt_CHWBt_phi_chi2.png", "plots/limit/CWt_CHWBt_phi_lim.png", "output/limit_out/CWt_CHWBt_phi.root",
            "C_{#tilde{W}}", "C_{H#tilde{W}B}", "#Delta#chi^{2}",
            "Calculated using #Delta#phi_{jj} distribution");
        }
        if ((std::string)(argv[1]) == (std::string)"oa"){
            //CHWt CHWBt phi
            limSet("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
                "h_Cov",                         "h_jj_dphi_true",               "h_jj_dphi_reco", 
                "jj_dphi_cut_h_CHWtil_NP",        "jj_dphi_cut_h_CHWtil_NP2",      -2.0, reso, 2.0,
                "jj_dphi_cut_h_CHWBtil_NP",           "jj_dphi_cut_h_CHWBtil_NP2",         -15.0, reso, 15.0,
                2, 0.68, 0.95,
                "plots/limit/CHWt_CHWBt_phi_chi2.png", "plots/limit/CHWt_CHWBt_phi_lim.png", "output/limit_out/CHWt_CHWBt_phi.root",
                "C_{H#tilde{W}}", "C_{H#tilde{W}B}", "#Delta#chi^{2}",
                "Calculated using #Delta#phi_{jj} distribution");
        }
        if ((std::string)(argv[1]) == (std::string)"ob"){
            //CWt CHWt phi
            limSet("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
                "h_Cov",                         "h_jj_dphi_true",               "h_jj_dphi_reco", 
                "jj_dphi_cut_h_CWtil_NP",        "jj_dphi_cut_h_CWtil_NP2",      -2.0, reso, 2.0,
                "jj_dphi_cut_h_CHWtil_NP",           "jj_dphi_cut_h_CHWtil_NP2",         -2.0, reso, 2.0,
                2, 0.68, 0.95,
                "plots/limit/CWt_CHWt_phi_chi2.png", "plots/limit/CWt_CHWt_phi_lim.png", "output/limit_out/CWt_CHWt_phi.root",
                "C_{#tilde{W}}", "C_{H#tilde{W}}", "#Delta#chi^{2}", 
                "Calculated using #Delta#phi_{jj} distribution");
        }
        if ((std::string)(argv[1]) == (std::string)"oc"){
            //CW CHW theta1
            limSet("output/unfold_out/at1.root", "output/rebin_out/at1.root",
                "h_Cov",                         "h_Angle_theta1_true",               "h_Angle_theta1_reco", 
                "Angle_theta1_cut_h_CW_NP",        "Angle_theta1_cut_h_CW_NP2",      -2.0, reso, 2.0,
                "Angle_theta1_cut_h_CHW_NP",           "Angle_theta1_cut_h_CHW_NP2",         -5.0, reso, 5.0,
                2, 0.68, 0.95,
                "plots/limit/CW_CHW_at1_chi2.png", "plots/limit/CW_CHW_at1_lim.png", "output/limit_out/CW_CHW_at1.root",
                "C_{W}", "C_{HW}", "#Delta#chi^{2}",
                "Calculated using #theta_{1} distribution");
        }
        if ((std::string)(argv[1]) == (std::string)"od"){
            //CW CHWB theta1
            limSet("output/unfold_out/at1.root", "output/rebin_out/at1.root",
                "h_Cov",                         "h_Angle_theta1_true",               "h_Angle_theta1_reco", 
                "Angle_theta1_cut_h_CW_NP",        "Angle_theta1_cut_h_CW_NP2",      -2.0, reso, 2.0,
                "Angle_theta1_cut_h_CHWB_NP",           "Angle_theta1_cut_h_CHWB_NP2",         -10.0, reso, 10.0,
                2, 0.68, 0.95,
                "plots/limit/CW_CHWB_at1_chi2.png", "plots/limit/CW_CHWB_at1_lim.png", "output/limit_out/CW_CHWB_at1.root",
                "C_{W}", "C_{HWB}", "#Delta#chi^{2}",
                "Calculated using #theta_{1} distribution");
        }
        if ((std::string)(argv[1]) == (std::string)"oe"){
            //CHW CHWB theta1
            limSet("output/unfold_out/at1.root", "output/rebin_out/at1.root",
                "h_Cov",                         "h_Angle_theta1_true",               "h_Angle_theta1_reco", 
                "Angle_theta1_cut_h_CHW_NP",        "Angle_theta1_cut_h_CHW_NP2",      -5.0, reso, 5.0,
                "Angle_theta1_cut_h_CHWB_NP",           "Angle_theta1_cut_h_CHWB_NP2",         -10.0, reso, 10.0,
                2, 0.68, 0.95,
                "plots/limit/CHW_CHWB_at1_chi2.png", "plots/limit/CHW_CHWB_at1_lim.png", "output/limit_out/CHW_CHWB_at1.root",
                "C_{HW}", "C_{HWB}", "#Delta#chi^{2}",
                "Calculated using #theta_{1} distribution");
        }

    //expected limits
        if ((std::string)(argv[1]) == (std::string)"es"){
            limSet("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
            "h_Cov_toy",                         "h_jj_dphi_true",               "h_jj_dphi_toy", 
            "jj_dphi_cut_h_CWtil_NP",        "jj_dphi_cut_h_CWtil_NP2",      -2.0, reso, 2.0,
            "jj_dphi_cut_h_CHWBtil_NP",           "jj_dphi_cut_h_CHWBtil_NP2",         -15.0, reso, 15.0,
            2, 0.68, 0.95,
            "plots/limit/e_CWt_CHWBt_phi_chi2.png", "plots/limit/e_CWt_CHWBt_phi_lim.png", "output/limit_out/e_CWt_CHWBt_phi.root",
            "C_{#tilde{W}}", "C_{H#tilde{W}B}", "#Delta#chi^{2}",
            "Calculated using #Delta#phi_{jj} distribution");
        }
        if ((std::string)(argv[1]) == (std::string)"ea"){
            //CHWt CHWBt phi
            limSet("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
                "h_Cov_toy",                         "h_jj_dphi_true",               "h_jj_dphi_toy", 
                "jj_dphi_cut_h_CHWtil_NP",        "jj_dphi_cut_h_CHWtil_NP2",      -2.0, reso, 2.0,
                "jj_dphi_cut_h_CHWBtil_NP",           "jj_dphi_cut_h_CHWBtil_NP2",         -15.0, reso, 15.0,
                2, 0.68, 0.95,
                "plots/limit/e_CHWt_CHWBt_phi_chi2.png", "plots/limit/e_CHWt_CHWBt_phi_lim.png", "output/limit_out/e_CHWt_CHWBt_phi.root",
                "C_{H#tilde{W}}", "C_{H#tilde{W}B}", "#Delta#chi^{2}",
                "Calculated using #Delta#phi_{jj} distribution");
        }
        if ((std::string)(argv[1]) == (std::string)"eb"){
            //CWt CHWt phi
            limSet("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
                "h_Cov_toy",                         "h_jj_dphi_true",               "h_jj_dphi_toy", 
                "jj_dphi_cut_h_CWtil_NP",        "jj_dphi_cut_h_CWtil_NP2",      -2.0, reso, 2.0,
                "jj_dphi_cut_h_CHWtil_NP",           "jj_dphi_cut_h_CHWtil_NP2",         -2.0, reso, 2.0,
                2, 0.68, 0.95,
                "plots/limit/e_CWt_CHWt_phi_chi2.png", "plots/limit/e_CWt_CHWt_phi_lim.png", "output/limit_out/e_CWt_CHWt_phi.root",
                "C_{#tilde{W}}", "C_{H#tilde{W}}", "#Delta#chi^{2}", 
                "Calculated using #Delta#phi_{jj} distribution");
        }
        if ((std::string)(argv[1]) == (std::string)"ec"){
            //CW CHW theta1
            limSet("output/unfold_out/at1.root", "output/rebin_out/at1.root",
                "h_Cov_toy",                         "h_Angle_theta1_true",               "h_Angle_theta1_toy", 
                "Angle_theta1_cut_h_CW_NP",        "Angle_theta1_cut_h_CW_NP2",      -2.0, reso, 2.0,
                "Angle_theta1_cut_h_CHW_NP",           "Angle_theta1_cut_h_CHW_NP2",         -5.0, reso, 5.0,
                2, 0.68, 0.95,
                "plots/limit/e_CW_CHW_at1_chi2.png", "plots/limit/e_CW_CHW_at1_lim.png", "output/limit_out/e_CW_CHW_at1.root",
                "C_{W}", "C_{HW}", "#Delta#chi^{2}",
                "Calculated using #theta_{1} distribution");
        }
        if ((std::string)(argv[1]) == (std::string)"ed"){
            //CW CHWB theta1
            limSet("output/unfold_out/at1.root", "output/rebin_out/at1.root",
                "h_Cov_toy",                         "h_Angle_theta1_true",               "h_Angle_theta1_toy", 
                "Angle_theta1_cut_h_CW_NP",        "Angle_theta1_cut_h_CW_NP2",      -2.0, reso, 2.0,
                "Angle_theta1_cut_h_CHWB_NP",           "Angle_theta1_cut_h_CHWB_NP2",         -10.0, reso, 10.0,
                2, 0.68, 0.95,
                "plots/limit/e_CW_CHWB_at1_chi2.png", "plots/limit/e_CW_CHWB_at1_lim.png", "output/limit_out/e_CW_CHWB_at1.root",
                "C_{W}", "C_{HWB}", "#Delta#chi^{2}",
                "Calculated using #theta_{1} distribution");
        }
        if ((std::string)(argv[1]) == (std::string)"ee"){
            //CHW CHWB theta1
            limSet("output/unfold_out/at1.root", "output/rebin_out/at1.root",
                "h_Cov_toy",                         "h_Angle_theta1_true",               "h_Angle_theta1_toy", 
                "Angle_theta1_cut_h_CHW_NP",        "Angle_theta1_cut_h_CHW_NP2",      -5.0, reso, 5.0,
                "Angle_theta1_cut_h_CHWB_NP",           "Angle_theta1_cut_h_CHWB_NP2",         -10.0, reso, 10.0,
                2, 0.68, 0.95,
                "plots/limit/e_CHW_CHWB_at1_chi2.png", "plots/limit/e_CHW_CHWB_at1_lim.png", "output/limit_out/e_CHW_CHWB_at1.root",
                "C_{HW}", "C_{HWB}", "#Delta#chi^{2}",
                "Calculated using #theta_{1} distribution");
        }
        if ((std::string)(argv[1]) == (std::string)"1d"){
        // std::cout<<"--------CHW----------"<<std::endl;
        //     limSet1d("output/unfold_out/at1.root", "output/rebin_out/at1.root",
        //         "h_Cov",                         "h_Angle_theta1_true",               "h_Angle_theta1_reco", 
        //         "Angle_theta1_cut_h_CHW_NP",        "Angle_theta1_cut_h_CHW_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
        //         "h_Cov",                         "h_jj_dphi_true",               "h_jj_dphi_reco", 
        //         "jj_dphi_cut_h_CHW_NP",        "jj_dphi_cut_h_CHW_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/at1.root", "output/rebin_out/at1.root",
        //         "h_Cov_toy",                         "h_Angle_theta1_true",               "h_Angle_theta1_toy", 
        //         "Angle_theta1_cut_h_CHW_NP",        "Angle_theta1_cut_h_CHW_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
        //         "h_Cov_toy",                         "h_jj_dphi_true",               "h_jj_dphi_toy", 
        //         "jj_dphi_cut_h_CHW_NP",        "jj_dphi_cut_h_CHW_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        // std::cout<<"---------------------"<<std::endl;

        // std::cout<<"---------CW----------"<<std::endl;
        //     limSet1d("output/unfold_out/at1.root", "output/rebin_out/at1.root",
        //         "h_Cov",                         "h_Angle_theta1_true",               "h_Angle_theta1_reco", 
        //         "Angle_theta1_cut_h_CW_NP",        "Angle_theta1_cut_h_CW_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
        //         "h_Cov",                         "h_jj_dphi_true",               "h_jj_dphi_reco", 
        //         "jj_dphi_cut_h_CW_NP",        "jj_dphi_cut_h_CW_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/at1.root", "output/rebin_out/at1.root",
        //         "h_Cov_toy",                         "h_Angle_theta1_true",               "h_Angle_theta1_toy", 
        //         "Angle_theta1_cut_h_CW_NP",        "Angle_theta1_cut_h_CW_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
        //         "h_Cov_toy",                         "h_jj_dphi_true",               "h_jj_dphi_toy", 
        //         "jj_dphi_cut_h_CW_NP",        "jj_dphi_cut_h_CW_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        // std::cout<<"---------------------"<<std::endl;

        std::cout<<"--------CHWB----------"<<std::endl;
            limSet1d("output/unfold_out/at1.root", "output/rebin_out/at1.root",
                "h_Cov",                         "h_Angle_theta1_true",               "h_Angle_theta1_reco", 
                "Angle_theta1_cut_h_CHWB_NP",        "Angle_theta1_cut_h_CHWB_NP2",      -20.0, reso1d, 20.0,
                1, 0.68, 0.95);
            limSet1d("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
                "h_Cov",                         "h_jj_dphi_true",               "h_jj_dphi_reco", 
                "jj_dphi_cut_h_CHWB_NP",        "jj_dphi_cut_h_CHWB_NP2",      -20.0, reso1d, 20.0,
                1, 0.68, 0.95);
            limSet1d("output/unfold_out/at1.root", "output/rebin_out/at1.root",
                "h_Cov_toy",                         "h_Angle_theta1_true",               "h_Angle_theta1_toy", 
                "Angle_theta1_cut_h_CHWB_NP",        "Angle_theta1_cut_h_CHWB_NP2",      -20.0, reso1d, 20.0,
                1, 0.68, 0.95);
            limSet1d("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
                "h_Cov_toy",                         "h_jj_dphi_true",               "h_jj_dphi_toy", 
                "jj_dphi_cut_h_CHWB_NP",        "jj_dphi_cut_h_CHWB_NP2",      -20.0, reso1d, 20.0,
                1, 0.68, 0.95);
        std::cout<<"---------------------"<<std::endl;
            
        // std::cout<<"---------CHWBt---------"<<std::endl;
        //     limSet1d("output/unfold_out/at1.root", "output/rebin_out/at1.root",
        //         "h_Cov",                         "h_Angle_theta1_true",               "h_Angle_theta1_reco", 
        //         "Angle_theta1_cut_h_CHWBtil_NP",        "Angle_theta1_cut_h_CHWBtil_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
        //         "h_Cov",                         "h_jj_dphi_true",               "h_jj_dphi_reco", 
        //         "jj_dphi_cut_h_CHWBtil_NP",        "jj_dphi_cut_h_CHWBtil_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/at1.root", "output/rebin_out/at1.root",
        //         "h_Cov_toy",                         "h_Angle_theta1_true",               "h_Angle_theta1_toy", 
        //         "Angle_theta1_cut_h_CHWBtil_NP",        "Angle_theta1_cut_h_CHWBtil_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
        //         "h_Cov_toy",                         "h_jj_dphi_true",               "h_jj_dphi_toy", 
        //         "jj_dphi_cut_h_CHWBtil_NP",        "jj_dphi_cut_h_CHWBtil_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        // std::cout<<"---------------------"<<std::endl;

        // std::cout<<"---------CHWt----------"<<std::endl;
        //     limSet1d("output/unfold_out/at1.root", "output/rebin_out/at1.root",
        //         "h_Cov",                         "h_Angle_theta1_true",               "h_Angle_theta1_reco", 
        //         "Angle_theta1_cut_h_CHWtil_NP",        "Angle_theta1_cut_h_CHWtil_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
        //         "h_Cov",                         "h_jj_dphi_true",               "h_jj_dphi_reco", 
        //         "jj_dphi_cut_h_CHWtil_NP",        "jj_dphi_cut_h_CHWtil_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/at1.root", "output/rebin_out/at1.root",
        //         "h_Cov_toy",                         "h_Angle_theta1_true",               "h_Angle_theta1_toy", 
        //         "Angle_theta1_cut_h_CHWtil_NP",        "Angle_theta1_cut_h_CHWtil_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
        //         "h_Cov_toy",                         "h_jj_dphi_true",               "h_jj_dphi_toy", 
        //         "jj_dphi_cut_h_CHWtil_NP",        "jj_dphi_cut_h_CHWtil_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        // std::cout<<"-----------------------"<<std::endl;

        // std::cout<<"---------CWt----------"<<std::endl;
        //     limSet1d("output/unfold_out/at1.root", "output/rebin_out/at1.root",
        //         "h_Cov",                         "h_Angle_theta1_true",               "h_Angle_theta1_reco", 
        //         "Angle_theta1_cut_h_CWtil_NP",        "Angle_theta1_cut_h_CWtil_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
        //         "h_Cov",                         "h_jj_dphi_true",               "h_jj_dphi_reco", 
        //         "jj_dphi_cut_h_CWtil_NP",        "jj_dphi_cut_h_CWtil_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/at1.root", "output/rebin_out/at1.root",
        //         "h_Cov_toy",                         "h_Angle_theta1_true",               "h_Angle_theta1_toy", 
        //         "Angle_theta1_cut_h_CWtil_NP",        "Angle_theta1_cut_h_CWtil_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        //     limSet1d("output/unfold_out/dphijj.root", "output/rebin_out/dPhiJJ.root",
        //         "h_Cov_toy",                         "h_jj_dphi_true",               "h_jj_dphi_toy", 
        //         "jj_dphi_cut_h_CWtil_NP",        "jj_dphi_cut_h_CWtil_NP2",      -10.0, reso1d, 10.0,
        //         1, 0.68, 0.95);
        // std::cout<<"----------------------"<<std::endl;
        }
}