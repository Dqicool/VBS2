#include"libs/genAna.h"
#include<TLegend.h>
#include<TLine.h>

void plotLim(   const char* fn_limo, const char* fn_lime, 
                const char* hn_lim1o, const char* hn_lim2o,
                const char* hn_lim1e, const char* hn_lim2e,
                const char* hn_asy,
                const char* xtitle, const char* ytitle, const char* fn_img){
    TFile* ine = TFile::Open(fn_lime);
    TFile* ino = TFile::Open(fn_limo);

    TH2* h_lim1e;
    ine->GetObject(hn_lim1e, h_lim1e);
    TH2* h_lim2e;
    ine->GetObject(hn_lim2e, h_lim2e);
    TH2* h_lim1o;
    ino->GetObject(hn_lim1o, h_lim1o);
    TH2* h_lim2o;
    ino->GetObject(hn_lim2o, h_lim2o);
	TH2* h_asy;
	if(((std::string)fn_limo).find((std::string)"phi") != std::string::npos){
		ino->GetObject(hn_asy, h_asy);
	}


    TCanvas c1("c1","c1",2000,2000);
    c1.SetMargin(0.15,0.15,0.15,0.15);
    h_lim1e->SetStats(0);
    h_lim1e->SetAxisRange(0,2,"Z");
    h_lim1e->SetContour(3);
    h_lim1e->GetXaxis()->SetTitle(xtitle);
    h_lim1e->GetYaxis()->SetTitle(ytitle);
    h_lim1e->SetLineStyle(kDotted);
    h_lim1e->SetLineWidth(5);
    h_lim1e->SetLineColor(kBlue);

    h_lim2e->SetStats(0);
    h_lim2e->SetAxisRange(0,2,"Z");
    h_lim2e->SetContour(3);
    h_lim2e->GetXaxis()->SetTitle(xtitle);
    h_lim2e->GetYaxis()->SetTitle(ytitle);
    //h_lim2->SetLineStyle(kDotted);
    h_lim2e->SetLineWidth(5);
    h_lim2e->SetLineColor(kBlue);

    h_lim1o->SetStats(0);
    h_lim1o->SetAxisRange(0,2,"Z");
    h_lim1o->SetContour(3);
    h_lim1o->GetXaxis()->SetTitle(xtitle);
    h_lim1o->GetYaxis()->SetTitle(ytitle);
    h_lim1o->SetLineStyle(kDotted);
    h_lim1o->SetLineWidth(5);
    h_lim1o->SetLineColor(kRed);

    h_lim2o->SetStats(0);
    h_lim2o->SetAxisRange(0,2,"Z");
    h_lim2o->SetContour(3);
    //h_lim2->SetLineStyle(kDotted);
    h_lim2o->SetLineWidth(5);
    h_lim2o->SetLineColor(kRed);

    

    h_lim2e->Draw("CONT3");
    h_lim1e->Draw("CONT3 SAME");
    h_lim2o->Draw("CONT3 SAME");
    h_lim1o->Draw("CONT3 SAME");
	if(((std::string)fn_limo).find((std::string)"phi") != std::string::npos){
		h_asy->SetStats(0);
		h_asy->GetZaxis()->SetTitle("#frac{#sigma_{+} - #sigma_{-}}{{#sigma_{+} - #sigma_{-}}");
		h_asy->SetAxisRange(-0.5,0.5,"Z");
		h_asy->Draw("COLZ SAME");
		h_lim2e->Draw("CONT3 SAME");
		h_lim1e->Draw("CONT3 SAME");
		h_lim2o->Draw("CONT3 SAME");
		h_lim1o->Draw("CONT3 SAME");
	}


    TLegend legend(0.6,0.70,0.84,0.84);
    legend.AddEntry(h_lim1o, "68% observed limit", "l");
    legend.AddEntry(h_lim2o, "95% observed limit", "l");
    legend.AddEntry(h_lim1e, "68% expected limit", "l");
    legend.AddEntry(h_lim2e, "95% expected limit", "l");
    legend.SetLineColorAlpha(kBlack, 0);
    legend.SetFillColorAlpha(kBlack, 0);
    legend.SetTextFont(22);
    legend.Draw();

	double xmin = h_lim2o->GetXaxis()->GetBinLowEdge(1);
	double xmax = h_lim2o->GetXaxis()->GetBinLowEdge(h_lim2o->GetXaxis()->GetNbins() + 1);
	double ymin = h_lim2o->GetYaxis()->GetBinLowEdge(1);
	double ymax = h_lim2o->GetYaxis()->GetBinLowEdge(h_lim2o->GetYaxis()->GetNbins() + 1);
	TLine horiline(xmin, 0, xmax, 0);
	TLine vertline(0, ymin, 0, ymax);
    horiline.SetLineWidth(1);
    horiline.Draw();
	vertline.SetLineWidth(1);
	vertline.Draw();

	TLatex texx;
    texx.SetTextAlign(11);
    texx.SetTextSize(0.03);
    texx.DrawLatex(0.02 * (xmax - xmin), 0.02 * (ymax - ymin), "#color[2]{SM}");

    c1.SaveAs(fn_img);
}

int main(){
    plotLim("output/limit_out/CHW_CHWB_at1.root", "output/limit_out/e_CHW_CHWB_at1.root", 
            "h_lim1", "h_lim2",
            "h_lim1", "h_lim2",
			"asy_hist",
            "C_{HW}", "C_{HWB}", "plots/limit/CHW_CHWB_at1.png");

    plotLim("output/limit_out/CW_CHW_at1.root", "output/limit_out/e_CW_CHW_at1.root", 
            "h_lim1", "h_lim2",
            "h_lim1", "h_lim2",
			"asy_hist",
            "C_{W}", "C_{HW}", "plots/limit/CW_CHW_at1.png");

    plotLim("output/limit_out/CW_CHWB_at1.root", "output/limit_out/e_CW_CHWB_at1.root", 
            "h_lim1", "h_lim2",
            "h_lim1", "h_lim2",
			"asy_hist",
            "C_{W}", "C_{HWB}", "plots/limit/CW_CHWB_at1.png");

    
    plotLim("output/limit_out/CWt_CHWt_phi.root", "output/limit_out/e_CWt_CHWt_phi.root", 
            "h_lim1", "h_lim2",
            "h_lim1", "h_lim2",
			"asy_hist",
            "C_{#tilde{W}}", "C_{H#tilde{W}}", "plots/limit/CWt_CHWt_phi.png");

    plotLim("output/limit_out/CWt_CHWBt_phi.root", "output/limit_out/e_CWt_CHWBt_phi.root", 
            "h_lim1", "h_lim2",
            "h_lim1", "h_lim2",
			"asy_hist",
            "C_{#tilde{W}}", "C_{H#tilde{W}B}", "plots/limit/CWt_CHWBt_phi.png");
            
    plotLim("output/limit_out/CHWt_CHWBt_phi.root", "output/limit_out/e_CHWt_CHWBt_phi.root", 
            "h_lim1", "h_lim2",
            "h_lim1", "h_lim2",
			"asy_hist",
            "C_{H#tilde{W}}", "C_{H#tilde{W}B}", "plots/limit/CHWt_CHWBt_phi.png");
}