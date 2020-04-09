#include "libs/genAna.h"
#define LUMI 139e3
#include<THStack.h>
#define IMAGEX 3000
#define IMAGEY 2000

//#define debug
using namespace std;
TH1D getDraw(const char* file, const char* filterr, const char* dist,  Color_t col)
{
    auto histoname = "h_" + (string)dist + "_" + (string)filterr;
    TFile *f = TFile::Open(file, "read");
    TH1D *h = (TH1D*)f->Get(&histoname[0]);
    h->SetFillColor(col);
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    return (*h);
}


void plotPropSig(const char* outfile, const char* dist, bool truth){
    ROOT::EnableImplicitMT(24);
    TH1::SetDefaultSumw2();
    //loading data from files
        std::vector<std::vector<std::string>> files;
        std::vector<Color_t> color_vec;
        std::vector<std::string> cata;
        std::vector<std::string> color_name;
        //SIGNAL
            std::vector<std::string> signals;
            signals.push_back("output/draw_out/364364.Sherpa_222_NNPDF30NNLO_lllljj_EW6_noHiggs.root");
            signals.push_back("output/draw_out/346228.PowhegPy8EG_NNPDF30_AZNLOCTEQ6L1_VBFH125_ZZ4lep_notau.root");

        //Higgs
            std::vector<std::string> higgses;
            #ifndef debug
            // higgses.push_back("output/draw_out/345038.PowhegPythia8EvtGen_NNPDF30_AZNLO_ZH125J_Zincl_MINLO.root");
            // higgses.push_back("output/draw_out/345039.PowhegPythia8EvtGen_NNPDF30_AZNLO_WpH125J_Wincl_MINLO.root");
            // higgses.push_back("output/draw_out/345040.PowhegPythia8EvtGen_NNPDF30_AZNLO_WmH125J_Wincl_MINLO.root");
            #endif
        //gg4l
            std::vector<std::string> gg4ls;
            #ifndef debug
            gg4ls.push_back("output/draw_out/345060.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_ZZ4l.root");
            gg4ls.push_back("output/draw_out/345706.Sherpa_222_NNPDF30NNLO_ggllll_130M4l.root");
            //gg4ls.push_back("output/draw_out/345708.Sherpa_222_NNPDF30NNLO_ggllllNoHiggs_0M4l130.root");
            #endif
        //TRIBOSON
            std::vector<std::string> tribosons;
            #ifndef debug
            // tribosons.push_back("output/draw_out/364243.Sherpa_222_NNPDF30NNLO_WWZ_4l2v_EW6.root");
            // tribosons.push_back("output/draw_out/364245.Sherpa_222_NNPDF30NNLO_WZZ_5l1v_EW6.root");
            // tribosons.push_back("output/draw_out/364247.Sherpa_222_NNPDF30NNLO_ZZZ_6l0v_EW6.root");
            // tribosons.push_back("output/draw_out/364248.Sherpa_222_NNPDF30NNLO_ZZZ_4l2v_EW6.root");
            #endif
        //qq4l
            std::vector<std::string> qq4ls;
            #ifndef debug
            qq4ls.push_back("output/draw_out/364250.Sherpa_222_NNPDF30NNLO_llll.root");
            qq4ls.push_back("output/draw_out/346340.PowhegPy8EG_A14NNPDF23_NNPDF30ME_ttH125_ZZ4l_allhad.root");
            qq4ls.push_back("output/draw_out/346341.PowhegPy8EG_A14NNPDF23_NNPDF30ME_ttH125_ZZ4l_semilep.root");
            qq4ls.push_back("output/draw_out/346342.PowhegPy8EG_A14NNPDF23_NNPDF30ME_ttH125_ZZ4l_dilep.root");
            #endif
        //WZllvl
            std::vector<std::string> WZs;
            #ifndef debug
            // // WZs.push_back("output/draw_out/361601.PowhegPy8EG_WZlvll_mll4.root");
            #endif
        //rest
            std::vector<std::string> rests;
            #ifndef debug
            // rests.push_back("output/draw_out/410142.Sherpa_NNPDF30NNLO_ttll_mll5.root");
            // // rests.push_back("output/draw_out/410472.PhPy8EG_A14_ttbar_hdamp258p75_dil.root");

            // //rests.push_back("output/draw_out/344295.Sherpa_Zee_4lMassFilter40GeV8GeV.root");
            // //rests.push_back("output/draw_out/344296.Sherpa_Zmumu_4lMassFilter40GeV8GeV.root");
            // // rests.push_back("output/draw_out/344297.Sherpa_Zee_3lPtFilter4GeV_4lMassVeto40GeV8GeV.root");
            // //rests.push_back("output/draw_out/344298.Sherpa_Zmumu_3lPtFilter4GeV_4lMassVeto40GeV8GeV.root");
            #endif
        files =         {WZs,       rests,      higgses,    tribosons,      signals,   gg4ls,       qq4ls};
        cata =          {"WZs",     "rests",    "higgses",  "tribosons",    "signals", "gg4ls",     "qq4ls"};
        color_vec =     {kOrange,   kOrange,    kOrange,    kOrange,        kRed,       kBlue,   kMagenta};
        color_name =    {"kOrange", "kOrange",  "kOrange",  "kOrange",      "kRed",     "kMagenta", "kBlue"};
        
    //stack
        string cut_st_n = (string)dist + "_cut_stack";
        string cr_st_n  = (string)dist + "_nct_stack";
        string sr_st_n  = (string)dist + "_sr_stack";

        string cut_h_n = (string)dist + "_cut_h";
        string cr_h_n  = (string)dist + "_cr_h";
        string sr_h_n  = (string)dist + "_sr_h";

        THStack* cut_stack = new THStack(&cut_st_n[0], "");
        THStack* cr_stack = new THStack(&cr_st_n[0], "");
        THStack* sr_stack = new THStack(&sr_st_n[0], "");

        std::vector<TH1D> cut_histo, cr_histo, sr_histo;
        

        for(uint i = 0; i < files.size(); i++){
            for (uint j = 0; j < (files[i]).size(); j++){
                if(truth){
                    cut_histo.push_back(getDraw(&(files[i][j])[0], "pass_cut_truthBorn", dist, color_vec[i]));
                    cr_histo.push_back(getDraw(&(files[i][j])[0], "pass_CR_truthBorn", dist, color_vec[i]));
                    sr_histo.push_back(getDraw(&(files[i][j])[0], "pass_SR_truthBorn", dist, color_vec[i]));
                }
                else{
                    cut_histo.push_back(getDraw(&(files[i][j])[0], "pass_cut", dist, color_vec[i]));
                    cr_histo.push_back(getDraw(&(files[i][j])[0], "pass_CR", dist, color_vec[i]));
                    sr_histo.push_back(getDraw(&(files[i][j])[0], "pass_SR", dist, color_vec[i]));
                }
            }
            cout<<cata[i] + " is in color " + color_name[i]<<endl;
        }
        
        auto nbins = (cut_histo[0]).GetNbinsX();
        auto xmin  = (cut_histo[0]).GetBinLowEdge(1);
        auto xmax =  (cut_histo[0]).GetBinLowEdge(nbins) + (cut_histo[0]).GetBinWidth(nbins);

        TH1D * h_inc_sr = new TH1D(&sr_h_n[0],"",nbins,xmin,xmax);
        TH1D * h_inc_cr = new TH1D(&cr_h_n[0],"",nbins,xmin,xmax);
        TH1D * h_inc_cut = new TH1D(&cut_h_n[0],"",nbins,xmin,xmax);
        
        for(uint i=0; i<cut_histo.size(); i++)
        {
            cut_stack->Add(&(cut_histo[i]));
            cr_stack->Add(&(cr_histo[i]));
            sr_stack->Add(&(sr_histo[i]));

            h_inc_cut->Add(&(cut_histo[i]));
            h_inc_cr->Add(&(cr_histo[i]));
            h_inc_sr->Add(&(sr_histo[i]));
        }
    //Store
        TFile * out = TFile::Open(outfile,"recreate");
        sr_stack->Write();
        cut_stack->Write();
        cr_stack->Write();

        h_inc_cut->Write();
        h_inc_cr->Write();
        h_inc_sr->Write();
        out->Close();

    //print stat
        cout<<"det\t"<<h_inc_cr->Integral()<<"\t\t"<<(cr_histo[0]).Integral() + (cr_histo[1]).Integral()<<"\t\t"<<h_inc_cr->Integral()/((cr_histo[0]).Integral() + (cr_histo[1]).Integral())<<":1"<<endl;
        cout<<"cut\t"<<h_inc_cut->Integral()<<"\t\t"<<(cut_histo[0]).Integral() + (cut_histo[1]).Integral()<<"\t\t"<<h_inc_cut->Integral()/((cut_histo[0]).Integral() + (cut_histo[1]).Integral())<<":1"<<endl;
        cout<<"sr\t"<<h_inc_sr->Integral()<<"\t\t"<<(sr_histo[0]).Integral() + (sr_histo[1]).Integral()<<"\t\t"<<h_inc_sr->Integral()/((sr_histo[0]).Integral() + (sr_histo[1]).Integral())<<":1"<<endl;;

    //drawstack
        //sr plot
            TCanvas c1("c1","",IMAGEX,IMAGEY);
            h_inc_sr->SetMarkerColor(kBlack);
            h_inc_sr->SetLineColor(kBlack);
            h_inc_sr->SetFillColor(kBlack);
            h_inc_sr->SetFillStyle(3017);
            //h_inc_sr->SetAxisRange(0,YMAX,"Y");
            sr_stack->Draw("hist");
            h_inc_sr->Draw("E2, same");
            auto sr_save_name = "plots/stack/sr_" + (string)dist + "_stack.png";
            c1.SaveAs(&sr_save_name[0]);
        //cr plot
            TCanvas c2("c2","",IMAGEX,IMAGEY);
            h_inc_cr->SetMarkerColor(kBlack);
            h_inc_cr->SetLineColor(kBlack);
            h_inc_cr->SetFillColor(kBlack);
            h_inc_cr->SetFillStyle(3017);
            //h_inc_nn->SetAxinnange(0,YMAX,"Y");
            cr_stack->Draw("hist");
            h_inc_cr->Draw("E2,same");
            auto cr_save_name = "plots/stack/nn_" + (string)dist + "_stack.png";
            c2.SaveAs(&cr_save_name[0]);
        //cut plot
            TCanvas c3("c3","",IMAGEX,IMAGEY);
            h_inc_cut->SetMarkerColor(kBlack);
            h_inc_cut->SetLineColor(kBlack);
            h_inc_cut->SetFillColor(kBlack);
            h_inc_cut->SetFillStyle(3017);
            //h_inc_cut->SetAxicutange(0,YMAX,"Y");
            cut_stack->Draw("hist");
            h_inc_cut->Draw("E2,same");
            auto cut_save_name = "plots/stack/cut_" + (string)dist + "_stack.png";
            c3.SaveAs(&cut_save_name[0]);
        // output
        double xs, errerr;
        xs = h_inc_cut->IntegralAndError(1,1000, errerr);
        cout<<"*************************"<<endl;
        cout<<"XS = "<<xs<<"\t"<<errerr<<endl;
        cout<<"*************************"<<endl;

}

void plotPropBkg(const char* outfile, const char* dist, bool truth){
    ROOT::EnableImplicitMT(24);
    TH1::SetDefaultSumw2();
    //loading data from files
        std::vector<std::vector<std::string>> files;
        std::vector<Color_t> color_vec;
        std::vector<std::string> cata;
        std::vector<std::string> color_name;
        //SIGNAL
            std::vector<std::string> signals;
            // signals.push_back("output/draw_out/364364.Sherpa_222_NNPDF30NNLO_lllljj_EW6_noHiggs.root");
            // signals.push_back("output/draw_out/346228.PowhegPy8EG_NNPDF30_AZNLOCTEQ6L1_VBFH125_ZZ4lep_notau.root");

        //Higgs
            std::vector<std::string> higgses;
            #ifndef debug
            higgses.push_back("output/draw_out/345038.PowhegPythia8EvtGen_NNPDF30_AZNLO_ZH125J_Zincl_MINLO.root");
            higgses.push_back("output/draw_out/345039.PowhegPythia8EvtGen_NNPDF30_AZNLO_WpH125J_Wincl_MINLO.root");
            higgses.push_back("output/draw_out/345040.PowhegPythia8EvtGen_NNPDF30_AZNLO_WmH125J_Wincl_MINLO.root");
            #endif
        //gg4l
            std::vector<std::string> gg4ls;
            #ifndef debug
            // gg4ls.push_back("output/draw_out/345060.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_ZZ4l.root");
            // gg4ls.push_back("output/draw_out/345706.Sherpa_222_NNPDF30NNLO_ggllll_130M4l.root");
            // //gg4ls.push_back("output/draw_out/345708.Sherpa_222_NNPDF30NNLO_ggllllNoHiggs_0M4l130.root");
            #endif
        //TRIBOSON
            std::vector<std::string> tribosons;
            #ifndef debug
            tribosons.push_back("output/draw_out/364243.Sherpa_222_NNPDF30NNLO_WWZ_4l2v_EW6.root");
            tribosons.push_back("output/draw_out/364245.Sherpa_222_NNPDF30NNLO_WZZ_5l1v_EW6.root");
            tribosons.push_back("output/draw_out/364247.Sherpa_222_NNPDF30NNLO_ZZZ_6l0v_EW6.root");
            tribosons.push_back("output/draw_out/364248.Sherpa_222_NNPDF30NNLO_ZZZ_4l2v_EW6.root");
            #endif
        //qq4l
            std::vector<std::string> qq4ls;
            #ifndef debug
            // qq4ls.push_back("output/draw_out/364250.Sherpa_222_NNPDF30NNLO_llll.root");
            // qq4ls.push_back("output/draw_out/346340.PowhegPy8EG_A14NNPDF23_NNPDF30ME_ttH125_ZZ4l_allhad.root");
            // qq4ls.push_back("output/draw_out/346341.PowhegPy8EG_A14NNPDF23_NNPDF30ME_ttH125_ZZ4l_semilep.root");
            // qq4ls.push_back("output/draw_out/346342.PowhegPy8EG_A14NNPDF23_NNPDF30ME_ttH125_ZZ4l_dilep.root");
            #endif
        //WZllvl
            std::vector<std::string> WZs;
            #ifndef debug
            // WZs.push_back("output/draw_out/361601.PowhegPy8EG_WZlvll_mll4.root");
            #endif
        //rest
            std::vector<std::string> rests;
            #ifndef debug
            rests.push_back("output/draw_out/410142.Sherpa_NNPDF30NNLO_ttll_mll5.root");
            // rests.push_back("output/draw_out/410472.PhPy8EG_A14_ttbar_hdamp258p75_dil.root");

            //rests.push_back("output/draw_out/344295.Sherpa_Zee_4lMassFilter40GeV8GeV.root");
            //rests.push_back("output/draw_out/344296.Sherpa_Zmumu_4lMassFilter40GeV8GeV.root");
            // rests.push_back("output/draw_out/344297.Sherpa_Zee_3lPtFilter4GeV_4lMassVeto40GeV8GeV.root");
            //rests.push_back("output/draw_out/344298.Sherpa_Zmumu_3lPtFilter4GeV_4lMassVeto40GeV8GeV.root");
            #endif
        files =         {WZs,       rests,      higgses,    tribosons,      signals,   gg4ls,       qq4ls};
        cata =          {"WZs",     "rests",    "higgses",  "tribosons",    "signals", "gg4ls",     "qq4ls"};
        color_vec =     {kOrange,   kOrange,    kOrange,    kOrange,        kRed,       kBlue,   kMagenta};
        color_name =    {"kOrange", "kOrange",  "kOrange",  "kOrange",      "kRed",     "kMagenta", "kBlue"};
        
    //stack
        string cut_st_n = (string)dist + "_cut_stack";
        string cr_st_n  = (string)dist + "_nct_stack";
        string sr_st_n  = (string)dist + "_sr_stack";

        string cut_h_n = (string)dist + "_cut_h";
        string cr_h_n  = (string)dist + "_cr_h";
        string sr_h_n  = (string)dist + "_sr_h";

        THStack* cut_stack = new THStack(&cut_st_n[0], "");
        THStack* cr_stack = new THStack(&cr_st_n[0], "");
        THStack* sr_stack = new THStack(&sr_st_n[0], "");

        std::vector<TH1D> cut_histo, cr_histo, sr_histo;
        

        for(uint i = 0; i < files.size(); i++){
            for (uint j = 0; j < (files[i]).size(); j++){
                if(truth){
                    cut_histo.push_back(getDraw(&(files[i][j])[0], "pass_cut_truthBorn", dist, color_vec[i]));
                    cr_histo.push_back(getDraw(&(files[i][j])[0], "pass_CR_truthBorn", dist, color_vec[i]));
                    sr_histo.push_back(getDraw(&(files[i][j])[0], "pass_SR_truthBorn", dist, color_vec[i]));
                }
                else{
                    cut_histo.push_back(getDraw(&(files[i][j])[0], "pass_cut", dist, color_vec[i]));
                    cr_histo.push_back(getDraw(&(files[i][j])[0], "pass_CR", dist, color_vec[i]));
                    sr_histo.push_back(getDraw(&(files[i][j])[0], "pass_SR", dist, color_vec[i]));
                }
            }
            cout<<cata[i] + " is in color " + color_name[i]<<endl;
        }
        
        auto nbins = (cut_histo[0]).GetNbinsX();
        auto xmin  = (cut_histo[0]).GetBinLowEdge(1);
        auto xmax =  (cut_histo[0]).GetBinLowEdge(nbins) + (cut_histo[0]).GetBinWidth(nbins);

        TH1D * h_inc_sr = new TH1D(&sr_h_n[0],"",nbins,xmin,xmax);
        TH1D * h_inc_cr = new TH1D(&cr_h_n[0],"",nbins,xmin,xmax);
        TH1D * h_inc_cut = new TH1D(&cut_h_n[0],"",nbins,xmin,xmax);
        
        for(uint i=0; i<cut_histo.size(); i++)
        {
            cut_stack->Add(&(cut_histo[i]));
            cr_stack->Add(&(cr_histo[i]));
            sr_stack->Add(&(sr_histo[i]));

            h_inc_cut->Add(&(cut_histo[i]));
            h_inc_cr->Add(&(cr_histo[i]));
            h_inc_sr->Add(&(sr_histo[i]));
        }
    //Store
        TFile * out = TFile::Open(outfile,"recreate");
        sr_stack->Write();
        cut_stack->Write();
        cr_stack->Write();

        h_inc_cut->Write();
        h_inc_cr->Write();
        h_inc_sr->Write();
        out->Close();

    //print stat
        cout<<"det\t"<<h_inc_cr->Integral()<<"\t\t"<<(cr_histo[0]).Integral() + (cr_histo[1]).Integral()<<"\t\t"<<h_inc_cr->Integral()/((cr_histo[0]).Integral() + (cr_histo[1]).Integral())<<":1"<<endl;
        cout<<"cut\t"<<h_inc_cut->Integral()<<"\t\t"<<(cut_histo[0]).Integral() + (cut_histo[1]).Integral()<<"\t\t"<<h_inc_cut->Integral()/((cut_histo[0]).Integral() + (cut_histo[1]).Integral())<<":1"<<endl;
        cout<<"sr\t"<<h_inc_sr->Integral()<<"\t\t"<<(sr_histo[0]).Integral() + (sr_histo[1]).Integral()<<"\t\t"<<h_inc_sr->Integral()/((sr_histo[0]).Integral() + (sr_histo[1]).Integral())<<":1"<<endl;;

    //drawstack
        //sr plot
            TCanvas c1("c1","",IMAGEX,IMAGEY);
            h_inc_sr->SetMarkerColor(kBlack);
            h_inc_sr->SetLineColor(kBlack);
            h_inc_sr->SetFillColor(kBlack);
            h_inc_sr->SetFillStyle(3017);
            //h_inc_sr->SetAxisRange(0,YMAX,"Y");
            sr_stack->Draw("hist");
            h_inc_sr->Draw("E2, same");
            auto sr_save_name = "plots/stack/sr_" + (string)dist + "_stack.png";
            c1.SaveAs(&sr_save_name[0]);
        //cr plot
            TCanvas c2("c2","",IMAGEX,IMAGEY);
            h_inc_cr->SetMarkerColor(kBlack);
            h_inc_cr->SetLineColor(kBlack);
            h_inc_cr->SetFillColor(kBlack);
            h_inc_cr->SetFillStyle(3017);
            //h_inc_nn->SetAxinnange(0,YMAX,"Y");
            cr_stack->Draw("hist");
            h_inc_cr->Draw("E2,same");
            auto cr_save_name = "plots/stack/nn_" + (string)dist + "_stack.png";
            c2.SaveAs(&cr_save_name[0]);
        //cut plot
            TCanvas c3("c3","",IMAGEX,IMAGEY);
            h_inc_cut->SetMarkerColor(kBlack);
            h_inc_cut->SetLineColor(kBlack);
            h_inc_cut->SetFillColor(kBlack);
            h_inc_cut->SetFillStyle(3017);
            //h_inc_cut->SetAxicutange(0,YMAX,"Y");
            cut_stack->Draw("hist");
            h_inc_cut->Draw("E2,same");
            auto cut_save_name = "plots/stack/cut_" + (string)dist + "_stack.png";
            c3.SaveAs(&cut_save_name[0]);
        // output
        double xs, errerr;
        xs = h_inc_cut->IntegralAndError(1,1000, errerr);
        cout<<"*************************"<<endl;
        cout<<"XS = "<<xs<<"\t"<<errerr<<endl;
        cout<<"*************************"<<endl;

}

int main()
{
    //signal
        cout<<"m4l:"<<endl;
        plotPropSig("output/stack_out/sigs/m4l.root", "llll_m", 0);
        cout<<"mjj:"<<endl;
        plotPropSig("output/stack_out/sigs/mjj.root", "jj_m", 0);
        cout<<"jj_delta_phi"<<endl;
        plotPropSig("output/stack_out/sigs/jjDelPhi.root", "jj_dphi", 0);

        plotPropSig("output/stack_out/sigs/ats.root", "Angle_theta_star", 0);
        plotPropSig("output/stack_out/sigs/aps.root", "Angle_phi_star", 0);
        plotPropSig("output/stack_out/sigs/ap0.root", "Angle_phi0", 0);
        plotPropSig("output/stack_out/sigs/ap1.root", "Angle_phi1", 0);
        plotPropSig("output/stack_out/sigs/at1.root", "Angle_theta1", 0);
        plotPropSig("output/stack_out/sigs/at2.root", "Angle_theta2", 0);


        cout<<"llll_truthBorn_m"<<endl;
        plotPropSig("output/stack_out/sigs/llll_true_m.root", "llll_truthBorn_m", 1);
            cout<<"jj_truthBorn_m"<<endl;
        plotPropSig("output/stack_out/sigs/jj_true_m.root", "jj_truthBorn_m" , 1);
            cout<<"llll_truthBorn_m"<<endl;
        plotPropSig("output/stack_out/sigs/jj_true_delphi.root", "jj_truthBorn_dphi" ,1);

        plotPropSig("output/stack_out/sigs/ats_true.root", "Angle_truthBorn_theta_star", 1);
        plotPropSig("output/stack_out/sigs/aps_true.root", "Angle_truthBorn_phi_star", 1);
        plotPropSig("output/stack_out/sigs/ap0_true.root", "Angle_truthBorn_phi0", 1);
        plotPropSig("output/stack_out/sigs/ap1_true.root", "Angle_truthBorn_phi1", 1);
        plotPropSig("output/stack_out/sigs/at1_true.root", "Angle_truthBorn_theta1", 1);
        plotPropSig("output/stack_out/sigs/at2_true.root", "Angle_truthBorn_theta2", 1);
    //bkgs
        cout<<"m4l:"<<endl;
        plotPropBkg("output/stack_out/bkgs/m4l.root", "llll_m", 0);
        cout<<"mjj:"<<endl;
        plotPropBkg("output/stack_out/bkgs/mjj.root", "jj_m", 0);
        cout<<"jj_delta_phi"<<endl;
        plotPropBkg("output/stack_out/bkgs/jjDelPhi.root", "jj_dphi", 0);

        plotPropBkg("output/stack_out/bkgs/ats.root", "Angle_theta_star", 0);
        plotPropBkg("output/stack_out/bkgs/aps.root", "Angle_phi_star", 0);
        plotPropBkg("output/stack_out/bkgs/ap0.root", "Angle_phi0", 0);
        plotPropBkg("output/stack_out/bkgs/ap1.root", "Angle_phi1", 0);
        plotPropBkg("output/stack_out/bkgs/at1.root", "Angle_theta1", 0);
        plotPropBkg("output/stack_out/bkgs/at2.root", "Angle_theta2", 0);


        cout<<"llll_truthBorn_m"<<endl;
        plotPropBkg("output/stack_out/bkgs/llll_true_m.root", "llll_truthBorn_m", 1);
            cout<<"jj_truthBorn_m"<<endl;
        plotPropBkg("output/stack_out/bkgs/jj_true_m.root", "jj_truthBorn_m" , 1);
            cout<<"llll_truthBorn_m"<<endl;
        plotPropBkg("output/stack_out/bkgs/jj_true_delphi.root", "jj_truthBorn_dphi" ,1);

        plotPropBkg("output/stack_out/bkgs/ats_true.root", "Angle_truthBorn_theta_star", 1);
        plotPropBkg("output/stack_out/bkgs/aps_true.root", "Angle_truthBorn_phi_star", 1);
        plotPropBkg("output/stack_out/bkgs/ap0_true.root", "Angle_truthBorn_phi0", 1);
        plotPropBkg("output/stack_out/bkgs/ap1_true.root", "Angle_truthBorn_phi1", 1);
        plotPropBkg("output/stack_out/bkgs/at1_true.root", "Angle_truthBorn_theta1", 1);
        plotPropBkg("output/stack_out/bkgs/at2_true.root", "Angle_truthBorn_theta2", 1);
}
