//#define debug
#include <TCanvas.h>
#include <ROOT/RDataFrame.hxx>
#include <cmath>
#define PRWLUMI (36.21 + 58.45 + 44.31) //fb-1

#define MJJ_X_MIN 0
#define MJJ_X_MAX 2000e3
#define MJJ_N_BIN 1000

#define M4L_X_MIN 0
#define M4L_X_MAX 1000e3
#define M4L_N_BIN 1000

#define PHI_X_MIN -M_PI
#define PHI_X_MAX  M_PI
#define PHI_N_BIN  10000

TH1D getHisto(ROOT::RDataFrame df, const char* filterr, const char* dist, const char* weight, int nbins, double xmin, double xmax)
{
    auto histoname = "h_" + (std::string)dist + "_" + (std::string)filterr;
    auto filteddf = df.Filter(filterr); 
    auto h = filteddf.Histo1D({&histoname[0],"",nbins,xmin,xmax},dist,weight);
    return (*h);
}

void draw(const char* in_file, const char* out_histo)
{
    ROOT::EnableImplicitMT();
    //Read files
        ROOT::RDataFrame df("SM4L_Nominal", in_file);

    //book histos
        TH1::SetDefaultSumw2();

        TH1D h_mjj_cut = getHisto(df, "pass_cut",       "jj_m", "weight", MJJ_N_BIN, MJJ_X_MIN, MJJ_X_MAX);
        TH1D h_mjj_sr  = getHisto(df, "pass_SR",        "jj_m", "weight", MJJ_N_BIN, MJJ_X_MIN, MJJ_X_MAX);
        TH1D h_mjj_cr = getHisto(df, "pass_CR",         "jj_m", "weight", MJJ_N_BIN, MJJ_X_MIN, MJJ_X_MAX);

        TH1D h_m4l_cut = getHisto(df, "pass_cut",       "llll_m", "weight", M4L_N_BIN, M4L_X_MIN, M4L_X_MAX);
        TH1D h_m4l_sr  = getHisto(df, "pass_SR",        "llll_m", "weight", M4L_N_BIN, M4L_X_MIN, M4L_X_MAX);
        TH1D h_m4l_cr = getHisto(df, "pass_CR",         "llll_m", "weight", M4L_N_BIN, M4L_X_MIN, M4L_X_MAX);

        TH1D h_phi_cut = getHisto(df, "pass_cut",       "jj_dphi", "weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
        TH1D h_phi_sr  = getHisto(df, "pass_SR",        "jj_dphi", "weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
        TH1D h_phi_cr = getHisto(df, "pass_CR",         "jj_dphi", "weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);

        TH1D h_ats_cut = getHisto(df, "pass_cut",       "Angle_theta_star", "weight",  PHI_N_BIN, 0, PHI_X_MAX);
        TH1D h_ats_sr  = getHisto(df, "pass_SR",        "Angle_theta_star", "weight",  PHI_N_BIN, 0, PHI_X_MAX);
        TH1D h_ats_cr  = getHisto(df, "pass_CR",         "Angle_theta_star", "weight", PHI_N_BIN, 0, PHI_X_MAX);

        TH1D h_aps_cut = getHisto(df, "pass_cut",       "Angle_phi_star", "weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX); //aps for angle phi star
        TH1D h_aps_sr  = getHisto(df, "pass_SR",        "Angle_phi_star", "weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
        TH1D h_aps_cr  = getHisto(df, "pass_CR",         "Angle_phi_star", "weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);

        TH1D h_p0_cut = getHisto(df, "pass_cut",       "Angle_phi0", "weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
        TH1D h_p0_sr  = getHisto(df, "pass_SR",        "Angle_phi0", "weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
        TH1D h_p0_cr  = getHisto(df, "pass_CR",         "Angle_phi0", "weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);

        TH1D h_p1_cut = getHisto(df, "pass_cut",       "Angle_phi1", "weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
        TH1D h_p1_sr  = getHisto(df, "pass_SR",        "Angle_phi1", "weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
        TH1D h_p1_cr  = getHisto(df, "pass_CR",         "Angle_phi1", "weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);

        TH1D h_t1_cut = getHisto(df, "pass_cut",       "Angle_theta1", "weight", PHI_N_BIN, 0, PHI_X_MAX);
        TH1D h_t1_sr  = getHisto(df, "pass_SR",        "Angle_theta1", "weight",  PHI_N_BIN, 0, PHI_X_MAX);
        TH1D h_t1_cr  = getHisto(df, "pass_CR",         "Angle_theta1", "weight",PHI_N_BIN, 0, PHI_X_MAX);

        TH1D h_t2_cut = getHisto(df, "pass_cut",       "Angle_theta2", "weight",  PHI_N_BIN, 0, PHI_X_MAX);
        TH1D h_t2_sr  = getHisto(df, "pass_SR",        "Angle_theta2", "weight", PHI_N_BIN, 0, PHI_X_MAX);
        TH1D h_t2_cr  = getHisto(df, "pass_CR",         "Angle_theta2", "weight", PHI_N_BIN, 0, PHI_X_MAX);


        //scale
        h_mjj_cut.Scale(1.0 / PRWLUMI);
        h_mjj_sr .Scale(1.0 / PRWLUMI);
        h_mjj_cr.Scale(1.0 / PRWLUMI);

        h_m4l_cut.Scale(1.0 / PRWLUMI);
        h_m4l_sr .Scale(1.0 / PRWLUMI);
        h_m4l_cr.Scale(1.0 / PRWLUMI);

        h_phi_cut.Scale(1.0 / PRWLUMI);
        h_phi_sr .Scale(1.0 / PRWLUMI);
        h_phi_cr .Scale(1.0 / PRWLUMI);

        h_ats_cut.Scale(1.0 / PRWLUMI);
        h_ats_sr .Scale(1.0 / PRWLUMI);
        h_ats_cr .Scale(1.0 / PRWLUMI);
        h_aps_cut.Scale(1.0 / PRWLUMI);
        h_aps_sr .Scale(1.0 / PRWLUMI);
        h_aps_cr .Scale(1.0 / PRWLUMI);
        h_p0_cut .Scale(1.0 / PRWLUMI);
        h_p0_sr  .Scale(1.0 / PRWLUMI);
        h_p0_cr  .Scale(1.0 / PRWLUMI);
        h_p1_cut .Scale(1.0 / PRWLUMI);
        h_p1_sr  .Scale(1.0 / PRWLUMI);
        h_p1_cr  .Scale(1.0 / PRWLUMI);
        h_t1_cut .Scale(1.0 / PRWLUMI);
        h_t1_sr  .Scale(1.0 / PRWLUMI);
        h_t1_cr  .Scale(1.0 / PRWLUMI);
        h_t2_cut .Scale(1.0 / PRWLUMI);
        h_t2_sr  .Scale(1.0 / PRWLUMI);
        h_t2_cr  .Scale(1.0 / PRWLUMI);



        
    //save
        TFile* out = TFile::Open(out_histo,"recreate");

        h_mjj_cut.Write();
        h_mjj_sr.Write();
        h_mjj_cr.Write();


        h_m4l_cut.Write();
        h_m4l_sr.Write();
        h_m4l_cr.Write();


        h_phi_cut.Write();
        h_phi_sr.Write();
        h_phi_cr.Write();

        h_ats_cut.Write();
        h_ats_sr .Write();
        h_ats_cr .Write();
        h_aps_cut.Write();
        h_aps_sr .Write();
        h_aps_cr .Write();
        h_p0_cut .Write();
        h_p0_sr  .Write();
        h_p0_cr  .Write();
        h_p1_cut .Write();
        h_p1_sr  .Write();
        h_p1_cr  .Write();
        h_t1_cut .Write();
        h_t1_sr  .Write();
        h_t1_cr  .Write();
        h_t2_cut .Write();
        h_t2_sr  .Write();
        h_t2_cr  .Write();


        if (df.HasColumn("fid_weight")){
            
            TH1D h_true_mjj_cut = getHisto(df, "pass_cut_truthBorn",       "jj_truthBorn_m", "fid_weight", MJJ_N_BIN, MJJ_X_MIN, MJJ_X_MAX);
            TH1D h_true_mjj_sr  = getHisto(df, "pass_SR_truthBorn",        "jj_truthBorn_m", "fid_weight", MJJ_N_BIN, MJJ_X_MIN, MJJ_X_MAX);
            TH1D h_true_mjj_cr  = getHisto(df, "pass_CR_truthBorn",        "jj_truthBorn_m", "fid_weight", MJJ_N_BIN, MJJ_X_MIN, MJJ_X_MAX);

            TH1D h_true_m4l_cut = getHisto(df, "pass_cut_truthBorn",       "llll_truthBorn_m", "fid_weight", M4L_N_BIN, M4L_X_MIN, M4L_X_MAX);
            TH1D h_true_m4l_sr  = getHisto(df, "pass_SR_truthBorn",        "llll_truthBorn_m", "fid_weight", M4L_N_BIN, M4L_X_MIN, M4L_X_MAX);
            TH1D h_true_m4l_cr  = getHisto(df, "pass_CR_truthBorn",        "llll_truthBorn_m", "fid_weight", M4L_N_BIN, M4L_X_MIN, M4L_X_MAX);

            TH1D h_true_phi_cut = getHisto(df, "pass_cut_truthBorn",       "jj_truthBorn_dphi", "fid_weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
            TH1D h_true_phi_sr  = getHisto(df, "pass_SR_truthBorn",        "jj_truthBorn_dphi", "fid_weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
            TH1D h_true_phi_cr  = getHisto(df, "pass_CR_truthBorn",        "jj_truthBorn_dphi", "fid_weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);

            TH1D h_true_ats_cut = getHisto(df, "pass_cut_truthBorn",       "Angle_truthBorn_theta_star", "fid_weight", PHI_N_BIN, 0, PHI_X_MAX);
            TH1D h_true_ats_sr  = getHisto(df, "pass_SR_truthBorn",        "Angle_truthBorn_theta_star", "fid_weight", PHI_N_BIN, 0, PHI_X_MAX);
            TH1D h_true_ats_cr  = getHisto(df, "pass_CR_truthBorn",        "Angle_truthBorn_theta_star", "fid_weight", PHI_N_BIN, 0, PHI_X_MAX);

            TH1D h_true_aps_cut = getHisto(df, "pass_cut_truthBorn",       "Angle_truthBorn_phi_star", "fid_weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX); //aps for angle phi star
            TH1D h_true_aps_sr  = getHisto(df, "pass_SR_truthBorn",        "Angle_truthBorn_phi_star", "fid_weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
            TH1D h_true_aps_cr  = getHisto(df, "pass_CR_truthBorn",        "Angle_truthBorn_phi_star", "fid_weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);

            TH1D h_true_p0_cut = getHisto(df, "pass_cut_truthBorn",       "Angle_truthBorn_phi0", "fid_weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
            TH1D h_true_p0_sr  = getHisto(df, "pass_SR_truthBorn",        "Angle_truthBorn_phi0", "fid_weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
            TH1D h_true_p0_cr  = getHisto(df, "pass_CR_truthBorn",        "Angle_truthBorn_phi0", "fid_weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);

            TH1D h_true_p1_cut = getHisto(df, "pass_cut_truthBorn",       "Angle_truthBorn_phi1", "fid_weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
            TH1D h_true_p1_sr  = getHisto(df, "pass_SR_truthBorn",        "Angle_truthBorn_phi1", "fid_weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);
            TH1D h_true_p1_cr  = getHisto(df, "pass_CR_truthBorn",        "Angle_truthBorn_phi1", "fid_weight", PHI_N_BIN, PHI_X_MIN, PHI_X_MAX);

            TH1D h_true_t1_cut = getHisto(df, "pass_cut_truthBorn",       "Angle_truthBorn_theta1", "fid_weight", PHI_N_BIN, 0, PHI_X_MAX);
            TH1D h_true_t1_sr  = getHisto(df, "pass_SR_truthBorn",        "Angle_truthBorn_theta1", "fid_weight", PHI_N_BIN, 0, PHI_X_MAX);
            TH1D h_true_t1_cr  = getHisto(df, "pass_CR_truthBorn",        "Angle_truthBorn_theta1", "fid_weight", PHI_N_BIN, 0, PHI_X_MAX);

            TH1D h_true_t2_cut = getHisto(df, "pass_cut_truthBorn",       "Angle_truthBorn_theta2", "fid_weight", PHI_N_BIN, 0, PHI_X_MAX);
            TH1D h_true_t2_sr  = getHisto(df, "pass_SR_truthBorn",        "Angle_truthBorn_theta2", "fid_weight", PHI_N_BIN, 0, PHI_X_MAX);
            TH1D h_true_t2_cr  = getHisto(df, "pass_CR_truthBorn",        "Angle_truthBorn_theta2", "fid_weight", PHI_N_BIN, 0, PHI_X_MAX);

            h_true_mjj_cut.Scale(1.0 / PRWLUMI);
            h_true_mjj_sr .Scale(1.0 / PRWLUMI);
            h_true_mjj_cr.Scale(1.0 / PRWLUMI);

            h_true_m4l_cut.Scale(1.0 / PRWLUMI);
            h_true_m4l_sr .Scale(1.0 / PRWLUMI);
            h_true_m4l_cr.Scale(1.0 / PRWLUMI);

            h_true_phi_cut.Scale(1.0 / PRWLUMI);
            h_true_phi_sr .Scale(1.0 / PRWLUMI);
            h_true_phi_cr.Scale(1.0 / PRWLUMI);

            h_true_ats_cut.Scale(1.0 / PRWLUMI);
            h_true_ats_sr .Scale(1.0 / PRWLUMI);
            h_true_ats_cr .Scale(1.0 / PRWLUMI);
            h_true_aps_cut.Scale(1.0 / PRWLUMI);
            h_true_aps_sr .Scale(1.0 / PRWLUMI);
            h_true_aps_cr .Scale(1.0 / PRWLUMI);
            h_true_p0_cut .Scale(1.0 / PRWLUMI);
            h_true_p0_sr  .Scale(1.0 / PRWLUMI);
            h_true_p0_cr  .Scale(1.0 / PRWLUMI);
            h_true_p1_cut .Scale(1.0 / PRWLUMI);
            h_true_p1_sr  .Scale(1.0 / PRWLUMI);
            h_true_p1_cr  .Scale(1.0 / PRWLUMI);
            h_true_t1_cut .Scale(1.0 / PRWLUMI);
            h_true_t1_sr  .Scale(1.0 / PRWLUMI);
            h_true_t1_cr  .Scale(1.0 / PRWLUMI);
            h_true_t2_cut .Scale(1.0 / PRWLUMI);
            h_true_t2_sr  .Scale(1.0 / PRWLUMI);
            h_true_t2_cr  .Scale(1.0 / PRWLUMI);


            h_true_mjj_cut.Write();
            h_true_mjj_sr.Write();
            h_true_mjj_cr.Write();


            h_true_m4l_cut.Write();
            h_true_m4l_sr.Write();
            h_true_m4l_cr.Write();


            h_true_phi_cut.Write();
            h_true_phi_sr.Write();
            h_true_phi_cr.Write();

            h_true_ats_cut.Write();
            h_true_ats_sr .Write();
            h_true_ats_cr .Write();
            h_true_aps_cut.Write();
            h_true_aps_sr .Write();
            h_true_aps_cr .Write();
            h_true_p0_cut .Write();
            h_true_p0_sr  .Write();
            h_true_p0_cr  .Write();
            h_true_p1_cut .Write();
            h_true_p1_sr  .Write();
            h_true_p1_cr  .Write();
            h_true_t1_cut .Write();
            h_true_t1_sr  .Write();
            h_true_t1_cr  .Write();
            h_true_t2_cut .Write();
            h_true_t2_sr  .Write();
            h_true_t2_cr  .Write();

        }

        out->Close();
}
#ifndef debug
int main(int argc, char** argv){
    const char* in_file = argv[1];
    const char* out_histo = argv[2];

    draw(in_file,out_histo);
}
#else
int main(){
    const char* in_file = "output/analyse_out/999_all/364250.Sherpa_222_NNPDF30NNLO_llll.root";
    const char* out_histo = "output/histo_out/999_all/364250.Sherpa_222_NNPDF30NNLO_llll.root";
    draw(in_file,out_histo);
}
#endif