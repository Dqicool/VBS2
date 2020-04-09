//claculation of interesting observables
#include "libs/genAna.h"

#define MT
//#define debug
//#define PRESERVEALL
using namespace qidong;

void analyse(char* in_file1, char* in_file2, char* in_file3, char* out_anaed_tree)
{
    #ifdef MT
    ROOT::EnableImplicitMT();
    #endif
    ROOT::RDataFrame dframe("SM4L_Nominal", {in_file1, in_file2, in_file3});
    //analyse
        auto ana=dframe.Define("dete",      getDete, {}).
                        Define("event",     getEvt, {"jet_pt", "jet_eta", "jet_phi", "jet_m", "lepton_pt", "lepton_eta", "lepton_phi", "lepton_m", "lepton_particleID", "lepton_charge","weight","dete" ,"lepton_d0raw", "lepton_d0sig", "lepton_z0sinTheta", "lepton_passesIsolation","jet_Jvt"}).
                        Define("jj_m",      findMJJ, {"event"}).
                        Define("jj_dphi",   findDPhiJJ, {"event"}).
                        Define("jj_dy",     findDYJJ, {"event"}).
                        Define("llll_m",    findM4l, {"event"}).
                        Define("z1_m",      findMZ1, {"event"}).
                        Define("z2_m",      findMZ2, {"event"}).
                        Define("pt_bala",   findPtBa, {"event"}).
                        Define("centr",     findCent, {"event"}).
                        Define("pass_cut",  passCut, {"event"}).
                        Define("pass_SR",   passSR, {"event"}).
                        Define("pass_CR",   passCR, {"event"}).
                        Define("passZZ", pass_z1z2, {"event"}).
                        Define("passSin", pass_sin, {"event"}).
                        Define("passJJ", pass_j1j2, {"event"}).
                        Define("Angle_theta_star", findAngleThetaStar, {"event"}).
                        Define("Angle_phi_star", findAnglePhiStar, {"event"}).
                        Define("Angle_phi0", findAnglePhi0, {"event"}).
                        Define("Angle_phi1", findAnglePhi1, {"event"}).
                        Define("Angle_theta1", findAngleTheta1, {"event"}).
                        Define("Angle_theta2", findAngleTheta2, {"event"});


        if(dframe.HasColumn("fid_weight")){
            ana =   ana.Define("truth", getTrue, {}).
                        Define("event_truthBorn",   getTrueEvt, {"jet_truthBorn_pt", "jet_truthBorn_eta", "jet_truthBorn_phi", "jet_truthBorn_m","lepton_truthBorn_pt", "lepton_truthBorn_eta", "lepton_truthBorn_phi", "lepton_truthBorn_m", "lepton_truthBorn_particleID", "lepton_truthBorn_charge", "fid_weight", "truth"}).
                        Define("jj_truthBorn_m",    findMJJ, {"event_truthBorn"}).
                        Define("jj_truthBorn_dphi", findDPhiJJ, {"event_truthBorn"}).
                        Define("jj_truthBorn_dy",   findDYJJ, {"event_truthBorn"}).
                        Define("llll_truthBorn_m",  findM4l, {"event_truthBorn"}).
                        Define("z1_truthBorn_m",    findMZ1, {"event_truthBorn"}).
                        Define("z2_truthBorn_m",    findMZ2, {"event_truthBorn"}).
                        Define("pt_truthBorn_bala", findPtBa, {"event_truthBorn"}).
                        Define("centr_truthBorn",   findCent, {"event_truthBorn"}).
                        Define("pass_cut_truthBorn",passCut, {"event_truthBorn"}).
                        Define("pass_SR_truthBorn", passSR, {"event_truthBorn"}).
                        Define("pass_CR_truthBorn", passCR, {"event_truthBorn"}).
                        Define("Angle_truthBorn_theta_star", findAngleThetaStar, {"event_truthBorn"}).
                        Define("Angle_truthBorn_phi_star", findAnglePhiStar, {"event_truthBorn"}).
                        Define("Angle_truthBorn_phi0", findAnglePhi0, {"event_truthBorn"}).
                        Define("Angle_truthBorn_phi1", findAnglePhi1, {"event_truthBorn"}).
                        Define("Angle_truthBorn_theta1", findAngleTheta1, {"event_truthBorn"}).
                        Define("Angle_truthBorn_theta2", findAngleTheta2, {"event_truthBorn"});
        }
    //save tree
        auto hehe = ana.GetColumnNames();
        if (ana.HasColumn("fid_weight")){
            hehe = { "weight", "jj_m", "jj_dphi", "jj_dy", "llll_m", "z1_m", "z2_m", "pt_bala", "centr", "pass_cut", "pass_SR", "pass_CR", 
                     "Angle_theta_star", "Angle_phi_star", "Angle_phi0", "Angle_phi1", "Angle_theta1","Angle_theta2",
                     "fid_weight", "jj_truthBorn_m", "jj_truthBorn_dphi", "jj_truthBorn_dy",  "llll_truthBorn_m", "z1_truthBorn_m",   "z2_truthBorn_m",   "pt_truthBorn_bala","centr_truthBorn",  "pass_cut_truthBorn", "pass_SR_truthBorn", "pass_CR_truthBorn",
                     "Angle_truthBorn_theta_star", "Angle_truthBorn_phi_star", "Angle_truthBorn_phi0", "Angle_truthBorn_phi1", "Angle_truthBorn_theta1", "Angle_truthBorn_theta2"};
            #ifndef PRESERVEALL
            ana.Filter("pass_cut == 1 || pass_cut_truthBorn == 1").Snapshot("SM4L_Nominal", out_anaed_tree, hehe);
            #else
            ana.Snapshot("SM4L_Nominal", out_anaed_tree, hehe);
            #endif
        }
        else{
            hehe = { "weight", "jj_m", "jj_dphi", "jj_dy", "llll_m", "z1_m", "z2_m", "pt_bala", "centr", "pass_cut", "pass_SR", "pass_CR", "Angle_theta_star", "Angle_phi_star", "Angle_phi0", "Angle_phi1", "Angle_theta1","Angle_theta2",};
            #ifndef PRESERVEALL
            ana.Filter("pass_cut == 1").Snapshot("SM4L_Nominal", out_anaed_tree, hehe);
            #else
            ana.Snapshot("SM4L_Nominal", out_anaed_tree, hehe);
            #endif
        }
}

#ifndef debug
int main(int argc, char** argv)
{
    char* in_file1 = argv[1];
    char* in_file2 = argv[2];
    char* in_file3 = argv[3];
    char* out_tree = argv[4];

    analyse(in_file1, in_file2, in_file3, out_tree);
}
#else
int main(){
    char* in_file1 = "data/mc16_13TeV.364364.Sherpa_222_NNPDF30NNLO_lllljj_EW6_noHiggs.deriv.DAOD_HIGG2D1.e6611_e5984_a875_r10201_r10210_p3654.root";
    char* in_file2 = "data/mc16_13TeV.364364.Sherpa_222_NNPDF30NNLO_lllljj_EW6_noHiggs.deriv.DAOD_HIGG2D1.e6611_e5984_a875_r10724_r10726_p3654.root";
    char* in_file3 = "data/mc16_13TeV.364364.Sherpa_222_NNPDF30NNLO_lllljj_EW6_noHiggs.deriv.DAOD_HIGG2D1.e6611_e5984_a875_r9364_r9315_p3654.root";
    char* out_tree = "debug.root";
    analyse(in_file1, in_file2, in_file3, out_tree);
}
#endif