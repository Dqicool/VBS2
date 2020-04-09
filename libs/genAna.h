#include <TROOT.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include "vector"
#include <TMath.h>
#include <iostream>
#include <string.h>
#include <ROOT/RDataFrame.hxx>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <Math/GenVector/PtEtaPhiM4Dfwd.h>
#include <Math/Vector4Dfwd.h>
#include<Math/Vector3Dfwd.h>
#include<Math/GenVector/Boost.h>
#include <TChain.h>
#include<TROOT.h>
#include<TStyle.h>
#include<TLatex.h>

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/EulerAngles.h"
#include "Math/GenVector/AxisAngle.h"
#include "Math/GenVector/Quaternion.h"
#include "Math/GenVector/LorentzRotation.h"
#include "Math/GenVector/Boost.h"
#include "Math/GenVector/Transform3D.h"
#include "Math/GenVector/Plane3D.h"
#include "Math/GenVector/VectorUtil.h"

#define Z_MASS 91.1876e3
#define GeV 1e3
namespace qidong{
    using namespace ROOT::Math;
    typedef std::vector<PtEtaPhiMVector> FourMomenta;

    class evt{
        public:
            class CP_angles{
                public:
                    float theta_star;
                    float phi_star;
                    float phi;
                    float phi1;
                    float theta1;
                    float theta2;
                    void veto(){
                        theta_star = -999;
                        phi_star = -999;
                        phi = -999;
                        phi1 = -999;
                        theta1 = -999;
                        theta2 = -999;
                    }
                    CP_angles(){
                        veto();
                    }
            };
            class results{
                public:
                    float m_jj;
                    float dphi_jj;
                    float dy_jj;
                    float m_4l;
                    float m_z1;
                    float m_z2;
                    float pt_bala;
                    float centra;
                    CP_angles angles;
                    
                    void veto(){
                        m_jj = -999;
                        m_jj = -999;
                        dphi_jj = -999;
                        dy_jj = -999;
                        m_4l = -999;
                        m_z1 = -999;
                        m_z2 = -999;
                        pt_bala = -999;
                        centra = -999;
                        angles.veto();
                    }
            };
            class Djet{
                public:
                    PtEtaPhiMVector mom;
                    float jvt;
                    bool is_truth;
                    Djet(float _pt, float _eta, float _phi, float _m, bool _is_truth, float _jvt=-999){
                        mom = PtEtaPhiMVector(_pt, _eta, _phi, _m);
                        if(_is_truth)   jvt = -999;
                        else            jvt = _jvt;
                    }
                    float abseta(){
                        return TMath::Abs(mom.eta());
                    }
                    float pt(){
                        return mom.pt();
                    }
                    float rapidity(){
                        return mom.Rapidity();
                    }
                    float E(){
                        return mom.energy();
                    }
            };
            class Dlepton{
                public:
                    PtEtaPhiMVector mom;
                    int id;
                    int charge;
                    bool is_truth;
                    char isolation;
                    float d0raw;
                    float d0sig;
                    float z0sintheta;
                    Dlepton( float  _pt,      float _eta,          float _phi,     float _m, 
                            int    _id,       float    _charge,       bool _is_truth,   
                            float _d0raw = -999,   
                            float _d0sig = -999,    
                            float _z0sintheta = -999,   
                            char   _isolation = '-'){
                        mom = PtEtaPhiMVector(_pt, _eta, _phi, _m);
                        id = _id;
                        charge = (int)_charge;
                        is_truth = _is_truth;
                        if(is_truth){
                            d0raw = -999;
                            d0sig = -999;
                            z0sintheta = -999;
                            isolation = '-';
                        }
                        else{
                            d0raw = _d0raw;
                            d0sig = _d0sig;
                            z0sintheta = _z0sintheta;
                            isolation = _isolation;
                        }
                    }
                    float abseta(){
                        return TMath::Abs(mom.eta());
                    }
                    float pt(){
                        return mom.pt();
                    }
                    int abspid(){
                        return id;
                    }
                    int pid(){
                        return id*charge;
                    }
            };

            typedef std::vector<Djet> Jets;
            typedef std::vector<Dlepton> Leptons;
            
            Jets                jets_raw;
            Leptons             leptons_raw;

            Jets                jets_pass;
            Leptons             leptons_pass;

            std::vector<Leptons> lepton_pair;
            FourMomenta         j1j2;
            FourMomenta         z1z2;

            bool                is_truth;
            float              weight;
            results             res;
            bool                pass_sin_par = 0;
            bool                pass_j1j2    = 0;
            bool                pass_z1z2    = 0;
            bool                pass_cut     = 0;
            bool                pass_SR      = 0;
            bool                pass_CR      = 0;
            bool                veto         = 0;

            Jets jetSel(Jets jets){
                auto siz = jets.size(); 
                Jets ret;
                for (uint i = 0; i < siz; i++){
                    if (jets[i].abseta() < 2.4 && jets[i].pt() > 30*GeV && jets[i].jvt > 0.6){
                        ret.push_back(jets[i]);
                    }
                    else if(jets[i].abseta() < 4.5 && jets[i].pt() > 40*GeV){
                        ret.push_back(jets[i]);
                    }
                } 
                return ret;
            }
            Leptons lepSel(Leptons leptons){
                auto siz = leptons.size();
                Leptons ret;
                for (uint i = 0;i < siz; i++)
                {
                    if (leptons[i].isolation == 1)
                    {
                        if(leptons[i].abspid() == 11){
                            if(leptons[i].pt() > 7*GeV && leptons[i].abseta() < 2.47 && TMath::Abs(leptons[i].d0sig) < 5 && TMath::Abs(leptons[i].z0sintheta) < 0.5){
                                ret.push_back(leptons[i]);
                            }
                        }
                        else if(leptons[i].abspid() == 13)
                        {
                            if(leptons[i].pt() > 7*GeV && leptons[i].abseta() < 2.7 && TMath::Abs(leptons[i].d0sig) < 3 && TMath::Abs(leptons[i].z0sintheta) < 0.5){
                                ret.push_back(leptons[i]);
                            }
                        }
                    }
                }
                return ret;
            }
            Jets jetTruSel(Jets jets){
                auto siz = jets.size(); 
                Jets ret;
                for (uint i = 0; i < siz; i++){
                    if (jets[i].abseta() < 2.4 && jets[i].pt() > 30*GeV){
                        ret.push_back(jets[i]);
                    }
                    else if(jets[i].abseta() < 4.5 && jets[i].pt() > 40*GeV){
                        ret.push_back(jets[i]);
                    }
                } 
                return ret;
            }
            Leptons lepTruSel(Leptons leptons){
                auto siz = leptons.size();
                Leptons ret;
                for (uint i = 0;i < siz; i++)
                {
                    //if (lepton_isTight[i] == 0)
                    {
                        if(leptons[i].abspid() == 11){
                            if(leptons[i].pt() > 7*GeV && leptons[i].abseta() < 2.47){
                                ret.push_back(leptons[i]);
                            }
                        }
                        else if(leptons[i].abspid() == 13)
                        {
                            if(leptons[i].pt() > 7*GeV && leptons[i].abseta() < 2.7){
                                ret.push_back(leptons[i]);
                            }
                        }
                    }
                }
                return ret;
            }
            Jets getJets(std::vector<float> _pts, std::vector<float> _etas, std::vector<float> _phis, std::vector<float> _ms, bool _is_truth, std::vector<float> _jvts = {}){
                Jets ret;
                for(uint i = 0;i < _pts.size(); i++){
                    ret.push_back(Djet(_pts[i], _etas[i], _phis[i], _ms[i], _is_truth, _jvts[i]));
                }
                return ret;
            }
            Jets getTrueJets(std::vector<float> _pts, std::vector<float> _etas, std::vector<float> _phis, std::vector<float> _ms, bool _is_truth){
                Jets ret;
                for(uint i = 0;i < _pts.size(); i++){
                ret.push_back(Djet(_pts[i], _etas[i], _phis[i], _ms[i], _is_truth));
                }
                return ret;
            }
            Leptons getLeps(std::vector<float> _pts,      std::vector<float> _etas,          std::vector<float> _phis,     std::vector<float> _ms, 
                                        std::vector<int>    _ids,       std::vector<float>    _charges,       bool _is_truth,   
                                        std::vector<float> _d0raws,   
                                        std::vector<float> _d0sigs,    
                                        std::vector<float> _z0sinthetas,   
                                        std::vector<char>   _isolations){
                Leptons ret;
                for(uint i = 0;i < _pts.size(); i++){
                    ret.push_back(Dlepton(_pts[i], _etas[i], _phis[i], _ms[i], 
                                        _ids[i], _charges[i], _is_truth, _d0raws[i], _d0sigs[i], _z0sinthetas[i], _isolations[i]));
                }
                return ret;
            }
            Leptons getTrueLeps(std::vector<float> _pts,      std::vector<float> _etas,          std::vector<float> _phis,     std::vector<float> _ms, 
                                        std::vector<int>    _ids,       std::vector<float>    _charges,       bool _is_truth){
                Leptons ret;
                for(uint i = 0;i < _pts.size(); i++){
                    ret.push_back(Dlepton(_pts[i], _etas[i], _phis[i], _ms[i], 
                                        _ids[i], _charges[i], _is_truth));
                }
                return ret;
            }
            static bool jetsPairECriterion(Jets a, Jets b){
                float Ea = 0, Eb = 0;
                for(Djet ai : a){
                    Ea += ai.E();
                }
                for(Djet bi : b){
                    Eb += bi.E();
                }
                return Ea > Eb;
            }
            static bool jetsPtCriterion(Djet a, Djet b){
                return a.pt() > b.pt();
            }
            FourMomenta getJ1J2(Jets jets){
                FourMomenta ret{};
                uint size = jets.size();
                if (size < 2) {}
                else{
                    std::vector<Jets> jets_pairs{};
                    for (uint i = 0; i < size; i++)
                    {
                        for (uint j = i + 1; j < size; j++)
                        {
                            if(jets[i].rapidity() * jets[j].rapidity() < 0)
                            {
                                Jets tmp_pairs{jets[i], jets[j]};
                                jets_pairs.push_back(tmp_pairs);
                            }
                        }
                    }
                    if(jets_pairs.size() > 0)
                    {
                        std::sort(jets_pairs.begin(), jets_pairs.end(), jetsPairECriterion);
                        std::sort(jets_pairs[0].begin(), jets_pairs[0].end(), jetsPtCriterion);
                        ret = {jets_pairs[0][0].mom, jets_pairs[0][1].mom};
                        if(jets_pairs[0][0].pt() < 40*GeV || jets_pairs[0][1].pt() < 30*GeV) {ret = {};}
                    }
                }
                return ret;
            }
            float getDeltaR(PtEtaPhiMVector a, PtEtaPhiMVector b){
                return std::sqrt((a.eta()-b.eta())*(a.eta()-b.eta()) + (a.phi()-b.phi())*(a.phi()-b.phi()));
            }
            std::vector<std::vector<Leptons>> getLeptonPairsComb(Leptons leps){
                std::vector<std::vector<Leptons>> ret;
                int size = leps.size();
                if(size < 4) ret = {};
                else{
                    for (int i = 0; i < size-1; i++){
                        for(int j = i+1; j < size; j++){
                            if (leps[i].pid() == -leps[j].pid()){
                                for (int k = i+1; k < size-1; k++){
                                    for (int l = k+1; l<size; l++){
                                        if (l != j && k != j){
                                            if (leps[k].pid() == -leps[l].pid()){
                                                auto delR_ij = getDeltaR(leps[i].mom, leps[j].mom);
                                                auto delR_kl = getDeltaR(leps[k].mom, leps[l].mom);
                                                if (delR_ij > 0.2 && delR_kl > 0.2){
                                                    if( leps[i].abspid() == leps[j].abspid() 
                                                        && ( leps[i].mom +  leps[j].mom).mass() > 10 
                                                        && ( leps[k].mom +  leps[l].mom).mass() > 10)
                                                    {
                                                        Leptons tmp1{leps[i], leps[j]};
                                                        Leptons tmp2{leps[k], leps[l]};
                                                        std::vector<Leptons> tmp{tmp1,tmp2};
                                                        ret.push_back(tmp);
                                                    }
                                                    else if(leps[i].abspid() != leps[k].abspid() ){
                                                        Leptons tmp1{leps[i], leps[j]};
                                                        Leptons tmp2{leps[k], leps[l]};
                                                        std::vector<Leptons> tmp{tmp1,tmp2};
                                                        ret.push_back(tmp);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                
                //std::cout<<ret.size()<<std::endl;
                return ret;
            }
            static bool lepsPairsMCriterion(Leptons a, Leptons b){
                float mass_dif_a = std::abs((a[0].mom + a[1].mom).mass() - Z_MASS);
                float mass_dif_b = std::abs((b[0].mom + b[1].mom).mass() - Z_MASS);
                return mass_dif_a < mass_dif_b;
            }
            static bool lepsPairsCombMCriterion(std::vector<Leptons> pairs_a, std::vector<Leptons> pairs_b){
                float comb_mass_diff_a =   std::abs((pairs_a[0][0].mom + pairs_a[0][1].mom).mass() - Z_MASS) +
                                            std::abs((pairs_a[1][0].mom + pairs_a[1][1].mom).mass() - Z_MASS);
                float comb_mass_diff_b =   std::abs((pairs_b[0][0].mom + pairs_b[0][1].mom).mass() - Z_MASS) +
                                            std::abs((pairs_b[1][0].mom + pairs_b[1][1].mom).mass() - Z_MASS);    
                return comb_mass_diff_a < comb_mass_diff_b;
            }

            std::vector<Leptons> getLepPair(std::vector<std::vector<Leptons>> lep_pairs_comb){
                std::vector<Leptons> ret{};
                if(lep_pairs_comb.size() >= 1){
                    for(uint i = 0; i < lep_pairs_comb.size(); i++){
                        std::sort(lep_pairs_comb[i].begin(), lep_pairs_comb[i].end(), lepsPairsMCriterion);
                    }
                    std::sort(lep_pairs_comb.begin(), lep_pairs_comb.end(), lepsPairsCombMCriterion);
                    ret = lep_pairs_comb[0];
                }
                return ret;
            }

            FourMomenta getZ1Z2(std::vector<Leptons> lep_pairs){
                FourMomenta ret{};
                if(lep_pairs.size() >= 1){
                    ret = {(lep_pairs[0][0].mom + lep_pairs[0][1].mom), 
                            (lep_pairs[1][0].mom + lep_pairs[1][1].mom)};
                }
                return ret;
            }
            float getPtBala (FourMomenta j1j2, FourMomenta z1z2){
                float ret = -999;
                if(j1j2.size() >= 2 && z1z2.size() >= 2){
                    float zzjj_sys_pt = (z1z2[0] + z1z2[1] + j1j2[0] + j1j2[1]).pt(); 
                    float zzjj_scaler_sum_pt =z1z2[0].pt()  + z1z2[1].pt() + j1j2[0].pt() + j1j2[1].pt();
                    ret = std::abs(zzjj_sys_pt / zzjj_scaler_sum_pt);
                }
                return ret;
            }
            float getDelPhiJJ(FourMomenta j1j2){
                float ret = -999;
                if(j1j2.size() >= 2){
                    if (j1j2[0].Rapidity() > j1j2[1].Rapidity())    
                        ret = j1j2[0].phi() - j1j2[1].phi();
                    else                
                        ret = j1j2[1].phi() - j1j2[0].phi();
                    
                    if(ret > M_PI){
                        ret -= 2*M_PI;
                    }
                    else if(ret < -M_PI){
                        ret += 2*M_PI;
                    }
                }
                return ret;
            }
            float getCentra (FourMomenta j1j2, FourMomenta z1z2){
                return TMath::Abs(((j1j2[0] + j1j2[1]).Rapidity() - (z1z2[0] + z1z2[1]).Rapidity()) / (j1j2[0].Rapidity() - j1j2[1].Rapidity()));     
            }

            CP_angles getAngles(std::vector<Leptons> lep_pairs, FourMomenta z1z2){
                auto starret = new CP_angles();
                auto ret = *starret;
                Boost zz_rest = Boost((z1z2[0] + z1z2[1]).BoostToCM());
                PtEtaPhiMVector q11 = zz_rest((lep_pairs[0][0].mom));
                PtEtaPhiMVector q12 = zz_rest(lep_pairs[0][1].mom);
                PtEtaPhiMVector q21 = zz_rest(lep_pairs[1][0].mom);
                PtEtaPhiMVector q22 = zz_rest(lep_pairs[1][1].mom);
                PtEtaPhiMVector q1  = zz_rest(z1z2[0]);
                PtEtaPhiMVector q2  = zz_rest(z1z2[1]);
                XYZVector n1hat  = q11.Vect().Cross(q12.Vect()) / std::sqrt((q11.Vect().Cross(q12.Vect())).mag2());
                XYZVector n2hat  = q21.Vect().Cross(q22.Vect()) / std::sqrt((q21.Vect().Cross(q22.Vect())).mag2());
                XYZVector nzhat(0.0,0.0,1.0);
                XYZVector nschat = nzhat.Cross(q1.Vect()) / std::sqrt(nzhat.Cross(q1.Vect()).Mag2());
                ret.theta_star = q1.theta();
                ret.phi_star   = q1.phi();
                ret.phi  = q1.Vect().Dot(n1hat.Cross(n2hat)) / std::abs(q1.Vect().Dot(n1hat.Cross(n2hat))) * (TMath::ACos(-n1hat.Dot(n2hat)));
                ret.phi1 = q1.Vect().Dot(n1hat.Cross(nschat)) / std::abs(q1.Vect().Dot(n1hat.Cross(nschat))) * (TMath::ACos(n1hat.Dot(nschat)));
                ret.theta1 = TMath::ACos(-(q2.Vect().Dot(q11.Vect())) / std::sqrt(q2.Vect().Mag2()) / std::sqrt(q11.Vect().Mag2()));
                ret.theta2 = TMath::ACos(-(q1.Vect().Dot(q21.Vect())) / std::sqrt(q1.Vect().Mag2()) / std::sqrt(q21.Vect().Mag2()));
                return ret;
            }
            
            evt(){
                jets_raw = {};
                leptons_raw = {};
                jets_pass = {};
                leptons_pass = {};
                lepton_pair = {};
                j1j2 = {};
                z1z2 = {};
                is_truth = 0;
                weight = -999;
                res.veto();
            }
            evt(std::vector<float> _jet_pts, std::vector<float> _jet_etas, std::vector<float> _jet_phis, std::vector<float> _jet_ms, 
                std::vector<float> _lep_pts,      std::vector<float> _lep_etas,          std::vector<float> _lep_phis,     std::vector<float> _lep_ms, 
                std::vector<int>    _lep_ids,     std::vector<float>    _lep_charges,   
                float _weight, bool _is_truth,
                std::vector<float> _lep_d0raws = {},   
                std::vector<float> _lep_d0sigs = {},    
                std::vector<float> _lep_z0sinthetas = {},   
                std::vector<char>   _lep_isolations = {},
                std::vector<float> _jet_jvts = {}){
                weight = _weight;
                is_truth = _is_truth;
                //single particle selection 
                    if (is_truth){
                        jets_raw = getTrueJets(_jet_pts, _jet_etas, _jet_phis, _jet_ms, _is_truth);
                        leptons_raw = getTrueLeps(_lep_pts, _lep_etas, _lep_phis, _lep_ms, _lep_ids, _lep_charges, _is_truth);
                        jets_pass = jetTruSel(jets_raw);
                        leptons_pass = lepTruSel(leptons_raw);
                    }
                    else{
                        jets_raw = getJets(_jet_pts, _jet_etas, _jet_phis, _jet_ms, _is_truth, _jet_jvts);
                        leptons_raw = getLeps(_lep_pts, _lep_etas, _lep_phis, _lep_ms, _lep_ids, _lep_charges, _is_truth, _lep_d0raws, _lep_d0sigs, _lep_z0sinthetas, _lep_isolations);
                        jets_pass = jetSel(jets_raw);
                        leptons_pass = lepSel(leptons_raw);
                    }
                    if(jets_pass.size() < 2 || leptons_pass.size() < 4){
                        res.veto();
                        veto = 1;
                    }
                //j1j2 selection
                    if (!veto){
                        pass_sin_par = 1;
                        j1j2 = getJ1J2(jets_pass);
                        if	(j1j2.size() < 2 || (j1j2[0] + j1j2[1]).mass() < 200*GeV || std::abs(j1j2[0].Rapidity() - j1j2[1].Rapidity()) < 2) {
                            res.veto();
                            veto = 1;
                        }  
                    }
                    
                //z1z2 sel
                    if(!veto){
                        pass_j1j2    = 1;
                        lepton_pair = getLepPair(getLeptonPairsComb(leptons_pass));
                        z1z2 = getZ1Z2(lepton_pair);
                        if	(z1z2.size() < 2 || z1z2[0].mass() < 70*GeV ||  z1z2[0].mass() > 110*GeV || z1z2[1].mass() < 70*GeV ||  z1z2[1].mass() > 110*GeV){
                            res.veto();
                            veto = 1;
                        }
                    }
                    

                //geather results
                    if(!veto){
                        pass_z1z2 = 1;
                        res.pt_bala = getPtBala(j1j2, z1z2);
                        res.m_jj    = (j1j2[0] + j1j2[1]).mass();
                        res.m_4l    = (z1z2[0] + z1z2[1]).mass();
                        res.dphi_jj = getDelPhiJJ(j1j2);
                        res.dy_jj   = TMath::Abs(j1j2[0].Rapidity() - j1j2[1].Rapidity());
                        res.m_z1    = z1z2[0].mass();
                        res.m_z2    = z1z2[1].mass();
                        res.centra  = getCentra(j1j2, z1z2);
                        if (std::abs(res.pt_bala) < 0.5){
                            pass_cut = 1;
                            if (res.centra < 0.4){ pass_SR = 1; }
                            else { pass_CR = 1; }
                        }
                        res.angles = getAngles(lepton_pair, z1z2);
                    }
                
            }
            
    };

    //For root
        bool getTrue(){ return 1;}
        bool getDete(){ return 0;}

        evt getTrueEvt(std::vector<float> _jet_pts, std::vector<float> _jet_etas, std::vector<float> _jet_phis, std::vector<float> _jet_ms, 
                std::vector<float> _lep_pts,      std::vector<float> _lep_etas,          std::vector<float> _lep_phis,     std::vector<float> _lep_ms, 
                std::vector<int>    _lep_ids,       std::vector<float>    _lep_charges,   
                float _weight, bool _is_truth){
            return evt(_jet_pts, _jet_etas, _jet_phis, _jet_ms, 
                       _lep_pts, _lep_etas, _lep_phis, _lep_ms, 
                       _lep_ids, _lep_charges,   
                       _weight,  _is_truth);
        }
        evt getEvt( std::vector<float> _jet_pts, std::vector<float> _jet_etas, std::vector<float> _jet_phis, std::vector<float> _jet_ms, 
                    std::vector<float> _lep_pts,      std::vector<float> _lep_etas,          std::vector<float> _lep_phis,     std::vector<float> _lep_ms, 
                    std::vector<int>    _lep_ids,       std::vector<float>    _lep_charges,   
                    float _weight, bool _is_truth,
                    std::vector<float> _lep_d0raws,   
                    std::vector<float> _lep_d0sigs,    
                    std::vector<float> _lep_z0sinthetas,   
                    std::vector<char>   _lep_isolations,
                    std::vector<float> _jet_jvts){
            return evt(_jet_pts, _jet_etas, _jet_phis, _jet_ms, 
                       _lep_pts, _lep_etas, _lep_phis, _lep_ms, 
                       _lep_ids, _lep_charges,   
                       _weight,  _is_truth,
                       _lep_d0raws,
                       _lep_d0sigs,
                       _lep_z0sinthetas,
                       _lep_isolations,
                       _jet_jvts);
        }
        float findAngleThetaStar(evt event){ return event.res.angles.theta_star; }
        float findAnglePhiStar(evt event){ return event.res.angles.phi_star; };
        float findAnglePhi0(evt event){ return event.res.angles.phi; };
        float findAnglePhi1(evt event){ return event.res.angles.phi1; };
        float findAngleTheta1(evt event){ return event.res.angles.theta1; };
        float findAngleTheta2(evt event){ return event.res.angles.theta2; };

        float findMJJ(evt event){ return event.res.m_jj; }
        float findDPhiJJ(evt event){ return event.res.dphi_jj;}
        float findDYJJ(evt event){ return event.res.dy_jj;}
        float findM4l(evt event){return event.res.m_4l;}
        float findMZ1(evt event){return event.res.m_z1;}
        float findMZ2(evt event){return event.res.m_z2;}
        float findPtBa(evt event){return event.res.pt_bala;}
        float findCent(evt event){return event.res.centra;}
        bool   passCut(evt event){return event.pass_cut;}
        bool   passSR(evt event){return event.pass_SR;}
        bool   passCR(evt event){return event.pass_CR;}
        bool   pass_sin(evt event) {return event.pass_sin_par;}
        bool   pass_j1j2(evt event) {return event.pass_j1j2;}
        bool   pass_z1z2(evt event) {return event.pass_z1z2;}
}

