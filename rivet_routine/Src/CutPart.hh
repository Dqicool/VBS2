#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Math/Vector4.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include <cmath>
#define Z_MASS 91.1876*Rivet::GeV

class CP_angles{
    public:
        double theta_star;
        double phi_star;
        double phi;
        double phi1;
        double theta1;
        double theta2;
};

Rivet::Particles lepSel(Rivet::Particles lepton)
{
    auto siz = lepton.size();
    Rivet::Particles ret;

    for (uint i = 0;i < siz; i++)
    {
        //if (lepton_isTight[i] == 0)
        {
            if(lepton[i].abspid() == 11){
                if(lepton[i].pt() > 7*Rivet::GeV && lepton[i].abseta() < 2.47){
                    ret.push_back(lepton[i]);
                }
            }
            else if(lepton[i].abspid() == 13)
            {
                if(lepton[i].pt() > 7*Rivet::GeV && lepton[i].abseta() < 2.7){
                    ret.push_back(lepton[i]);
                }
            }
        }
    }
    
    return ret;
}

Rivet::FourMomenta jetSel(Rivet::Jets jet)
{
    auto siz = jet.size(); 
    Rivet::FourMomenta ret;
    for (uint i = 0; i < siz; i++){
        if (jet[i].abseta() < 2.4 && jet[i].pt() > 30*Rivet::GeV){
            ret.push_back(jet[i].momentum());
        }
        else if(jet[i].abseta() < 4.5 && jet[i].pt() > 40*Rivet::GeV){
            ret.push_back(jet[i].momentum());
        }
    } 
    return ret;
}

double getDeltaR(Rivet::FourMomentum a,Rivet::FourMomentum b){
    return std::sqrt((a.eta()-b.eta())*(a.eta()-b.eta()) + (a.phi()-b.phi())*(a.phi()-b.phi()));
}


bool jetsPairECriterion(Rivet::FourMomenta a, Rivet::FourMomenta b){
    double Ea = 0, Eb = 0;
    for(Rivet::FourMomentum ai : a){
        Ea += ai.E();
    }
    for(Rivet::FourMomentum bi : b){
        Eb += bi.E();
    }
    return Ea > Eb;
}

bool jetsPtCriterion(Rivet::FourMomentum a, Rivet::FourMomentum b){
    return a.pt() > b.pt();
}

Rivet::FourMomenta getJ1J2(Rivet::FourMomenta jets){
    Rivet::FourMomenta ret{};
    uint size = jets.size();
    if (size < 2) {}
    else{
        std::vector<Rivet::FourMomenta> jets_pairs{};
        for (uint i = 0; i < size; i++)
        {
            for (uint j = i + 1; j < size; j++)
            {
                if(jets[i].rapidity() * jets[j].rapidity() < 0)
                {
                    Rivet::FourMomenta tmp_pairs{jets[i], jets[j]};
                    jets_pairs.push_back(tmp_pairs);
                }
            }
        }
        if(jets_pairs.size() > 0)
        {
            std::sort(jets_pairs.begin(), jets_pairs.end(), jetsPairECriterion);
            std::sort(jets_pairs[0].begin(), jets_pairs[0].end(), jetsPtCriterion);
            ret = jets_pairs[0];
            if(jets_pairs[0][0].pt() < 40*Rivet::GeV || jets_pairs[0][1].pt() < 30*Rivet::GeV) {ret = {};}
        }
    }
    return ret;
}

std::vector<std::vector<Rivet::Particles>> getLeptonPairsComb(Rivet::Particles leps){
    std::vector<std::vector<Rivet::Particles>> ret;
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
                                    auto delR_ij = getDeltaR(leps[i], leps[j]);
                                    auto delR_kl = getDeltaR(leps[k], leps[l]);
                                    if (delR_ij > 0.2 && delR_kl > 0.2){
                                        if( leps[i].abspid() == leps[j].abspid() 
                                            && ( leps[i].momentum() +  leps[j].momentum()).mass() > 10 
                                            && ( leps[k].momentum() +  leps[l].momentum()).mass() > 10)
                                        {
                                            Rivet::Particles tmp1{leps[i], leps[j]};
                                            Rivet::Particles tmp2{leps[k], leps[l]};
                                            std::vector<Rivet::Particles> tmp{tmp1,tmp2};
                                            ret.push_back(tmp);
                                        }
                                        else if(leps[i].abspid() != leps[k].abspid() ){
                                            Rivet::Particles tmp1{leps[i], leps[j]};
                                            Rivet::Particles tmp2{leps[k], leps[l]};
                                            std::vector<Rivet::Particles> tmp{tmp1,tmp2};
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

bool lepsPairsMCriterion(Rivet::Particles a, Rivet::Particles b){
    double mass_dif_a = std::abs((a[0].momentum() + a[1].momentum()).mass() - Z_MASS);
    double mass_dif_b = std::abs((b[0].momentum() + b[1].momentum()).mass() - Z_MASS);
    return mass_dif_a < mass_dif_b;
}

bool lepsPairsCombMCriterion(std::vector<Rivet::Particles> pairs_a, std::vector<Rivet::Particles> pairs_b){
    double comb_mass_diff_a =   std::abs((pairs_a[0][0].momentum() + pairs_a[0][1].momentum()).mass() - Z_MASS) +
                                std::abs((pairs_a[1][0].momentum() + pairs_a[1][1].momentum()).mass() - Z_MASS);
    double comb_mass_diff_b =   std::abs((pairs_b[0][0].momentum() + pairs_b[0][1].momentum()).mass() - Z_MASS) +
                                std::abs((pairs_b[1][0].momentum() + pairs_b[1][1].momentum()).mass() - Z_MASS);    
    return comb_mass_diff_a < comb_mass_diff_b;
}

std::vector<Rivet::Particles> getLepPair(std::vector<std::vector<Rivet::Particles>> lep_pairs_comb){
    std::vector<Rivet::Particles> ret{};
    if(lep_pairs_comb.size() >= 1){
        for(uint i = 0; i < lep_pairs_comb.size(); i++){
            std::sort(lep_pairs_comb[i].begin(), lep_pairs_comb[i].end(), lepsPairsMCriterion);
        }
        std::sort(lep_pairs_comb.begin(), lep_pairs_comb.end(), lepsPairsCombMCriterion);
        ret = lep_pairs_comb[0];
    }
    return ret;
}

Rivet::FourMomenta getZ1Z2(std::vector<Rivet::Particles> lep_pairs){
    Rivet::FourMomenta ret{};
    if(lep_pairs.size() >= 1){
        ret = {(lep_pairs[0][0].mom() + lep_pairs[0][1].mom()), 
                (lep_pairs[1][0].mom() + lep_pairs[1][1].mom())};
    }
    return ret;
}

double getPtBala (Rivet::FourMomenta j1j2, Rivet::FourMomenta z1z2){
    double ret = -999;
    if(j1j2.size() >= 2 && z1z2.size() >= 2){
        double zzjj_sys_pt = (z1z2[0] + z1z2[1] + j1j2[0] + j1j2[1]).pt(); 
        double zzjj_scaler_sum_pt =z1z2[0].pt()  + z1z2[1].pt() + j1j2[0].pt() + j1j2[1].pt();
        ret = std::abs(zzjj_sys_pt / zzjj_scaler_sum_pt);
    }
    return ret;
}

double getJJDelPhi(Rivet::FourMomenta j1j2){
    double ret = -999;
    if(j1j2.size() >= 2){
        if (j1j2[0].rapidity() > j1j2[1].rapidity())    
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

CP_angles getAngles(std::vector<Rivet::Particles> lep_pairs, Rivet::FourMomenta z1z2){
    auto starret = new CP_angles();
    auto ret = *starret;
    auto zz_rest = Rivet::LorentzTransform::mkFrameTransformFromBeta((z1z2[0] + z1z2[1]).betaVec());
    auto q11 = zz_rest.transform(lep_pairs[0][0].momentum());
    auto q12 = zz_rest.transform(lep_pairs[0][1].momentum());
    auto q21 = zz_rest.transform(lep_pairs[1][0].momentum());
    auto q22 = zz_rest.transform(lep_pairs[1][1].momentum());
    auto q1  = zz_rest.transform(z1z2[0]);
    auto q2  = zz_rest.transform(z1z2[1]);
    auto n1hat  = q11.vector3().cross(q12.vector3()) / std::sqrt((q11.vector3().cross(q12.vector3())).mod2());
    auto n2hat  = q21.vector3().cross(q22.vector3()) / std::sqrt((q21.vector3().cross(q22.vector3())).mod2());
    Rivet::Vector3 nzhat(0.0, 0.0, 1.0);
    auto nschat = nzhat.cross(q1.vector3()) / std::sqrt(nzhat.cross(q1.vector3()).mod2());
    ret.theta_star = q1.theta();
    ret.phi_star   = q1.phi(Rivet::MINUSPI_PLUSPI);
    ret.phi  = q1.vector3().dot(n1hat.cross(n2hat)) / std::abs(q1.vector3().dot(n1hat.cross(n2hat))) * (std::acos(-n1hat.dot(n2hat)));
    ret.phi1 = q1.vector3().dot(n1hat.cross(nschat)) / std::abs(q1.vector3().dot(n1hat.cross(nschat))) * (std::acos(n1hat.dot(nschat)));
    ret.theta1 = std::acos(-(q2.vector3().dot(q11.vector3())) / std::sqrt(q2.vector3().mod2()) / std::sqrt(q11.vector3().mod2()));
    ret.theta2 = std::acos(-(q1.vector3().dot(q21.vector3())) / std::sqrt(q1.vector3().mod2()) / std::sqrt(q21.vector3().mod2()));
    return ret;
}