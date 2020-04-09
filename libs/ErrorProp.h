#include <iostream>
#include <stdlib.h>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TF1.h"

#ifndef _ERRORPROP_
#define _ERRORPROP_
//6 Error Propagation functions
    double ErrAPlusB(double ErrA, double ErrB, double CorreAB)
    {
        double Err2 = pow(ErrA,2)+pow(ErrB,2)+ CorreAB * 2*ErrA*ErrB; 
        return TMath::Sqrt(Err2);
    }
    double ErrAMinusB(double ErrA, double ErrB, double  CorreAB)
    {
        double Err2 = pow(ErrA,2) + pow(ErrB,2) - CorreAB * 2*ErrA*ErrB; 
        return TMath::Sqrt(Err2);
    }
    double ErrAMultB(double C, double A, double B, double ErrA, double ErrB,double  CorreAB)
    {
        double f2 = C*C;
        double AErrA2 = pow(ErrA/A,2);
        double BErrB2 = pow(ErrB/B,2);
        double ErrAB = ErrA*ErrB/(B*A);
        double Err2 = f2*(AErrA2 + BErrB2 + CorreAB*2*ErrAB);
        return TMath::Sqrt(Err2);
    }
    double ErrADiviB(double C, double A, double B, double ErrA, double ErrB, double  CorreAB)
    {
        double f2 = C*C;
        double AErrA2 = pow(ErrA/A,2);
        double BErrB2 = pow(ErrB/B,2);
        double ErrAB = ErrA/A*ErrB/B;
        double Err2 = f2*(AErrA2+BErrB2-CorreAB * 2 * ErrAB);
        return TMath::Sqrt(Err2);
    }
    

    void GetSigCov(double *ComCov, double *SigCov, int NumComPar, int NumSigPar)
    {
        for( int i=0; i<NumSigPar; i++)  //column
        {
            for (int j = 0; j<NumSigPar;j++) //row
            {
                SigCov[i*NumSigPar+j] = ComCov[i*NumComPar+j];
            }
        }
    }
#endif