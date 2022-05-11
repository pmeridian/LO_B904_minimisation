#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom.h"
#include "TMath.h"
#include "TError.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TTree.h"
#include <iostream>

TH1F* h_fit=0;
TTree* ori_tree=0;
int min_bin=0;
int max_bin=0;


TH1F* smeared_scaled_histo(TTree* tree,TH1F* hFit, const double* pars)
{
  TH1F* hOut=(TH1F*)hFit->Clone("hOut");
  hOut->Reset();

  float charge_sig[100];
  tree->SetBranchAddress("charge_sig",&charge_sig);
  int n_entries=tree->GetEntries();
  /* std::cout << n_entries << std::endl; */
  std::cout << pars[0] << "," << pars[1] << "," << pars[2] << std::endl;
  for(int j=0;j<50;++j)
    {
      for(int i=0;i<n_entries;++i)
	{
	  tree->GetEntry(i);
	  float val=charge_sig[0]*pars[1]+gRandom->Gaus(0,pars[2]);
	  hOut->Fill(val,pars[0]);
	}
    }
  hOut->Smooth();
  hOut->Scale(1/50.);
  return hOut;
}

double chi2_func(const double* pars)
{
  double chi2=0;
  TH1F* hOut=smeared_scaled_histo(ori_tree,h_fit,pars);
  for(int i=min_bin;i<=max_bin;++i)
    {
      chi2+=TMath::Power(h_fit->GetBinContent(i)-hOut->GetBinContent(i),2)/(h_fit->GetBinError(i)*h_fit->GetBinError(i)+hOut->GetBinError(i)*hOut->GetBinError(i));
    }
  std::cout << chi2 << std::endl;
  delete hOut;
  return chi2;
}

int NumericalMinimization(TTree* tree_ori,TH1F* fit_h,double xMin,double xMax)
{
  ori_tree=tree_ori;
  h_fit=fit_h;
  min_bin=fit_h->FindBin(xMin);
  max_bin=fit_h->FindBin(xMax);
  
  ROOT::Math::Minimizer* minimum =
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  
  //  minimum->SetPrecision(0.1); 
  minimum->SetPrintLevel(1);
  ROOT::Math::Functor f(&chi2_func,3);
  double step[3] = {1,0.1,100};
  minimum->SetFunction(f);
  minimum->SetStrategy(0);
  minimum->SetTolerance(10);
  minimum->SetVariable(0,"norm",5.8, step[0]);
  minimum->SetVariableLimits(0, 5., 7.);
  minimum->SetVariable(1,"scale",0.77, step[1]);
  minimum->SetVariableLimits(1, 0.6, 1.2);
  minimum->SetVariable(2,"smear",230, step[2]);
  minimum->SetVariableLimits(2, 150, 250);
  minimum->Minimize();

  return 0;
}
