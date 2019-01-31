// Arguments:
// KevinRecipe.cpp("12C_FSI_ON.root", "12C_FSI_OFF.root", "40Ar_FSI_ON.root");
// Where the arguments points to files that have been processed by fsianal.cpp (or runit.sh)
const double Shift = -999;

double WeightSecondOrder(double xvar, double a0, double a1, double a2) {

  double w=a0+a1*xvar+a2*xvar*xvar;
  //if (w < 0) return 0.0;
  return a0+a1*xvar+a2*xvar*xvar;
}

double SolveA0(double alpha, double beta, double gamma, double delta, double eps, double mu, double sig) {
  double a0 = 0.0;

  double alpha2 = alpha*alpha;
  double beta2 = beta*beta;
  double gamma2 = gamma*gamma;
  double delta2 = delta*delta;
  double eps2 = eps*eps;
  double mu2 = mu*mu;
  double sig2 = sig*sig;

  double numerator = -1.0*(-1.0*alpha*delta2 + alpha*eps*gamma - beta*eps*mu + delta*gamma*mu + 2.0*beta*delta*mu2 - alpha*beta*delta*mu2 - 2.0*gamma2*mu2 + alpha*gamma2*mu2 + beta*delta*sig - gamma2*sig);
  double denominator = alpha*delta2 + beta2*eps - 2*beta*delta*gamma - alpha*eps*gamma + gamma*gamma2;
  if (denominator == 0) return 0.0;

  return numerator/denominator;
}

double SolveA1(double alpha, double beta, double gamma, double delta, double eps, double mu, double sig) {
  double a1 = 0.0;

  double alpha2 = alpha*alpha;
  double beta2 = beta*beta;
  double gamma2 = gamma*gamma;
  double delta2 = delta*delta;
  double eps2 = eps*eps;
  double mu2 = mu*mu;
  double sig2 = sig*sig;

  double numerator = -1.0*(-1.0*alpha*beta*eps + alpha*delta*gamma + alpha*eps*mu - gamma2*mu - 2.0*alpha*delta*mu2 + alpha2*delta*mu2 + 2.0*beta*gamma*mu2 - alpha*beta*gamma*mu2 - alpha*delta*sig + beta*gamma*sig);
  double denominator = alpha*delta2 + beta2*eps - 2.0*beta*delta*gamma - alpha*eps*gamma + gamma*gamma2;
  if (denominator == 0) return 0;

  return numerator/denominator;
}

double SolveA2(double alpha, double beta, double gamma, double delta, double eps, double mu, double sig) {
  double a2 = 0.0;

  double alpha2 = alpha*alpha;
  double beta2 = beta*beta;
  double gamma2 = gamma*gamma;
  double delta2 = delta*delta;
  double eps2 = eps*eps;
  double mu2 = mu*mu;
  double sig2 = sig*sig;

  double numerator = -1.0*(alpha*beta*delta - alpha*gamma2 - alpha*delta*mu + beta*gamma*mu - 2*beta2*mu2 + alpha*beta2*mu2 + 2*alpha*gamma*mu2 - alpha2*gamma*mu2 - beta2*sig + alpha*gamma*sig);
  double denominator = alpha*delta2 + beta2*eps - 2*beta*delta*gamma - alpha*eps*gamma + gamma*gamma2;
  if (denominator == 0) return 0;

  return numerator/denominator;
}

double SolveA0_New(double alpha, double beta, double gamma, double delta, double eps, double mu, double sig) {
  double alpha2 = alpha*alpha;
  double beta2 = beta*beta;
  double gamma2 = gamma*gamma;
  double delta2 = delta*delta;
  double eps2 = eps*eps;
  double mu2 = mu*mu;
  double sig2 = sig*sig;

  double numerator = -1.0*(-1.0*alpha*delta2 + alpha*eps*gamma - alpha*beta*eps*mu + alpha*delta*gamma*mu + alpha*beta*delta*mu2 - alpha*gamma2*mu2 + alpha*beta*delta*sig - alpha*gamma2*sig);
  double denominator = alpha*delta2 + beta2*eps - 2*beta*delta*gamma - alpha*eps*gamma + gamma*gamma2;
  if (denominator == 0) {
    std::cout << "AAAARGH A0 INFINITE" << std::endl;
    return -999.99;
  }

  return numerator/denominator;
}

double SolveA1_New(double alpha, double beta, double gamma, double delta, double eps, double mu, double sig) {

  double alpha2 = alpha*alpha;
  double beta2 = beta*beta;
  double gamma2 = gamma*gamma;
  double delta2 = delta*delta;
  double eps2 = eps*eps;
  double mu2 = mu*mu;
  double sig2 = sig*sig;

  double numerator = -1.0*(-1.0*alpha*beta*eps + alpha*delta*gamma + alpha2*eps*mu - alpha*gamma2*mu - alpha2*delta*mu2 + alpha*beta*gamma*mu2 - alpha2*delta*sig + alpha*beta*gamma*sig);
  double denominator = alpha*delta2 + beta2*eps - 2*beta*delta*gamma - alpha*eps*gamma + gamma*gamma2;
  if (denominator == 0) {
    std::cout << "AAAARGH A1 INFINITE" << std::endl;
    return -999.99;
  }

  return numerator/denominator;
}

double SolveA2_New(double alpha, double beta, double gamma, double delta, double eps, double mu, double sig) {

  double alpha2 = alpha*alpha;
  double beta2 = beta*beta;
  double gamma2 = gamma*gamma;
  double delta2 = delta*delta;
  double eps2 = eps*eps;
  double mu2 = mu*mu;
  double sig2 = sig*sig;

  double numerator = -1.0*(alpha*beta*delta - alpha*gamma2 - alpha2*delta*mu + alpha*beta*gamma*mu - alpha*beta2*mu2 + alpha2*gamma*mu2 - alpha*beta2*sig + alpha2*gamma*sig);
  double denominator = alpha*delta2 + beta2*eps - 2*beta*delta*gamma - alpha*eps*gamma + gamma*gamma2;
  if (denominator == 0) {
    std::cout << "AAAARGH A1 INFINITE" << std::endl;
    return -999.99;
  }

  return numerator/denominator;
}

TH2D** SolveSystem(TH2D* MuHist, TH2D* VarHist, TH2D* ZeroHist, TH2D* FirstHist, TH2D* SecondHist, TH2D* ThirdHist, TH2D* FourthHist) {

  TH2D *A0Hist = (TH2D*)MuHist->Clone("A0Hist");
  A0Hist->Reset();
  A0Hist->SetDirectory(0);
  A0Hist->SetNameTitle("A0Hist","A0Hist");

  TH2D *A1Hist = (TH2D*)MuHist->Clone("A1Hist");
  A1Hist->Reset();
  A1Hist->SetDirectory(0);
  A1Hist->SetNameTitle("A1Hist","A1Hist");

  TH2D *A2Hist = (TH2D*)MuHist->Clone("A2Hist");
  A2Hist->Reset();
  A2Hist->SetDirectory(0);
  A2Hist->SetNameTitle("A2Hist","A2Hist");

  // Loop over each q0, q3 bin and solve for a0, a1 and a2
  for (int i = 0; i < MuHist->GetXaxis()->GetNbins(); ++i) {
    for (int j = 0; j < MuHist->GetYaxis()->GetNbins(); ++j) {
      // Mean
      double mu = MuHist->GetBinContent(i+1,j+1);
      // Variance
      double var = VarHist->GetBinContent(i+1,j+1);

      // Integral[x^0 * f(x)dx]
      double zeroth = ZeroHist->GetBinContent(i+1, j+1);
      // Integral[x^1 * f(x)dx]
      double first = FirstHist->GetBinContent(i+1, j+1);
      // Integral[x^2 * f(x)dx]
      double second = SecondHist->GetBinContent(i+1, j+1);
      // Integral[x^3 * f(x)dx]
      double third = ThirdHist->GetBinContent(i+1, j+1);
      // Integral[x^4 * f(x)dx]
      double fourth = FourthHist->GetBinContent(i+1, j+1);

      // Now solve the system
      //double a0 = SolveA0(zeroth, first, second, third, fourth, mu, var);
      //double a1 = SolveA1(zeroth, first, second, third, fourth, mu, var);
      //double a2 = SolveA2(zeroth, first, second, third, fourth, mu, var);
      double a0 = SolveA0_New(zeroth, first, second, third, fourth, mu, var);
      double a1 = SolveA1_New(zeroth, first, second, third, fourth, mu, var);
      double a2 = SolveA2_New(zeroth, first, second, third, fourth, mu, var);

      //std::cout << a0 << " " << a1 << " " << a2 << std::endl;
      //std::cout << zeroth << " " << first << " " << second << " " << third << " " << fourth << " " << mu << " " << var << std::endl;

      A0Hist->SetBinContent(i+1, j+1, a0);
      A1Hist->SetBinContent(i+1, j+1, a1);
      A2Hist->SetBinContent(i+1, j+1, a2);
    }
  }

  SetLimits(A0Hist);
  SetLimits(A1Hist);
  SetLimits(A2Hist);

  TH2D **arr = new TH2D[3];
  arr[0] = A0Hist;
  arr[1] = A1Hist;
  arr[2] = A2Hist;

  return arr;
}

void SetLimits(TH2D* Hist) {
  // find the maximum
  double maximum = 0.0;
  double minimum = 999;
  for (int i = 0; i < Hist->GetXaxis()->GetNbins(); ++i) {
    for (int j = 0; j < Hist->GetYaxis()->GetNbins(); ++j) {
      double tempmax = Hist->GetBinContent(i+1, j+1);
      if (tempmax > maximum) maximum = tempmax;
      if (tempmax < minimum) minimum = tempmax;
    }
  }
  Hist->SetMaximum(maximum);
  Hist->SetMinimum(minimum);
}

TH2D* SubtractAndAddHist(std::string filename1, std::string filename2, std::string filename3, double shift, std::string HistName) {

  std::cout << "12C FSI ON, 12C FSI OFF, 40Ar FSI ON: " << filename1 << ", " << filename2 << ", " << filename3 << std::endl;

  TFile *file1 = new TFile(filename1.c_str(), "OPEN");
  TFile *file2 = new TFile(filename2.c_str(), "OPEN");
  TFile *file3 = new TFile(filename3.c_str(), "OPEN");

  // One is FSI on C
  TH2D *one = (TH2D*)file1->Get(HistName.c_str())->Clone(Form("%s_One", HistName.c_str()));
  one->SetDirectory(0);
  //one->SetName((std::string(one->GetName())+"_One").c_str());
  //SetLimits(one);

  // Two is FSI off C
  TH2D *two = (TH2D*)file2->Get(HistName.c_str())->Clone(Form("%s_Two", HistName.c_str()));
  two->SetDirectory(0);
  //two->SetName((std::string(two->GetName())+"_Two").c_str());
  //SetLimits(two);

  // Three is FSI on Ar
  TH2D *three = (TH2D*)file3->Get(HistName.c_str())->Clone(Form("%s_Three", HistName.c_str()));
  three->SetDirectory(0);
  //two->SetName((std::string(two->GetName())+"_Two").c_str());
  //SetLimits(two);

  TH2D *NewHist = one->Clone(Form("NewHist_%s", HistName.c_str()));
  NewHist->SetDirectory(0);
  NewHist->Reset();

  for (int i = 0; i < one->GetXaxis()->GetNbins(); ++i) {
    for (int j = 0; j < one->GetYaxis()->GetNbins(); ++j) {
      // This is for every q0, q3 bin
      // For 12C FSI on
      double content = one->GetBinContent(i+1, j+1);
      // For 12C FSI off
      double content2 = two->GetBinContent(i+1, j+1);
      // For 40Ar FSI on
      double content3 = three->GetBinContent(i+1, j+1);

      double diff = content-content2;
      //one->SetBinContent(i+1, j+1, content2+diff*shift);
      NewHist->SetBinContent(i+1, j+1, content3 + diff*shift);
    }
  }
  //SetLimits(one);
  NewHist->SetTitle(Form("%s_Subtract%.2f_%s_Add_%s", filename1.c_str(), shift, filename2.c_str(), filename3.c_str()));

  //one->SetName(Form("%s_Shift%.2f", one->GetName(), shift));
  NewHist->GetZaxis()->SetTitle(Form("%s difference (^{12}C FSI on-^{12}CFSI off)*%.2f + ^{40}Ar FSI On", NewHist->GetZaxis()->GetTitle(), shift));

  delete two;
  delete one;
  delete three;
  file1->Close();
  file2->Close();
  file3->Close();

  return NewHist;
}

void KevinRecipe(std::string file1, std::string file2, std::string file3) {

  std::cout << "FSI ON C, FSI OFF C, FSI ON AR" << std::endl;
  TH1::AddDirectory(false);

  std::string MeanName = "MeanPlot";
  double meanshift = 0.3;
  // Get the shifted mean histogram
  TH2D *MeanHist = SubtractAndAddHist(file1, file2, file3, meanshift, MeanName);

  std::string VarName = "RMSPlot";
  double varshift = 0.3;
  // Get the shifted variance histogram
  TH2D *VarHist = SubtractAndAddHist(file1, file2, file3, varshift, VarName);

  //std::string SkewName = "SkewPlot";
  //TH2D *SkewHist = SubtractAndAddHist(file1, file2, Shift, SkewName);

  // Open the FSI ON C file
  TFile *FileFSIon = new TFile(file1.c_str(), "OPEN");

  // Open the FSI OFF C file
  TFile *FileFSIoff = new TFile(file2.c_str(), "OPEN");

  // Open the FSI ON Ar file
  TFile *FileFSIonAr = new TFile(file3.c_str(), "OPEN");

  ///////////////////
  std::string ZeroName = "ZerothPlot";
  TH2D *ZeroHist = (TH2D*)FileFSIonAr->Get(ZeroName.c_str())->Clone("ZerothPlot");
  ZeroHist->SetDirectory(0);
  std::string FirstName = "FirstPlot";
  TH2D *FirstHist = (TH2D*)FileFSIonAr->Get(FirstName.c_str())->Clone("FirstPlot");
  FirstHist->SetDirectory(0);
  std::string SecondName = "SecondPlot";
  TH2D *SecondHist = (TH2D*)FileFSIonAr->Get(SecondName.c_str())->Clone("SecondPlot");
  SecondHist->SetDirectory(0);
  std::string ThirdName = "ThirdPlot";
  TH2D *ThirdHist = (TH2D*)FileFSIonAr->Get(ThirdName.c_str())->Clone("ThirdPlot");
  ThirdHist->SetDirectory(0);
  std::string FourthName = "FourthPlot";
  TH2D *FourthHist = (TH2D*)FileFSIonAr->Get(FourthName.c_str())->Clone("FourthPlot");
  FourthHist->SetDirectory(0);

  // Now have the difference in Mean (0th moment), Variance (1st moment), Skew (2nd moment)
  //TFile *Output = new TFile(Form("%s_Subtract%.2f_%s", file1.c_str(), Shift, file2.c_str()), "RECREATE");
  std::string canvname = file1+"FULL_validation_new_";
  stringstream ss;
  ss << "MeanShift_" << meanshift << "_VarShift_" << varshift;
  canvname += ss.str();

  TFile *Output = new TFile((canvname+".root").c_str(), "RECREATE");
  Output->cd();
  MeanHist->Write();
  VarHist->Write();

  // Solve the system
  TH2D **arr = SolveSystem(MeanHist, VarHist, ZeroHist, FirstHist, SecondHist, ThirdHist, FourthHist);
  for (int i = 0; i < 3; ++i) {
    arr[i]->Write();
  }

  // Now do the validaiton
  // We have the A0, A1 and A2 in arr

  // Get the Carbon FSI On
  TH3D *FSIPlot_C = (TH3D*)FileFSIon->Get("q0q3Eav")->Clone("C12FSIOn_q0q3Eav");
  FSIPlot_C->SetDirectory(0);

  TH3D *NoFSIPlot_C = (TH3D*)FileFSIoff->Get("q0q3Eav")->Clone("C12FSIOff_q0q3Eav");
  NoFSIPlot_C->SetDirectory(0);

  // Get the Argon FSI On
  TH3D *FSIPlot = (TH3D*)FileFSIonAr->Get("q0q3Eav")->Clone("Ar40FSIOn_q0q3Eav");
  FSIPlot->SetDirectory(0);

  TCanvas *canv = new TCanvas("MetaData", "MetaData", 1024, 1024);
  canv->SetRightMargin(canv->GetRightMargin()*0.8);
  canv->SetLeftMargin(canv->GetLeftMargin()*1.1);
  canv->SetBottomMargin(canv->GetBottomMargin()*0.8);
  canv->Print((canvname+".pdf[").c_str());

  //TLegend *leg = new TLegend(0.05, 0.4, 0.4, 0.95);
  TLegend *leg = new TLegend(0., 0., 1., 1.);
  leg->AddEntry("", Form("^{12}C FSI On: %s", file1.c_str()), "");
  leg->AddEntry("", Form("^{12}C FSI Off: %s", file2.c_str()), "");
  leg->AddEntry("", Form("^{40}Ar FSI On: %s", file3.c_str()), "");
  leg->AddEntry("", Form("Mean Shift: %.2f", meanshift), "");
  leg->AddEntry("", Form("Var Shift: %.2f", varshift), "");
  canv->Clear();
  leg->Draw();
  canv->Print((canvname+".pdf").c_str());
  canv->Write();
  delete leg;

  canv->Clear();

  TH3D* WeightPlot = FSIPlot->Clone("WeightPlot");
  WeightPlot->SetDirectory(0);
  WeightPlot->Reset();
  WeightPlot->SetTitle(Form("FSI Weighting, #mu_{shift}=%.2f #sigma^{2}_{shift}=%.2f", meanshift, varshift));

  int NegativeCounts = 0;
  int TotalCounts = 0;
  // Loop over each bin
  for (int i = 0; i < NoFSIPlot_C->GetXaxis()->GetNbins(); ++i) {
    for (int j = 0; j < NoFSIPlot_C->GetYaxis()->GetNbins(); ++j) {
      TH1D *FSIPlot1D_C = FSIPlot_C->ProjectionZ(Form("proj_fsiC_%i_%i", i+1, j+1), i+1, i+1, j+1, j+1);
      TH1D *FSIPlot1D_Ar = FSIPlot->ProjectionZ(Form("proj_fsiAr_%i_%i", i+1, j+1), i+1, i+1, j+1, j+1);
      TH1D *NoFSIPlot1D_C = NoFSIPlot_C->ProjectionZ(Form("proj_nofsiC_%i_%i", i+1, j+1), i+1, i+1, j+1, j+1);

      double q3lo = FSIPlot->GetXaxis()->GetBinLowEdge(i+1);
      double q3hi = FSIPlot->GetXaxis()->GetBinLowEdge(i+2);
      double q0lo = FSIPlot->GetYaxis()->GetBinLowEdge(j+1);
      double q0hi = FSIPlot->GetYaxis()->GetBinLowEdge(j+2);
      FSIPlot1D_Ar->SetTitle(Form("E_{av}/q_{0}: q_{0}=%.2f-%.2f, |q_{3}|=%.2f-%.2f", q0lo, q0hi, q3lo, q3hi));

      TH1D *FSIPlot1D_Weighted = (TH1D*)FSIPlot1D_C->Clone(Form("%s_clone", FSIPlot1D_C->GetName()));
      FSIPlot1D_Weighted->Reset();
      // Get the coefficients for this q0,q3
      double a0 = arr[0]->GetBinContent(i+1, j+1);
      double a1 = arr[1]->GetBinContent(i+1, j+1);
      double a2 = arr[2]->GetBinContent(i+1, j+1);

      TH1D *WeightingFunction = (TH1D*)FSIPlot1D_C->Clone(Form("Weighting_%i_%i", i+1, j+1));
      WeightingFunction->Reset();
      WeightingFunction->GetYaxis()->SetTitle("Weighting function");
      WeightingFunction->SetTitle(Form("Weighting %s", FSIPlot1D_Ar->GetTitle()));

      // Weighting is a function of Eav/q0 -> Loop over axis
      for (int k = 0; k < FSIPlot1D_Weighted->GetXaxis()->GetNbins(); ++k) {
        TotalCounts++;
        double xvar = FSIPlot1D_Weighted->GetXaxis()->GetBinCenter(k+1);
        double weight = WeightSecondOrder(xvar, a0, a1, a2);
        if (a0 == -999.99 || a1 == -999.99 || a2 == -999.99) {
          weight = 1.0;
        }
        if (weight < 0) NegativeCounts++;

        double content = FSIPlot1D_Ar->GetBinContent(k+1)*weight;
        if (weight == 0 && FSIPlot1D_C->Integral() > 0 && FSIPlot1D_Ar->Integral() > 0) {
          std::cout << "Found zero weight: " << std::endl;
          std::cout << "q0, q3: " << q0lo << "-" << q0hi << ", " << q3lo << "-" << q3hi << std::endl;
          std::cout << "Eav/q0: " << xvar << " a0: " << a0 << " a1: " << a1 << " a2: " << a2 << std::endl;
          std::cout << "FSIPlot1D_C: " << FSIPlot1D_C->Integral() << std::endl;
          std::cout << "FSIPlot1D_Ar: " << FSIPlot1D_Ar->Integral() << std::endl;
          std::cout << "FSIPlot1D_Weighted: " << FSIPlot1D_Weighted->Integral() << std::endl;
        }
        WeightingFunction->SetBinContent(k+1, weight);
        WeightPlot->SetBinContent(i+1, j+1, k+1, weight);
        FSIPlot1D_Weighted->SetBinContent(k+1, content);
      }

      FSIPlot1D_C->SetLineWidth(2);
      NoFSIPlot1D_C->SetLineWidth(2);
      FSIPlot1D_Ar->SetLineWidth(2);
      FSIPlot1D_Weighted->SetLineWidth(2);

      FSIPlot1D_C->SetLineColor(kBlack);
      NoFSIPlot1D_C->SetLineColor(kBlack);
      NoFSIPlot1D_C->SetLineStyle(kDashed);
      FSIPlot1D_Ar->SetLineWidth(2);
      FSIPlot1D_Ar->SetLineColor(kBlue);
      FSIPlot1D_Weighted->SetLineColor(kRed);
      FSIPlot1D_Weighted->SetLineStyle(kDashed);

      if (FSIPlot1D_C->Integral() == 0 && 
          FSIPlot1D_Ar->Integral() == 0 && 
          FSIPlot1D_Weighted->Integral() == 0) continue;

      double mutarget=(FSIPlot1D_C->GetMean()-NoFSIPlot1D_C->GetMean())*meanshift+ FSIPlot1D_Ar->GetMean();
      double vartarget = (FSIPlot1D_C->GetRMS()*FSIPlot1D_C->GetRMS()- NoFSIPlot1D_C->GetRMS()*NoFSIPlot1D_C->GetRMS())*varshift+ FSIPlot1D_Ar->GetRMS()*FSIPlot1D_Ar->GetRMS();

      // Draw a legend
      TLegend *leg = new TLegend(0.3, 0.5, 0.9, 0.9);
      leg->AddEntry(FSIPlot1D_Weighted, Form("#splitline{FSI weighted: #int=%.3e, #mu=%.3f, #sigma^{2}=%.3f}{Target: #int=%.3e, #mu=%.3f, #sigma^{2}=%.3f}",
                                                                                                  FSIPlot1D_Weighted->Integral(), 
                                                                                                  FSIPlot1D_Weighted->GetMean(), 
                                                                                                  FSIPlot1D_Weighted->GetRMS()*FSIPlot1D_Weighted->GetRMS(), 
                                                                                                  FSIPlot1D_Ar->Integral(),
                                                                                                  mutarget,
                                                                                                  vartarget), "l");

      leg->AddEntry(FSIPlot1D_Ar, Form("w/ FSI ^{40}Ar #int=%.3e, #mu=%.3f, #sigma^{2}=%.3f",  FSIPlot1D_Ar->Integral(), 
                                                                                    FSIPlot1D_Ar->GetMean(), 
                                                                                    FSIPlot1D_Ar->GetRMS()*
                                                                                    FSIPlot1D_Ar->GetRMS()), 
                                                                                    "l");

      leg->AddEntry(FSIPlot1D_C,Form("w/ FSI ^{12}C #int=%.3e, #mu=%.3f, #sigma^{2}=%.3f",  FSIPlot1D_C->Integral(), 
                                                                                      FSIPlot1D_C->GetMean(), 
                                                                                      FSIPlot1D_C->GetRMS()*
                                                                                      FSIPlot1D_C->GetRMS()), 
                                                                                      "l"); 

      leg->AddEntry(NoFSIPlot1D_C,Form("w/o FSI ^{12}C #int=%.3e, #mu=%.3f, #sigma^{2}=%.3f",  NoFSIPlot1D_C->Integral(), 
                                                                                      NoFSIPlot1D_C->GetMean(), 
                                                                                      NoFSIPlot1D_C->GetRMS()*
                                                                                      NoFSIPlot1D_C->GetRMS()), 
                                                                                      "l"); 

      // Get the maximum
      double maximum = FSIPlot1D_C->GetMaximum();
      if (maximum < FSIPlot1D_Ar->GetMaximum()) maximum = FSIPlot1D_Ar->GetMaximum();
      if (maximum < FSIPlot1D_Weighted->GetMaximum()) maximum = FSIPlot1D_Weighted->GetMaximum();
      if (maximum < NoFSIPlot1D_C->GetMaximum()) maximum = NoFSIPlot1D_C->GetMaximum();
      maximum *= 1.3;
      // Find minimum
      double minimum = 0.0;
      if (minimum > FSIPlot1D_Ar->GetMinimum()) minimum = FSIPlot1D_Ar->GetMinimum();
      if (minimum > FSIPlot1D_Weighted->GetMinimum()) minimum = FSIPlot1D_Weighted->GetMinimum();
      if (minimum > NoFSIPlot1D_C->GetMinimum()) minimum = NoFSIPlot1D_C->GetMinimum();
      minimum *= 1.3;
      FSIPlot1D_Ar->GetYaxis()->SetRangeUser(minimum, maximum);

      canv->Clear();
      FSIPlot1D_Ar->Draw();
      FSIPlot1D_C->Draw("same");
      NoFSIPlot1D_C->Draw("same");
      FSIPlot1D_Weighted->Draw("same");
      leg->Draw("same");
      canv->Print((canvname+".pdf").c_str());

      canv->Clear();
      WeightingFunction->Draw();
      //leg->Draw("same");
      TLegend *leg2 = new TLegend(0.12, 0.8, 0.88, 0.92);
      leg2->AddEntry(WeightingFunction, Form("w(x)=%.2f + %.2fx + %.2fx^{2}", a0, a1, a2), "");
      leg2->Draw("same");
      canv->Print((canvname+".pdf").c_str());

      delete FSIPlot1D_Weighted;
      delete FSIPlot1D_C;
      delete NoFSIPlot1D_C;
      delete FSIPlot1D_Ar;
      delete WeightingFunction;
      delete leg;
      delete leg2;
    }
  }
  Output->cd();
  WeightPlot->Write();
  canv->Print((canvname+".pdf]").c_str());

  std::cout << "Found " << NegativeCounts << "/" << TotalCounts << " (" << double(NegativeCounts)/double(TotalCounts)*100.0 << "%) " << "bad events" << std::endl;

  Output->Close();
  std::cout << "Wrote " << Output->GetName() << std::endl;
  std::cout << "Wrote " << canvname << ".pdf" << std::endl;
  std::cout << "Done" << std::endl;
}
