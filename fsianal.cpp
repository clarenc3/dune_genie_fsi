// Takes one argument
// Runs on a NUISANCE flat-tree format file
// Produces differential cross sections in q0, q3, Eav/q0
//
//
int MyPalette[255];
int kDarkRainBow = 55;

// Make the TPalette
void SetCustomPalette() {
  // Uncomment these to use the blue->white->red palette (good for correlation matrices)
  TColor::InitializeColors();
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.25, 1.00, 1.00, 0.50 };
  Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00 };
  Double_t blue[NRGBs]  = { 0.50, 1.00, 1.00, 0.25, 0.00 };

  int start = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont); 
  for (int i = 0; i < NCont; i++) {
    MyPalette[i] = start+i;
  }
}

double GetMeMyIntegral(TH1D *Hist, int order) {
  double sum = 0.0;
  for (int i = 0; i < Hist->GetXaxis()->GetNbins(); ++i) {
    /// Only go up until 1.2 Eav/q0
    double x = 1.0;
    // Multiple the order
    int counter = 0;
    for (int j = 0; j < order; ++j) {
      counter++;
      x *= Hist->GetXaxis()->GetBinCenter(i+1);
    }
    double func = Hist->GetBinContent(i+1);
    sum += x*func;
  }

  /*
  // Let's run a debug
  // Get the counts to make a pdf
  double integral = 0.0;
  for (int i = 0; i < Hist->GetXaxis()->GetNbins(); ++i) {
    integral += Hist->GetBinContent(i+1);
  }

  // Get the mean
  double mu = Hist->GetMean();
  // Get the standard deviation
  double stddev = sqrt(sum/integral - mu*mu);

  // Can compare to the mean and std dev from the ROOT functions
  if (order == 1 && mu != Hist->GetMean()) { 
    std::cerr << "EEK" << std::endl;
    std::cerr << "calc: " << mu << " mean: " << Hist->GetMean() << std::endl;
    throw;
  } else if (order == 2 && stddev != Hist->GetRMS()) { 
    std::cerr << "EEK" << std::endl;
    std::cerr << "calc: " << stddev << " stddev: " << Hist->GetRMS() << std::endl;
    //throw;
  }
  */

  return sum;
}

void fsianal(std::string filename) {
  SetCustomPalette();
  TFile *file = new TFile(filename.c_str(), "OPEN");
  TTree *tree = (TTree*)file->Get("FlatTree_VARS");
  tree->SetBranchStatus("*", false);

  int mode;
  int PDGLep;
  bool cc;
  float q0;
  float q3;
  float Eav;
  bool isCCinc;
  float Enu;

  // set up the branches
  tree->SetBranchStatus("q0", true);
  tree->SetBranchAddress("q0", &q0);
  tree->SetBranchStatus("q3", true);
  tree->SetBranchAddress("q3", &q3);
  tree->SetBranchStatus("Eav", true);
  tree->SetBranchAddress("Eav", &Eav);

  tree->SetBranchStatus("Enu_true", true);
  tree->SetBranchAddress("Enu_true", &Enu);

  tree->SetBranchStatus("flagCCINC", true);
  tree->SetBranchAddress("flagCCINC", &isCCinc);
  tree->SetBranchStatus("cc", true);
  tree->SetBranchAddress("cc", &cc);
  tree->SetBranchStatus("Mode", true);
  tree->SetBranchAddress("Mode", &mode);
  tree->SetBranchStatus("PDGLep", true);
  tree->SetBranchAddress("PDGLep", &PDGLep);


  tree->SetBranchStatus("fScaleFactor", true);
  double MinScale = tree->GetMinimum("fScaleFactor");
  double MaxScale = tree->GetMaximum("fScaleFactor");
  if (MinScale != MaxScale) {
    std::cerr << "Different scalefactor for different events" << std::endl;
    throw;
  }
  std::cout << "Scalefactor: " << MinScale << std::endl;

  const int nbins = 100;
  TH3D *plot = new TH3D("q0q3Eav", "q0q3Eav", nbins, 0, 2.0, nbins, 0, 2.0, 24, 0, 1.0);
  plot->GetXaxis()->SetTitle("q_{3} (GeV)");
  plot->GetYaxis()->SetTitle("q_{0} (GeV)");
  plot->GetZaxis()->SetTitle("E_{av}/q_{0}");

  TH2D *plot2d = new TH2D("q0q3", "q0q3", nbins, 0, 2.0, nbins, 0, 2.0);
  plot2d->GetXaxis()->SetTitle("q_{3} (GeV)");
  plot2d->GetYaxis()->SetTitle("q_{0} (GeV)");
  plot2d->GetZaxis()->SetTitle("d#sigma/dq_{3}dq_{0} (cm^{2}/nucleon/GeV/GeV)");

  const int nModesTotal = 6;
  TH3D *plotMode[nModesTotal];
  for (int i = 0; i < nModesTotal; ++i) {
    plotMode[i] = new TH3D(Form("q0q3Eav_%i", i), Form("q0q3Eav_%i", i), nbins, 0, 2.0, nbins, 0, 2.0, 24, 0, 1.0);
  }

  int ncc = 0;
  int nPDGLep = 0;
  int nModes = 0;
  int nCCinc = 0;

  int nEvents = tree->GetEntries();
  for (int i = 0; i < nEvents; ++i) {
    tree->GetEntry(i);
    //if (mode < 31 && mode != 16 && Enu > 3.9 && Enu < 4.2) {
    if (mode < 31 && mode != 16) {
      plot->Fill(q3, q0, Eav/q0);
      plot2d->Fill(q3, q0);
      // CCQE
      if (mode == 1) plotMode[0]->Fill(q3, q0, Eav/q0);
      // 2p2h
      else if (mode == 2) plotMode[1]->Fill(q3, q0, Eav/q0);
      // CC1pi+ or 1pi0
      else if (mode == 11 || mode == 12 || mode == 13) plotMode[2]->Fill(q3, q0, Eav/q0);
      // CC multi-pi
      else if (mode == 21) plotMode[3]->Fill(q3, q0, Eav/q0);
      // CC DIS
      else if (mode == 26) plotMode[4]->Fill(q3, q0, Eav/q0);
      else plotMode[5]->Fill(q3, q0, Eav/q0);
    }
    if (cc) ncc++;
    if (mode < 31) nModes++;
    if (PDGLep == 13) nPDGLep++;
    if (isCCinc) nCCinc++;
  }

  plot->Scale(MinScale);
  plot2d->Scale(MinScale);
  for (int i = 0; i < nModesTotal; ++i) {
    plotMode[i]->Scale(MinScale);
  }

  std::string canvname = filename;
  //std::string canvname = "Flat.root";
  canvname = canvname.substr(0, canvname.find(".root"));
  canvname += "_UpTo1GeV_NoIntCut";
  //canvname += "_Flat";
  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetLeftMargin(canv->GetLeftMargin()*1.25);
  canv->cd();
  canv->Print((canvname+".pdf[").c_str());

  TH2D *MeanPlot = new TH2D("MeanPlot", "MeanPlot", plot->GetXaxis()->GetNbins(), plot->GetXaxis()->GetBinLowEdge(1), plot->GetXaxis()->GetBinLowEdge(plot->GetXaxis()->GetNbins()+1), 
      plot->GetYaxis()->GetNbins(), plot->GetYaxis()->GetBinLowEdge(1), plot->GetYaxis()->GetBinLowEdge(plot->GetYaxis()->GetNbins()+1));
  MeanPlot->GetXaxis()->SetTitle(plot->GetXaxis()->GetTitle());
  MeanPlot->GetYaxis()->SetTitle(plot->GetYaxis()->GetTitle());
  MeanPlot->GetZaxis()->SetTitle("E_{av}/q_{0} Mean");
  TH2D *RMSPlot = new TH2D("RMSPlot", "RMSPlot", plot->GetXaxis()->GetNbins(), plot->GetXaxis()->GetBinLowEdge(1), plot->GetXaxis()->GetBinLowEdge(plot->GetXaxis()->GetNbins()+1), 
      plot->GetXaxis()->GetNbins(), plot->GetYaxis()->GetBinLowEdge(1), plot->GetYaxis()->GetBinLowEdge(plot->GetYaxis()->GetNbins()+1));
  RMSPlot->GetXaxis()->SetTitle(plot->GetXaxis()->GetTitle());
  RMSPlot->GetYaxis()->SetTitle(plot->GetYaxis()->GetTitle());
  RMSPlot->GetZaxis()->SetTitle("E_{av}/q_{0} Variance");

  TH2D *SkewPlot = new TH2D("SkewPlot", "SkewPlot", plot->GetXaxis()->GetNbins(), plot->GetXaxis()->GetBinLowEdge(1), plot->GetXaxis()->GetBinLowEdge(plot->GetXaxis()->GetNbins()+1), 
      plot->GetXaxis()->GetNbins(), plot->GetYaxis()->GetBinLowEdge(1), plot->GetYaxis()->GetBinLowEdge(plot->GetYaxis()->GetNbins()+1));
  SkewPlot->GetXaxis()->SetTitle(plot->GetXaxis()->GetTitle());
  SkewPlot->GetYaxis()->SetTitle(plot->GetYaxis()->GetTitle());
  SkewPlot->GetZaxis()->SetTitle("E_{av}/q_{0} Skew");

  TH2D *ZerothPlot = new TH2D("ZerothPlot", "ZerothPlot", plot->GetXaxis()->GetNbins(), plot->GetXaxis()->GetBinLowEdge(1), plot->GetXaxis()->GetBinLowEdge(plot->GetXaxis()->GetNbins()+1), 
      plot->GetXaxis()->GetNbins(), plot->GetYaxis()->GetBinLowEdge(1), plot->GetYaxis()->GetBinLowEdge(plot->GetYaxis()->GetNbins()+1));
  ZerothPlot->GetXaxis()->SetTitle(plot->GetXaxis()->GetTitle());
  ZerothPlot->GetYaxis()->SetTitle(plot->GetYaxis()->GetTitle());
  ZerothPlot->GetZaxis()->SetTitle("E_{av}/q_{0} Zeroth");

  TH2D *FirstPlot = new TH2D("FirstPlot", "FirstPlot", plot->GetXaxis()->GetNbins(), plot->GetXaxis()->GetBinLowEdge(1), plot->GetXaxis()->GetBinLowEdge(plot->GetXaxis()->GetNbins()+1), 
      plot->GetXaxis()->GetNbins(), plot->GetYaxis()->GetBinLowEdge(1), plot->GetYaxis()->GetBinLowEdge(plot->GetYaxis()->GetNbins()+1));
  FirstPlot->GetXaxis()->SetTitle(plot->GetXaxis()->GetTitle());
  FirstPlot->GetYaxis()->SetTitle(plot->GetYaxis()->GetTitle());
  FirstPlot->GetZaxis()->SetTitle("E_{av}/q_{0} First");

  TH2D *SecondPlot = new TH2D("SecondPlot", "SecondPlot", plot->GetXaxis()->GetNbins(), plot->GetXaxis()->GetBinLowEdge(1), plot->GetXaxis()->GetBinLowEdge(plot->GetXaxis()->GetNbins()+1), 
      plot->GetXaxis()->GetNbins(), plot->GetYaxis()->GetBinLowEdge(1), plot->GetYaxis()->GetBinLowEdge(plot->GetYaxis()->GetNbins()+1));
  SecondPlot->GetXaxis()->SetTitle(plot->GetXaxis()->GetTitle());
  SecondPlot->GetYaxis()->SetTitle(plot->GetYaxis()->GetTitle());
  SecondPlot->GetZaxis()->SetTitle("E_{av}/q_{0} Second");

  TH2D *ThirdPlot = new TH2D("ThirdPlot", "ThirdPlot", plot->GetXaxis()->GetNbins(), plot->GetXaxis()->GetBinLowEdge(1), plot->GetXaxis()->GetBinLowEdge(plot->GetXaxis()->GetNbins()+1), 
      plot->GetXaxis()->GetNbins(), plot->GetYaxis()->GetBinLowEdge(1), plot->GetYaxis()->GetBinLowEdge(plot->GetYaxis()->GetNbins()+1));
  ThirdPlot->GetXaxis()->SetTitle(plot->GetXaxis()->GetTitle());
  ThirdPlot->GetYaxis()->SetTitle(plot->GetYaxis()->GetTitle());
  ThirdPlot->GetZaxis()->SetTitle("E_{av}/q_{0} Third");

  TH2D *FourthPlot = new TH2D("FourthPlot", "FourthPlot", plot->GetXaxis()->GetNbins(), plot->GetXaxis()->GetBinLowEdge(1), plot->GetXaxis()->GetBinLowEdge(plot->GetXaxis()->GetNbins()+1), 
      plot->GetXaxis()->GetNbins(), plot->GetYaxis()->GetBinLowEdge(1), plot->GetYaxis()->GetBinLowEdge(plot->GetYaxis()->GetNbins()+1));
  FourthPlot->GetXaxis()->SetTitle(plot->GetXaxis()->GetTitle());
  FourthPlot->GetYaxis()->SetTitle(plot->GetYaxis()->GetTitle());
  FourthPlot->GetZaxis()->SetTitle("E_{av}/q_{0} Fourth");

  TH1D *ProjectedModes[nModesTotal];
  TLegend *leg = new TLegend(0., 0., 1.0, 1.0);
  for (int i = 0; i < nModesTotal; ++i) {
    ProjectedModes[i] = new TH1D(Form("%i",i), Form("%i",i),1,0,1);
    if (i==0) ProjectedModes[i]->SetLineColor(kBlue-7);
    else if (i==1) ProjectedModes[i]->SetLineColor(kCyan-9);
    else if (i==2) ProjectedModes[i]->SetLineColor(kRed-7);
    else if (i==3) ProjectedModes[i]->SetLineColor(kGreen+2);
    else if (i==4) ProjectedModes[i]->SetLineColor(kMagenta-6);
    else if (i==5) ProjectedModes[i]->SetLineColor(kYellow+1);
    ProjectedModes[i]->SetLineWidth(2);
    ProjectedModes[i]->SetLineStyle(kDashed);
    std::string title;
    if (i == 0) title = "CCQE";
    else if (i == 1) title = "CC 2p2h";
    else if (i == 2) title = "CC1#pi^{+/-/0} Res";
    else if (i == 3) title = "CC Multi-#pi";
    else if (i == 4) title = "CC DIS";
    else if (i == 5) title = "CC Other";
    leg->AddEntry(ProjectedModes[i], title.c_str(), "l");
  }

  canv->cd();
  leg->Draw();
  canv->Print((canvname+".pdf").c_str());

  leg->SetX1NDC(0.60);
  leg->SetY1NDC(0.40);
  leg->SetX2NDC(0.90);
  leg->SetY2NDC(0.88);

  // Now have our plot, get the mean and RMS for every q3, q0 bin
  for (int i = 0; i < plot->GetXaxis()->GetNbins(); ++i) {
    for (int j = 0; j < plot->GetYaxis()->GetNbins(); ++j) {
      TH1D* one_d = plot->ProjectionZ(Form("pz_%i_%i", i+1, j+1), i+1, i+1, j+1, j+1);

      //if (one_d->Integral()/plot->Integral()*100. < 0.01) continue;
      if (one_d->Integral() <= 0) continue;

      double q3lo = plot->GetXaxis()->GetBinLowEdge(i+1);
      double q3hi = plot->GetXaxis()->GetBinLowEdge(i+2);
      double q0lo = plot->GetYaxis()->GetBinLowEdge(j+1);
      double q0hi = plot->GetYaxis()->GetBinLowEdge(j+2);

      one_d->SetTitle(Form("E_{av}/q_{0}: q_{0}=%.2f-%.2f, |q_{3}|=%.2f-%.2f", q0lo, q0hi, q3lo, q3hi));

      // Make the mode histograms
      for (int k = 0; k < nModesTotal; ++k) {
        ProjectedModes[k] = plotMode[k]->ProjectionZ(Form("mode%i_pz_%i_%i", k, i+1, j+1), i+1, i+1, j+1, j+1);
        ProjectedModes[k]->SetLineWidth(2);
        ProjectedModes[k]->SetLineStyle(kDashed);
        if (k==0) ProjectedModes[k]->SetLineColor(kBlue-7);
        else if (k==1) ProjectedModes[k]->SetLineColor(kCyan-9);
        else if (k==2) ProjectedModes[k]->SetLineColor(kRed-7);
        else if (k==3) ProjectedModes[k]->SetLineColor(kGreen+2);
        else if (k==4) ProjectedModes[k]->SetLineColor(kMagenta-6);
        else if (k==5) ProjectedModes[k]->SetLineColor(kYellow+1);
      }

      // Make some lines
      double mean = one_d->GetMean();
      double rms = one_d->GetRMS();
      // Actually want the variance here!
      rms = rms*rms;
      double skew = one_d->GetSkewness();

      MeanPlot->SetBinContent(i+1, j+1, mean);
      RMSPlot->SetBinContent(i+1, j+1, rms);
      SkewPlot->SetBinContent(i+1, j+1, skew);

      // Calculate the different order integrals: GetMeMyIntegral = Integral(x^n f(x))
      double zerothcalc = one_d->Integral();
      double firstcalc = GetMeMyIntegral(one_d, 1);
      double secondcalc = GetMeMyIntegral(one_d, 2);
      double thirdcalc = GetMeMyIntegral(one_d, 3);
      double fourthcalc = GetMeMyIntegral(one_d, 4);
      ZerothPlot->SetBinContent(i+1, j+1, zerothcalc);
      FirstPlot->SetBinContent(i+1, j+1, firstcalc);
      SecondPlot->SetBinContent(i+1, j+1, secondcalc);
      ThirdPlot->SetBinContent(i+1, j+1, thirdcalc);
      FourthPlot->SetBinContent(i+1, j+1, fourthcalc);

      TLegend *leg2 = new TLegend(0.15, 0.6, 0.4, 0.9);
      leg2->AddEntry("", Form("#mu=%.2f", mean), "");
      leg2->AddEntry("", Form("#sigma^{2}=%.2f", rms), "");
      leg2->AddEntry("", Form("Skew=%.2f", skew), "");

      canv->Clear();
      TH1D *fill = new TH1D("fill", "fill", 500, 0, 1.0);
      for (int k = 0; k < fill->GetNbinsX()+1; ++k) {
        if (fill->GetBinCenter(k+1) >= mean-sqrt(rms) && fill->GetBinCenter(k+1) <= mean+sqrt(rms)) fill->SetBinContent(k+1, one_d->GetMaximum());
      }
      fill->SetFillStyle(3002);
      fill->SetFillColor(kGray);
      fill->SetLineColor(kRed);
      one_d->Draw();
      fill->Draw("same");
      one_d->Draw("same");
      for (int k = 0; k < nModesTotal; ++k) {
        if (ProjectedModes[k]->Integral()>0) ProjectedModes[k]->Draw("same");
      }
      leg->Draw("same");
      leg2->Draw("same");
      if (one_d->Integral() == 0) {
        delete fill;
        delete leg2;
        continue;
      }
      //canv->Print((canvname+".pdf").c_str());
      delete fill;
      delete leg2;
    }
  }
  /*
  std::cout << "Found " << ncc << " with cc flag" << std::endl;
  std::cout << "Found " << nModes << " with mode flag" << std::endl;
  std::cout << "Found " << nPDGLep << " with PDGLep flag" << std::endl;
  std::cout << "Found " << nCCinc << " with nCCinc flag" << std::endl;
  std::cout << "Total: " << nEvents << std::endl;
  */

  gStyle->SetPalette(55);
  // Now do the Mean and RMS
  canv->Clear();
  MeanPlot->SetMaximum(1.0);
  MeanPlot->SetMinimum(0);
  MeanPlot->Draw("colz");
  canv->SetRightMargin(canv->GetRightMargin()*1.7);
  canv->Print((canvname+".pdf").c_str());

  canv->Clear();
  RMSPlot->SetMaximum(0.10);
  RMSPlot->SetMinimum(0);
  RMSPlot->Draw("colz");
  canv->Print((canvname+".pdf").c_str());

  gStyle->SetPalette(255, MyPalette);
  canv->Clear();
  SkewPlot->SetMaximum(2);
  SkewPlot->SetMinimum(-2);
  SkewPlot->Draw("colz");
  canv->Print((canvname+".pdf").c_str());

  gStyle->SetPalette(kDarkRainBow);

  canv->Clear();
  ZerothPlot->Draw("colz");
  canv->Print((canvname+".pdf").c_str());

  canv->Clear();
  FirstPlot->Draw("colz");
  canv->Print((canvname+".pdf").c_str());

  canv->Clear();
  SecondPlot->Draw("colz");
  canv->Print((canvname+".pdf").c_str());

  canv->Clear();
  ThirdPlot->Draw("colz");
  canv->Print((canvname+".pdf").c_str());

  canv->Clear();
  FourthPlot->Draw("colz");
  canv->Print((canvname+".pdf").c_str());

  TFile *OutputFile = new TFile((canvname+"_Mean_RMS_Skew.root").c_str(), "RECREATE");
  plot->Write();
  for (int i = 0; i < nModesTotal; ++i) plotMode[i]->Write();
  plot2d->Write();
  MeanPlot->Write();
  RMSPlot->Write();
  SkewPlot->Write();
  ZerothPlot->Write();
  FirstPlot->Write();
  SecondPlot->Write();
  ThirdPlot->Write();
  FourthPlot->Write();

  std::cout << "Finished and wrote to " << OutputFile->GetName() << " and " << canvname << ".pdf" << std::endl;
  OutputFile->Close();

  canv->Print((canvname+".pdf]").c_str());
}
