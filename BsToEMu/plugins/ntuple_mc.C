#define ntuple_mc_cxx
#include "ntuple_mc.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ntuple_mc::Loop(){
  TH1F *h1_delPhi1 = new TH1F("h1_delPhi1", "#Delta#phi of same sign charges", 30, 0, 3.2);
  TH1F *h1_delPhi2 = new TH1F("h1_delPhi2", "#Delta#phi of opposite sign charges", 30, 0, 3.2);
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if( ele_pt -> size() > 1 ) cout <<  mu_pt -> size() << '\t' << ele_pt -> size() << "There is something to change...\n";
    
    for( int i1 = 0; i1 < mu_pt -> size(); i1++ ){
      float delPhi = fabs(mu_pt -> at(i1) - ele_pt -> at(i1));
      while(delPhi > M_PI) delPhi = 2*M_PI - delPhi;
      
      if( mu_charge -> at(i1) * ele_charge->at(i1) > 0) h1_delPhi1 -> Fill(delPhi);
      else if( mu_charge -> at(i1) * ele_charge->at(i1) < 0 ) h1_delPhi2 -> Fill(delPhi);
      else cout << "Something else need to be tested...\n";
    }
  }
  gStyle -> SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1", "c1", 650, 700);
  auto leg = new TLegend(.13,.62,.57,.83);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->AddEntry(h1_delPhi1,"Same sign");
  leg->AddEntry(h1_delPhi2,"Opposite sign","L");

  h1_delPhi1 -> SetLineWidth(2);
  h1_delPhi2 -> SetLineWidth(2);
  h1_delPhi2 -> SetLineColor(kRed);
  
  h1_delPhi2 -> Draw();
  h1_delPhi1 -> Draw("same");
  cout << h1_delPhi1 -> Integral() << '\t' << h1_delPhi2 -> Integral() << endl;
  leg -> Draw();
}
void ntuple_mc::Loop2(){
  TFile *fscale = new TFile("../../june/2/scale_factors_fullBpark.root");
  TH2 *hist_scale_factor = (TH2*)fscale->Get("hist_scale_factor");
  //auto f00 = TFile::Open("Bs2EMu_mc_blind_5p.root", "RECREATE");
  auto f00 = TFile::Open("Bs2EMu_mc_minimal.root", "RECREATE");
  TTree *ntuple = new TTree("ntuple", "MC" );
  float mass_, mcorr_, pt_, eta_, phi_, dr_mu_, dr_ele_, vtx_chi2_, vtx_prob_, cos2d_, lxy_, lxy_err_, lxy_sig_, mu_pt_, mu_eta_, mu_phi_, mu_e_, mu_trkIsolation_, mu_id_soft_mva_, mu_dxy_sig_, ele_pt_, ele_eta_, ele_phi_, ele_e_, delEta_, delPhi_, ele_trkIsolation_, pt1byM_, pt2byM_, triggerSF_, charge_mult_, avg_pt2_;
  //char test_wp[10];
  ntuple -> Branch( "mass", &mass_, "mass/F");
  ntuple -> Branch( "mcorr", &mcorr_, "mcorr/F");
  ntuple -> Branch( "pt", &pt_, "pt/F");
  ntuple -> Branch( "eta", &eta_, "eta/F");
  ntuple -> Branch( "phi", &phi_, "phi/F");
  ntuple -> Branch( "dr_mu", &dr_mu_, "dr_mu/F");
  ntuple -> Branch( "dr_ele", &dr_ele_, "dr_ele/F");
  ntuple -> Branch( "vtx_chi2", &vtx_chi2_, "vtx_chi2/F");
  ntuple -> Branch( "vtx_prob", &vtx_prob_, "vtx_prob/F");
  ntuple -> Branch( "cos2d", &cos2d_, "cos2d/F");
  ntuple -> Branch( "lxy", &lxy_, "lxy/F");
  ntuple -> Branch( "lxy_err", &lxy_err_, "lxy_err/F");
  ntuple -> Branch( "lxy_sig", &lxy_sig_, "lxy_sig/F");
  ntuple -> Branch( "mu_pt", &mu_pt_, "mu_pt/F");
  ntuple -> Branch( "mu_eta", &mu_eta_, "mu_eta/F");
  ntuple -> Branch( "mu_phi", &mu_phi_, "mu_phi/F");
  ntuple -> Branch( "mu_e", &mu_e_, "mu_e/F");
  ntuple -> Branch( "mu_trkIsolation", &mu_trkIsolation_, "mu_trkIsolation/F");
  ntuple -> Branch( "mu_id_soft_mva", &mu_id_soft_mva_, "mu_id_soft_mva/F");
  ntuple -> Branch( "mu_dxy_sig", &mu_dxy_sig_, "mu_dxy_sig/F");
  ntuple -> Branch( "ele_pt", &ele_pt_, "ele_pt/F");
  ntuple -> Branch( "ele_eta", &ele_eta_, "ele_eta/F");
  ntuple -> Branch( "ele_phi", &ele_phi_, "ele_phi/F");
  ntuple -> Branch( "ele_e", &ele_e_, "ele_e/F");
  ntuple -> Branch( "delEta", &delEta_, "delEta/F");
  ntuple -> Branch( "delPhi", &delPhi_, "delPhi/F");
  ntuple -> Branch( "ele_trkIsolation", &ele_trkIsolation_, "ele_trkIsolation/F");
  ntuple -> Branch( "pt1byM", &pt1byM_, "pt1byM/F");
  ntuple -> Branch( "pt2byM", &pt2byM_, "pt2byM/F");
  ntuple -> Branch( "triggerSF", &triggerSF_, "triggerSF/F");
  ntuple -> Branch( "charge_mult", &charge_mult_, "charge_mult/F");
  ntuple -> Branch( "avg_pt2", &avg_pt2_, "avg_pt2/F");
  int xx=0, xy=0;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // There are 62 events in data where 2 sets satisfying same critera in one even.
    // Also in all those, pT of muon is same. I am selecting the one with high vtx probability
        int i1, i2 = -999; double vtxprob = -999;
    for(i1 = 0; i1 < vtx_prob -> size(); i1++){
      if( vtxprob < vtx_prob -> at(i1) ) vtxprob = vtx_prob -> at(i1);
      i2 = i1;
    }
    if( vtxprob < 0 ) continue;
    
    mass_ = mass -> at(i2);
    //if( mass_ > 4.8 && mass_ < 5.8 ) continue;
    pt_ = pt -> at(i2);
    eta_ = eta -> at(i2);
    phi_ = phi -> at(i2);
    dr_mu_ = dr_mu -> at(i2);
    dr_ele_ = dr_ele -> at(i2);
    vtx_chi2_ = vtx_chi2 -> at(i2);
    vtx_prob_ = vtx_prob -> at(i2);
    cos2d_ = cos2d -> at(i2);
    lxy_ = lxy -> at(i2);
    lxy_err_ = lxy_err -> at(i2);
    lxy_sig_ = lxy_sig -> at(i2);
    mu_pt_ = mu_pt -> at(i2);
    mu_eta_ = mu_eta -> at(i2);
    mu_phi_ = mu_phi -> at(i2);
    mu_e_ = mu_e -> at(i2);
    mu_trkIsolation_ = mu_trkIsolation -> at(i2);
    mu_id_soft_mva_ = mu_id_soft_mva -> at(i2);
    mu_dxy_sig_ = mu_dxy_sig->at(i2);
    ele_pt_ = ele_pt -> at(i2);
    ele_eta_ = ele_eta -> at(i2);
    ele_phi_ = ele_phi -> at(i2);
    ele_e_ = ele_e -> at(i2);
    ele_trkIsolation_ = ele_trkIsolation -> at(i2);
    delEta_ = fabs(mu_eta->at(i2) - ele_eta->at(i2));
    delPhi_ = fabs(mu_phi -> at(i2) - ele_phi -> at(i2));
    while(delPhi_ > M_PI) delPhi_ = 2*M_PI - delPhi_;
    float pt1 = mu_pt_; float pt2 = ele_pt_;
    if( pt1 < pt2 ){ pt1 = pt2; pt2 = mu_pt_; }
    pt1byM_ = pt1 / mass_;
    pt2byM_ = pt2 / mass_;
    charge_mult_ = 0;
    avg_pt2_ = 0;
    int ci;
    for(ci=0; ci<charge_pt->size(); ci++){
      if( charge_pt ->at(ci) > 3 ){
	charge_mult_++;
	avg_pt2_ += charge_pt->at(ci) * charge_pt->at(ci);
      }
    }
    if(ci) avg_pt2_ /= ci;

    triggerSF_ = hist_scale_factor->GetBinContent(hist_scale_factor->FindBin(mu_pt->at(i2), fabs(mu_dxy_sig_)));
    //triggerSF_ = 1.0;
    //if(jentry % 20 == 5)
    ntuple->Fill();
  }
  ntuple -> Write();
  f00->Write();
  delete f00;
}
void ntuple_mc::Loop3(){
  TH1F *h1s_npv          = new TH1F("h1s_npv",     "# of PV",         50,    0, 50 );
  TH1F *h1s_pt           = new TH1F("h1s_pt",       "Bs pT",          50,    0, 80 );
  TH1F *h1s_eta          = new TH1F("h1s_eta",      "Bs #eta",        50,   -3, 3  );
  //TH1F *h1s_rap          = new TH1F("h1s_rap",      "Bs rapidity",    50,   -3, 3 );
  TH1F *h1s_phi          = new TH1F("h1s_phi",      "Bs #phi",        50, -3.2, 3.2);
  TH1F *h1s_vtx_chi2     = new TH1F("h1s_vtx_chi2", "vtx #chi^{2}",   50,    0, 40 );
  TH1F *h1s_vtx_prob     = new TH1F("h1s_vtx_prob", "vtx prob",       50,    0,  1 );
  TH1F *h1s_cos2d        = new TH1F("h1s_cos2d",    "cos#alpha_{2D}",      50, 0.99, 1);
  //TH1F *h1s_alpha3d      = new TH1F("h1s_alpha3d",  "#alpha_{3D}",         50,   0, 3.2);
  TH1F *h1s_lxy          = new TH1F("h1s_lxy",      "L_{xy}",              50,   0, 1);
  TH1F *h1s_lxy_err      = new TH1F("h1s_lxy_err",  "#sigma(L_{xy})",      50,   0, 1);
  TH1F *h1s_lxy_sig      = new TH1F("h1s_lxy_sig",  "L_{xy}/#sigma(L_{xy})", 100, 0, 50);
  //TH1F *h1s_l3d          = new TH1F("h1s_l3d",      "L_{3D}",              50,   0, 1);
  //TH1F *h1s_l3d_err      = new TH1F("h1s_l3d_err",  "#sigma(L_{3D})",      50, 0, 1);
  //TH1F *h1s_l3d_sig      = new TH1F("h1s_l3d_sig",  "significance(L_{3D})",50, 0, 50);
  TH1F *h1s_mu_pt        = new TH1F("h1s_mu_pt", "pT_{#mu}", 50, 0, 40);
  TH1F *h1s_mu_eta       = new TH1F("h1s_mu_eta", "#eta_{#mu}", 50, -3, 3);
  //TH1F *h1s_mu_rap       = new TH1F("h1s_mu_rap",    "#mu rapidity",    50,   -3, 3 );
  TH1F *h1s_mu_phi       = new TH1F("h1s_mu_phi", "#phi_{#mu}", 50, -3.2, 3.2);
  //TH1F *h1s_mu_dxy       = new TH1F("h1s_mu_dxy", "d_{xy}(#mu)", 50, -1, 1);
  TH1F *h1s_mu_dxy       = new TH1F("h1s_mu_dxy", "d_{xy}(#mu)", 50, -0.05, 0.05);
  TH1F *h1s_mu_dxy_e     = new TH1F("h1s_dxy_e", "#sigma(dxy#mu)", 50, 0, 1);
  TH1F *h1s_mu_dxy_sig   = new TH1F("h1s_dxy_sig", "d_{xy}/#sigma(d_{xy})(#mu)", 50, -50, 50);
  //TH1F *h1s_mu_dz        = new TH1F("h1s_mu_dz", "dz(#mu)", 50, -1, 1);
  TH1F *h1s_mu_dz        = new TH1F("h1s_mu_dz", "dz(#mu)", 50, -0.05, 0.05);
  TH1F *h1s_mu_dz_e      = new TH1F("h1s_mu_dz_e", "#sigma(dz)(#mu)", 50, 0, 1);
  TH1F *h1s_mu_dz_sig    = new TH1F("h1s_mu_dz_sig", "d_{z}/#sigma(d_{z})(#mu)", 50, -50, 50);
  TH1F *h1s_mu_trkIso    = new TH1F("h1s_mu_trkIso", "#mu trk Isolation", 50, 0, 1);
  TH1F *h1s_mu_trkIso2   = new TH1F("h1s_mu_trkIso2", "#mu trk Isolation", 51, 0, 1.02);
  TH1F *h1s_mu_trk_chi2  = new TH1F("h1s_mu_trk_chi2", "#mu trk #chi^2", 50, 0, 150);
  TH1F *h1s_mu_trk_ndof  = new TH1F("h1s_mu_trk_ndof", "#mu trk nDOF", 50, 0, 50);
  TH1F *h1s_ele_pt       = new TH1F("h1s_ele_pt", "pT(e)", 50, 0, 40);
  TH1F *h1s_ele_eta      = new TH1F("h1s_ele_eta", "#eta(e)", 50, -3, 3);
  //TH1F *h1s_ele_rap      = new TH1F("h1s_ele_rap",   "electron rapidity",    50,   -3, 3 );
  TH1F *h1s_ele_phi      = new TH1F("h1s_ele_phi", "#phi(e)", 50, -3.2, 3.2);
  //TH1F *h1s_ele_dxy      = new TH1F("h1s_ele_dxy", "dxy(e)", 50, -1, 1);
  TH1F *h1s_ele_dxy      = new TH1F("h1s_ele_dxy", "dxy(e)", 50, -0.05, 0.05);
  TH1F *h1s_ele_dxy_e    = new TH1F("h1s_ele_dxy_e", "#sigma(dxy)(e)", 50, 0, 0.1);
  TH1F *h1s_ele_dxy_sig  = new TH1F("h1s_ele_dxy_sig", "dxy/#sigma(dxy)(e)", 50, -50, 50);
  //TH1F *h1s_ele_dz       = new TH1F("h1s_ele_dz", "dz(e)", 50, -1, 1);
  TH1F *h1s_ele_dz       = new TH1F("h1s_ele_dz", "dz(e)", 50, -0.05, 0.05);
  TH1F *h1s_ele_dz_e     = new TH1F("h1s_ele_dz_e", "#sigma(dz)(e)", 50, 0, 0.1);
  TH1F *h1s_ele_dz_sig   = new TH1F("h1s_ele_dz_sig", "dz/#sigma(dz)(e)", 50, -50, 50);
  TH1F *h1s_ele_trkIso   = new TH1F("h1s_ele_trkIso", "I_{e} = #frac{#vec{p_{T}}(e)}{p_{T}(e)+#sum_{trk}#vec{p_{T}}(trk)}", 50, 0, 1);
  TH1F *h1s_ele_trkIso2  = new TH1F("h1s_ele_trkIso2", "I_{e} = #frac{#vec{p_{T}}(e)}{p_{T}(e)+#sum_{trk}#vec{p_{T}}(trk)}", 51, 0, 1.02);
  TH1F *h1s_ele_trk_chi2 = new TH1F("h1s_ele_trk_chi2", "e trk #chi^2", 50, 0, 100);
  TH1F *h1s_ele_trk_ndof = new TH1F("h1s_ele_trk_ndof", "e trk nDOF", 45, 0, 45);
  TH1F *h1s_delEta       = new TH1F("h1s_delEta", "#Delta#eta", 50, 0, 2);
  TH1F *h1s_delPhi = new TH1F("h1s_delPhi", "#Delta#phi", 30, 0, 3.2);

  TH2 *dr_s = new TH2F("dr_s", "#DeltaR with B-candidate (data)", 50, 0, 1.2, 50, 0, 1.2);
  //TH2 *h2s_emu_dz    = new TH2F("h2s_emu_dz", "dz_{mu} vs dz_{e}", 50, 0, 1, 50, 0, 1);
  TH2 *h2s_emu_dz    = new TH2F("h2s_emu_dz", "dz_{mu} vs dz_{e}", 50, 0, 0.05, 50, 0, 0.05);
  //TH2 *h2s_cos2d_alpha3d = new TH2F("h2s_cos2d_alpha3d","cos2d vs #alpha_{3D}", 50, 0.99, 1.001,50, 0, 3);
  //TH2 *h2s_cos2d_alpha3d = new TH2F("h2s_cos2d_alpha3d","cos2d vs #alpha_{3D}", 50, 0.996, 1.001,50, 0, 0.3);
  //TH2 *h2s_lxy_l3d       = new TH2F("h2s_lxy_l3d",  "L_{xy} vs L_{3D}", 50, 0, 1, 50, 0, 1);
  //TH2 *h2s_lxy_l3d_sig   = new TH2F("h2s_lxy_l3d_sig",  "sig(L_{xy}) vs sig(L_{3D})", 50, 0, 50, 50, 0, 50);
  TH2 *h2s_lumi_pv       = new TH2F("h2s_lumi_pv", "lumi vs PV", 50, 0, 50, 50, 0, 3500 );
  TH2 *h2s_eta_emu       = new TH2F("h2s_eta_emu", "#eta(e) vs #eta(#mu)", 50, -4, 4, 50, -4, 4 );
  TH2 *h2s_phi_emu       = new TH2F("h2s_phi_emu", "#phi(e) vs #phi(#mu)", 50, -4.2, 4.2, 50, -4.2, 4.2 );

  TH1F *h1o_npv          = new TH1F("h1o_npv",     "# of PV",         50,    0, 50 );
  TH1F *h1o_pt           = new TH1F("h1o_pt",       "Bs pT",          50,    0, 80 );
  TH1F *h1o_eta          = new TH1F("h1o_eta",      "Bs #eta",        50,   -3, 3  );
  //TH1F *h1o_rap          = new TH1F("h1o_rap",      "Bs rapidity",    50,   -3, 3 );
  TH1F *h1o_phi          = new TH1F("h1o_phi",      "Bs #phi",        50, -3.2, 3.2);
  TH1F *h1o_vtx_chi2     = new TH1F("h1o_vtx_chi2", "vtx #chi^{2}",   50,    0, 40 );
  TH1F *h1o_vtx_prob     = new TH1F("h1o_vtx_prob", "vtx prob",       50,    0,  1 );
  TH1F *h1o_cos2d        = new TH1F("h1o_cos2d",    "cos#alpha_{2D}",      50, 0.99, 1);
  //TH1F *h1o_alpha3d      = new TH1F("h1o_alpha3d",  "#alpha_{3D}",         50,   0, 3.2);
  TH1F *h1o_lxy          = new TH1F("h1o_lxy",      "L_{xy}",              50,   0, 1);
  TH1F *h1o_lxy_err      = new TH1F("h1o_lxy_err",  "#sigma(L_{xy})",      50,   0, 1);
  TH1F *h1o_lxy_sig      = new TH1F("h1o_lxy_sig",  "L_{xy}/#sigma(L_{xy})", 100, 0, 50);
  //TH1F *h1o_l3d          = new TH1F("h1o_l3d",      "L_{3D}",              50,   0, 1);
  //TH1F *h1o_l3d_err      = new TH1F("h1o_l3d_err",  "#sigma(L_{3D})",      50, 0, 1);
  //TH1F *h1o_l3d_sig      = new TH1F("h1o_l3d_sig",  "significance(L_{3D})",50, 0, 50);
  TH1F *h1o_mu_pt        = new TH1F("h1o_mu_pt", "pT_{#mu}", 50, 0, 40);
  TH1F *h1o_mu_eta       = new TH1F("h1o_mu_eta", "#eta_{#mu}", 50, -3, 3);
  //TH1F *h1o_mu_rap       = new TH1F("h1o_mu_rap",    "#mu rapidity",    50,   -3, 3 );
  TH1F *h1o_mu_phi       = new TH1F("h1o_mu_phi", "#phi_{#mu}", 50, -3.2, 3.2);
  //TH1F *h1o_mu_dxy       = new TH1F("h1o_mu_dxy", "d_{xy}(#mu)", 50, -1, 1);
  TH1F *h1o_mu_dxy       = new TH1F("h1o_mu_dxy", "d_{xy}(#mu)", 50, -0.05, 0.05);
  TH1F *h1o_mu_dxy_e     = new TH1F("h1o_dxy_e", "#sigma(dxy#mu)", 50, 0, 1);
  TH1F *h1o_mu_dxy_sig   = new TH1F("h1o_dxy_sig", "d_{xy}/#sigma(d_{xy})(#mu)", 50, -50, 50);
  //TH1F *h1o_mu_dz        = new TH1F("h1o_mu_dz", "dz(#mu)", 50, -1, 1);
  TH1F *h1o_mu_dz        = new TH1F("h1o_mu_dz", "dz(#mu)", 50, -0.05, 0.05);
  TH1F *h1o_mu_dz_e      = new TH1F("h1o_mu_dz_e", "#sigma(dz)(#mu)", 50, 0, 1);
  TH1F *h1o_mu_dz_sig    = new TH1F("h1o_mu_dz_sig", "d_{z}/#sigma(d_{z})(#mu)", 50, -50, 50);
  TH1F *h1o_mu_trkIso    = new TH1F("h1o_mu_trkIso", "#mu trk Isolation", 50, 0, 1);
  TH1F *h1o_mu_trkIso2   = new TH1F("h1o_mu_trkIso2", "#mu trk Isolation", 51, 0, 1.02);
  TH1F *h1o_mu_trk_chi2  = new TH1F("h1o_mu_trk_chi2", "#mu trk #chi^2", 50, 0, 150);
  TH1F *h1o_mu_trk_ndof  = new TH1F("h1o_mu_trk_ndof", "#mu trk nDOF", 50, 0, 50);
  TH1F *h1o_ele_pt       = new TH1F("h1o_ele_pt", "pT(e)", 50, 0, 40);
  TH1F *h1o_ele_eta      = new TH1F("h1o_ele_eta", "#eta(e)", 50, -3, 3);
  //TH1F *h1o_ele_rap      = new TH1F("h1o_ele_rap",   "electron rapidity",    50,   -3, 3 );
  TH1F *h1o_ele_phi      = new TH1F("h1o_ele_phi", "#phi(e)", 50, -3.2, 3.2);
  //TH1F *h1o_ele_dxy      = new TH1F("h1o_ele_dxy", "dxy(e)", 50, -1, 1);
  TH1F *h1o_ele_dxy      = new TH1F("h1o_ele_dxy", "dxy(e)", 50, -0.05, 0.05);
  TH1F *h1o_ele_dxy_e    = new TH1F("h1o_ele_dxy_e", "#sigma(dxy)(e)", 50, 0, 0.1);
  TH1F *h1o_ele_dxy_sig  = new TH1F("h1o_ele_dxy_sig", "dxy/#sigma(dxy)(e)", 50, -50, 50);
  //TH1F *h1o_ele_dz       = new TH1F("h1o_ele_dz", "dz(e)", 50, -1, 1);
  TH1F *h1o_ele_dz       = new TH1F("h1o_ele_dz", "dz(e)", 50, -0.05, 0.05);
  TH1F *h1o_ele_dz_e     = new TH1F("h1o_ele_dz_e", "#sigma(dz)(e)", 50, 0, 0.1);
  TH1F *h1o_ele_dz_sig   = new TH1F("h1o_ele_dz_sig", "dz/#sigma(dz)(e)", 50, -50, 50);
  TH1F *h1o_ele_trkIso   = new TH1F("h1o_ele_trkIso", "I_{e} = #frac{#vec{p_{T}}(e)}{p_{T}(e)+#sum_{trk}#vec{p_{T}}(trk)}", 50, 0, 1);
  TH1F *h1o_ele_trkIso2  = new TH1F("h1o_ele_trkIso2", "I_{e} = #frac{#vec{p_{T}}(e)}{p_{T}(e)+#sum_{trk}#vec{p_{T}}(trk)}", 51, 0, 1.02);
  TH1F *h1o_ele_trk_chi2 = new TH1F("h1o_ele_trk_chi2", "e trk #chi^2", 50, 0, 100);
  TH1F *h1o_ele_trk_ndof = new TH1F("h1o_ele_trk_ndof", "e trk nDOF", 45, 0, 45);
  TH1F *h1o_delEta       = new TH1F("h1o_delEta", "#Delta#eta", 50, 0, 2);
  TH1F *h1o_delPhi = new TH1F("h1o_delPhi", "#Delta#phi", 30, 0, 3.2);

  TH2 *dr_o = new TH2F("dr_o", "#DeltaR with B-candidate (data)", 50, 0, 1.2, 50, 0, 1.2);
  //TH2 *h2o_emu_dz    = new TH2F("h2o_emu_dz", "dz_{mu} vs dz_{e}", 50, 0, 1, 50, 0, 1);
  TH2 *h2o_emu_dz    = new TH2F("h2o_emu_dz", "dz_{mu} vs dz_{e}", 50, 0, 0.05, 50, 0, 0.05);
  //TH2 *h2o_cos2d_alpha3d = new TH2F("h2o_cos2d_alpha3d","cos2d vs #alpha_{3D}", 50, 0.99, 1.001,50, 0, 3);
  //TH2 *h2o_cos2d_alpha3d = new TH2F("h2o_cos2d_alpha3d","cos2d vs #alpha_{3D}", 50, 0.996, 1.001,50, 0, 0.3);
  //TH2 *h2o_lxy_l3d       = new TH2F("h2o_lxy_l3d",  "L_{xy} vs L_{3D}", 50, 0, 1, 50, 0, 1);
  //TH2 *h2o_lxy_l3d_sig   = new TH2F("h2o_lxy_l3d_sig",  "sig(L_{xy}) vs sig(L_{3D})", 50, 0, 50, 50, 0, 50);
  TH2 *h2o_lumi_pv       = new TH2F("h2o_lumi_pv", "lumi vs PV", 50, 0, 50, 50, 0, 3500 );
  TH2 *h2o_eta_emu       = new TH2F("h2o_eta_emu", "#eta(e) vs #eta(#mu)", 50, -4, 4, 50, -4, 4 );
  TH2 *h2o_phi_emu       = new TH2F("h2o_phi_emu", "#phi(e) vs #phi(#mu)", 50, -4.2, 4.2, 50, -4.2, 4.2 );

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if( !hlt_path ){ cout << "HLT not fired.\n"; continue;}
    int i1, i2=0;
    if(t_pt -> size() == 0) continue;
    cout << t_pt -> size() << '\t' << mu_pt -> size() << endl;
    double deta = mu_eta->at(0)-t_eta->at(0);
    double dphi = fabs(mu_phi->at(0)-t_phi->at(0));
    if(dphi > M_PI) dphi =  2*M_PI - dphi;/*
    double t_dr = sqrt(deta*deta+dphi*dphi);
    for(i1 = 1; i1 < mu_eta->size(); i1++){
      double deta1 = mu_eta->at(i1)-t_eta->at(0);
      double dphi1 = fabs(mu_phi->at(1)-t_phi->at(0));
      if(dphi1 > M_PI) dphi1 =  2*M_PI - dphi1;
      double t_dr1 = sqrt(deta1*deta1+dphi1*dphi1);
      if(t_dr > t_dr1){
        i2=i1; t_dr=t_dr1;deta=deta1;dphi=dphi1;
        cout << "Something tried to change.\n";
      }
    }
    if(mu_pt -> size() == 0) continue;
    else if( mu_pt -> size() > 1 ) cout << "Need to filter something...\n";
    float delEta = fabs(mu_eta->at(i2) - ele_eta->at(i2));
    float delPhi = fabs(mu_pt -> at(0) - ele_pt -> at(0));
    while(delPhi > M_PI) delPhi = 2*M_PI - delPhi;*/
    
    //if( mu_charge -> at(i2) == ele_charge->at(i2) ) {
    //  h1s_npv -> Fill(npv);
    //  h2s_lumi_pv       -> Fill(npv, lumi);
    
      /*h1s_pt           -> Fill(pt->at(i2));
      h1s_eta          -> Fill(eta->at(i2));
      //h1s_rap          -> Fill(rap->at(i2));
      h1s_phi          -> Fill(phi->at(i2));
      h1s_vtx_chi2     -> Fill(vtx_chi2->at(i2));
      h1s_vtx_prob     -> Fill(vtx_prob->at(i2));
      h1s_cos2d        -> Fill(cos2d->at(i2));
      //h1s_alpha3d      -> Fill(alpha3d->at(i2));
      h1s_lxy          -> Fill(lxy->at(i2));
      h1s_lxy_err      -> Fill(lxy_err->at(i2));
      h1s_lxy_sig      -> Fill(lxy_sig->at(i2));
      //h1s_l3d          -> Fill(l3d->at(i2));
      //h1s_l3d_err      -> Fill(l3d_err->at(i2));
      //h1s_l3d_sig      -> Fill(l3d_sig->at(i2));
      h1s_mu_pt        -> Fill(mu_pt->at(i2));
      h1s_mu_eta       -> Fill(mu_eta->at(i2));
      //h1s_mu_rap       -> Fill(mu_rap->at(i2));
      h1s_mu_phi       -> Fill(mu_phi->at(i2));
      h1s_mu_dxy       -> Fill(mu_dxy->at(i2));
      h1s_mu_dxy_e     -> Fill(mu_dxy_e->at(i2));
      h1s_mu_dxy_sig   -> Fill(mu_dxy_sig->at(i2));
      h1s_mu_dz        -> Fill(mu_dz->at(i2));
      h1s_mu_dz_e      -> Fill(mu_dz_e->at(i2));
      h1s_mu_dz_sig    -> Fill(mu_dz_sig->at(i2));
      h1s_mu_trkIso    -> Fill(mu_trkIsolation->at(i2));
      h1s_mu_trkIso2   -> Fill(mu_trkIsolation->at(i2));
      h1s_mu_trk_chi2  -> Fill(mu_trk_chi2->at(i2));
      h1s_mu_trk_ndof  -> Fill(mu_trk_ndof->at(i2));
      h1s_ele_pt       -> Fill(ele_pt->at(i2));
      h1s_ele_eta      -> Fill(ele_eta->at(i2));
      //h1s_ele_rap      -> Fill(ele_rap->at(i2));
      h1s_ele_phi      -> Fill(ele_phi->at(i2));
      h1s_ele_dxy      -> Fill(ele_dxy->at(i2));
      h1s_ele_dxy_e    -> Fill(ele_dxy_e->at(i2));
      h1s_ele_dxy_sig  -> Fill(ele_dxy_sig->at(i2));
      h1s_ele_dz       -> Fill(ele_dz->at(i2));
      h1s_ele_dz_e     -> Fill(ele_dz_e->at(i2));
      h1s_ele_dz_sig   -> Fill(ele_dz_sig->at(i2));
      h1s_ele_trkIso   -> Fill(ele_trkIsolation->at(i2));
      h1s_ele_trkIso2  -> Fill(ele_trkIsolation->at(i2));
      h1s_ele_trk_chi2 -> Fill(ele_trk_chi2->at(i2));
      h1s_ele_trk_ndof -> Fill(ele_trk_ndof->at(i2));
      h1s_delEta   -> Fill(delEta);
      h1s_delPhi       -> Fill(delPhi);
      dr_s           -> Fill(dr_mu->at(i2), dr_ele->at(i2));
      h2s_emu_dz    -> Fill(mu_dz->at(i2), ele_dz->at(i2));
      //h2s_cos2d_alpha3d -> Fill(cos2d-> at(i2), alpha3d->at(i2));
      //h2s_lxy_l3d       -> Fill(lxy->at(i2), l3d->at(i2));
      //h2s_lxy_l3d_sig   -> Fill(lxy_sig->at(i2), l3d_sig->at(i2));
      h2s_eta_emu       -> Fill(mu_eta->at(i2), ele_eta->at(i2));
      h2s_phi_emu       -> Fill(mu_phi->at(i2), ele_phi->at(i2));*/
      
      //} else if( mu_charge -> at(0) == -ele_charge->at(0) ) {
      //h1o_npv -> Fill(npv);
      //h2o_lumi_pv       -> Fill(npv, lumi);
    
      /*h1o_pt           -> Fill(pt->at(i2));
      h1o_eta          -> Fill(eta->at(i2));
      //h1o_rap          -> Fill(rap->at(i2));
      h1o_phi          -> Fill(phi->at(i2));
      h1o_vtx_chi2     -> Fill(vtx_chi2->at(i2));
      h1o_vtx_prob     -> Fill(vtx_prob->at(i2));
      h1o_cos2d        -> Fill(cos2d->at(i2));
      //h1o_alpha3d      -> Fill(alpha3d->at(i2));
      h1o_lxy          -> Fill(lxy->at(i2));
      h1o_lxy_err      -> Fill(lxy_err->at(i2));
      h1o_lxy_sig      -> Fill(lxy_sig->at(i2));
      //h1o_l3d          -> Fill(l3d->at(i2));
      //h1o_l3d_err      -> Fill(l3d_err->at(i2));
      //h1o_l3d_sig      -> Fill(l3d_sig->at(i2));
      h1o_mu_pt        -> Fill(mu_pt->at(i2));
      h1o_mu_eta       -> Fill(mu_eta->at(i2));
      //h1o_mu_rap       -> Fill(mu_rap->at(i2));
      h1o_mu_phi       -> Fill(mu_phi->at(i2));
      h1o_mu_dxy       -> Fill(mu_dxy->at(i2));
      h1o_mu_dxy_e     -> Fill(mu_dxy_e->at(i2));
      h1o_mu_dxy_sig   -> Fill(mu_dxy_sig->at(i2));
      h1o_mu_dz        -> Fill(mu_dz->at(i2));
      h1o_mu_dz_e      -> Fill(mu_dz_e->at(i2));
      h1o_mu_dz_sig    -> Fill(mu_dz_sig->at(i2));
      h1o_mu_trkIso    -> Fill(mu_trkIsolation->at(i2));
      h1o_mu_trkIso2   -> Fill(mu_trkIsolation->at(i2));
      h1o_mu_trk_chi2  -> Fill(mu_trk_chi2->at(i2));
      h1o_mu_trk_ndof  -> Fill(mu_trk_ndof->at(i2));
      h1o_ele_pt       -> Fill(ele_pt->at(i2));
      h1o_ele_eta      -> Fill(ele_eta->at(i2));
      //h1o_ele_rap      -> Fill(ele_rap->at(i2));
      h1o_ele_phi      -> Fill(ele_phi->at(i2));
      h1o_ele_dxy      -> Fill(ele_dxy->at(i2));
      h1o_ele_dxy_e    -> Fill(ele_dxy_e->at(i2));
      h1o_ele_dxy_sig  -> Fill(ele_dxy_sig->at(i2));
      h1o_ele_dz       -> Fill(ele_dz->at(i2));
      h1o_ele_dz_e     -> Fill(ele_dz_e->at(i2));
      h1o_ele_dz_sig   -> Fill(ele_dz_sig->at(i2));
      h1o_ele_trkIso   -> Fill(ele_trkIsolation->at(i2));
      h1o_ele_trkIso2  -> Fill(ele_trkIsolation->at(i2));
      h1o_ele_trk_chi2 -> Fill(ele_trk_chi2->at(i2));
      h1o_ele_trk_ndof -> Fill(ele_trk_ndof->at(i2));
      h1o_delEta   -> Fill(delEta);
      h1o_delPhi       -> Fill(delPhi);
      dr_o           -> Fill(dr_mu->at(i2), dr_ele->at(i2));
      h2o_emu_dz    -> Fill(mu_dz->at(i2), ele_dz->at(i2));
      //h2o_cos2d_alpha3d -> Fill(cos2d-> at(i2), alpha3d->at(i2));
      //h2o_lxy_l3d       -> Fill(lxy->at(i2), l3d->at(i2));
      //h2o_lxy_l3d_sig   -> Fill(lxy_sig->at(i2), l3d_sig->at(i2));
      h2o_eta_emu       -> Fill(mu_eta->at(i2), ele_eta->at(i2));
      h2o_phi_emu       -> Fill(mu_phi->at(i2), ele_phi->at(i2));*/

      //} else cout << "Something else need to be tested...\n";

  }
  /*auto f00 = TFile::Open("Bs2EMu_data_hist_gm.root", "RECREATE");

  h1s_npv          -> Write();
  h1s_pt           -> Write();
  h1s_eta          -> Write();
  //h1s_rap          -> Write();
  h1s_phi          -> Write();
  h1s_vtx_chi2     -> Write();
  h1s_vtx_prob     -> Write();
  h1s_cos2d        -> Write();
  //h1s_alpha3d      -> Write();
  h1s_lxy          -> Write();
  h1s_lxy_err      -> Write();
  h1s_lxy_sig      -> Write();
  //h1s_l3d          -> Write();
  //h1s_l3d_err      -> Write();
  //h1s_l3d_sig      -> Write();
  h1s_mu_pt        -> Write();
  h1s_mu_eta       -> Write();
  //h1s_mu_rap       -> Write();
  h1s_mu_phi       -> Write();
  h1s_mu_dxy       -> Write();
  h1s_mu_dxy_e     -> Write();
  h1s_mu_dxy_sig   -> Write();
  h1s_mu_dz        -> Write();
  h1s_mu_dz_e      -> Write();
  h1s_mu_dz_sig    -> Write();
  h1s_mu_trkIso    -> Write();
  h1s_mu_trkIso2   -> Write();
  h1s_mu_trk_chi2  -> Write();
  h1s_mu_trk_ndof  -> Write();
  h1s_ele_pt       -> Write();
  h1s_ele_eta      -> Write();
  //h1s_ele_rap      -> Write();
  h1s_ele_phi      -> Write();
  h1s_ele_dxy      -> Write();
  h1s_ele_dxy_e    -> Write();
  h1s_ele_dxy_sig  -> Write();
  h1s_ele_dz       -> Write();
  h1s_ele_dz_e     -> Write();
  h1s_ele_dz_sig   -> Write();
  h1s_ele_trkIso   -> Write();
  h1s_ele_trkIso2  -> Write();
  h1s_ele_trk_chi2 -> Write();
  h1s_ele_trk_ndof -> Write();
  h1s_delEta       -> Write();
  h1s_delPhi       -> Write();
    
  dr_s              -> Write();
  h2s_emu_dz        -> Write();
  //h2s_cos2d_alpha3d -> Write();
  //h2s_lxy_l3d       -> Write();
  //h2s_lxy_l3d_sig   -> Write();
  h2s_lumi_pv       -> Write();
  h2s_eta_emu       -> Write();
  h2s_phi_emu       -> Write();

  h1o_npv          -> Write();
  h1o_pt           -> Write();
  h1o_eta          -> Write();
  //h1o_rap          -> Write();
  h1o_phi          -> Write();
  h1o_vtx_chi2     -> Write();
  h1o_vtx_prob     -> Write();
  h1o_cos2d        -> Write();
  //h1o_alpha3d      -> Write();
  h1o_lxy          -> Write();
  h1o_lxy_err      -> Write();
  h1o_lxy_sig      -> Write();
  //h1o_l3d          -> Write();
  //h1o_l3d_err      -> Write();
  //h1o_l3d_sig      -> Write();
  h1o_mu_pt        -> Write();
  h1o_mu_eta       -> Write();
  //h1o_mu_rap       -> Write();
  h1o_mu_phi       -> Write();
  h1o_mu_dxy       -> Write();
  h1o_mu_dxy_e     -> Write();
  h1o_mu_dxy_sig   -> Write();
  h1o_mu_dz        -> Write();
  h1o_mu_dz_e      -> Write();
  h1o_mu_dz_sig    -> Write();
  h1o_mu_trkIso    -> Write();
  h1o_mu_trkIso2   -> Write();
  h1o_mu_trk_chi2  -> Write();
  h1o_mu_trk_ndof  -> Write();
  h1o_ele_pt       -> Write();
  h1o_ele_eta      -> Write();
  //h1o_ele_rap      -> Write();
  h1o_ele_phi      -> Write();
  h1o_ele_dxy      -> Write();
  h1o_ele_dxy_e    -> Write();
  h1o_ele_dxy_sig  -> Write();
  h1o_ele_dz       -> Write();
  h1o_ele_dz_e     -> Write();
  h1o_ele_dz_sig   -> Write();
  h1o_ele_trkIso   -> Write();
  h1o_ele_trkIso2  -> Write();
  h1o_ele_trk_chi2 -> Write();
  h1o_ele_trk_ndof -> Write();
  h1o_delEta       -> Write();
  h1o_delPhi       -> Write();
    
  dr_o              -> Write();
  h2o_emu_dz        -> Write();
  //h2o_cos2d_alpha3d -> Write();
  //h2o_lxy_l3d       -> Write();
  //h2o_lxy_l3d_sig   -> Write();
  h2o_lumi_pv       -> Write();
  h2o_eta_emu       -> Write();
  h2o_phi_emu       -> Write();*/

}
