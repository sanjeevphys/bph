#define ntuple_data_cxx
#include "ntuple_data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ntuple_data::Loop(){
  TH1F *h1_trg_dr   = new TH1F("h1_trg_dr", "dr(#mu, trg)", 50, 0, 5);
  TH1F *h1_trg_dr2  = new TH1F("h1_trg_dr2", "dr(#mu, trg)", 50, 0, 1);
  TH1F *h1_trg_dr3  = new TH1F("h1_trg_dr3", "dr(#mu, trg)", 50, 0, 0.1);
  TH1F *h1_trg_dr4  = new TH1F("h1_trg_dr4", "dr(#mu, trg)", 50, 0, 0.01);
  TH1F *h1_trg_dpT  = new TH1F("h1_trg_dpT", "#Delta pT/pT", 50, 0, 10 );
  TH1F *h1_trg_dpT2 = new TH1F("h1_trg_dpT2", "#Delta pT/pT", 50, 0, 1 );
  TH1F *h1_trg_dpT3 = new TH1F("h1_trg_dpT3", "#Delta pT/pT", 50, 0, 0.1 );

  TH1F *ha_trg_dr1  = new TH1F("ha_trg_dr1", "dr(#mu, trg)", 50, 0, 5);
  TH1F *ha_trg_dr2  = new TH1F("ha_trg_dr2", "dr(#mu, trg)", 50, 0, 1);
  TH1F *ha_trg_dr3  = new TH1F("ha_trg_dr3", "dr(#mu, trg)", 50, 0, 0.1);
  TH1F *ha_trg_dr4  = new TH1F("ha_trg_dr4", "dr(#mu, trg)", 50, 0, 0.01);
  TH1F *hb_trg_dr1  = new TH1F("hb_trg_dr1", "dr(#mu, trg)", 50, 0, 5);
  TH1F *hb_trg_dr2  = new TH1F("hb_trg_dr2", "dr(#mu, trg)", 50, 0, 1);
  TH1F *hb_trg_dr3  = new TH1F("hb_trg_dr3", "dr(#mu, trg)", 50, 0, 0.1);
  TH1F *hb_trg_dr4  = new TH1F("hb_trg_dr4", "dr(#mu, trg)", 50, 0, 0.01);
  
  TH1F *h1_pt1byM = new TH1F("h1_pt1byM", "pT1/m", 50, 0, 20);
  TH1F *h1_pt2byM = new TH1F("h1_pt2byM", "pT2/m", 50, 0, 10);
  TH1F *h1_delEta = new TH1F("h1_delEta", "#Delta#eta", 50, 0, 2);
  TH1F *h1_delPhi = new TH1F("h1_delPhi", "#Delta#phi", 50, 0, 2);

  TH1F *h1_vtx_prob     = new TH1F("h1_vtx_prob", "vtx prob",       50,    0,  1 );
  TH1F *h1_vtx_prob2    = new TH1F("h1_vtx_prob2", "vtx prob",   50,    0,  0.05 );
  TH1F *h1_lxy_err      = new TH1F("h1_lxy_err",  "#sigma(L_{xy})",      50,   0, 1);
  TH1F *h1_lxy_err2     = new TH1F("h1_lxy_err2",  "#sigma(L_{xy})",    50,   0, 0.1);
  TH1F *h1_lxy_sig      = new TH1F("h1_lxy_sig",  "L_{xy}/#sigma(L_{xy})", 100, 0, 50);
  TH1F *h1_cos2d        = new TH1F("h1_cos2d",    "cos#alpha_{2D}",      50, 0.99, 1);

  TH1F *h1_mu_trkProb  = new TH1F("h1_mu_trkProb", "muon trk prob", 50, 0, 1);
  TH1F *h1_mu_trk_ndof = new TH1F("h1_mu_trk_ndof", "mu trk nDOF", 45, 0, 45);
  TH1F *h1_mu_soft_mva = new TH1F("h1_mu_soft_mva", "mva muon", 50, -1.01, 1.01 );
  TH1F *h1_mu_dxy_sig  = new TH1F("h1_mu_dxy_sig", "dxy/#sigma(dxy)(#mu)", 50, -50, 50);
  TH1F *h1_mu_dxy_sig2 = new TH1F("h1_mu_dxy_sig2", "dxy/#sigma(dxy)(#mu)", 50, -10, 10);

  TH1F *ha_mu_trkProb  = new TH1F("ha_mu_trkProb", "muon trk prob", 50, 0, 1);
  TH1F *ha_mu_trk_ndof = new TH1F("ha_mu_trk_ndof", "mu trk nDOF", 45, 0, 45);
  TH1F *ha_mu_soft_mva = new TH1F("ha_mu_soft_mva", "mva muon", 50, -1.01, 1.01 );
  TH1F *ha_mu_dxy_sig  = new TH1F("ha_mu_dxy_sig", "dxy/#sigma(dxy)(#mu)", 50, -50, 50);
  TH1F *ha_mu_dxy_sig2 = new TH1F("ha_mu_dxy_sig2", "dxy/#sigma(dxy)(#mu)", 50, -10, 10);

  TH1F *h1_ele_trkProb  = new TH1F("h1_ele_trkProb", "Electron trk prob", 50, 0, 1);
  TH1F *h1_ele_trk_ndof = new TH1F("h1_ele_trk_ndof", "e trk nDOF", 45, 0, 45);
  TH1F *h1_ele_id_iso   = new TH1F("h1_ele_id_iso", "mva iso electron", 50, -1.01, 1.01);
  TH1F *h1_ele_id_noIso = new TH1F("h1_ele_id_noIso", "mva noIso electron", 50, -1.01, 1.01 ) ;
  TH1F *h1_ele_dxy_sig  = new TH1F("h1_ele_dxy_sig", "dxy/#sigma(dxy)(e)", 50, -50, 50);
  TH1F *h1_ele_dxy_sig2 = new TH1F("h1_ele_dxy_sig2", "dxy/#sigma(dxy)(e)", 50, -10, 10);
  TH1F *h1_ele_sigmaietaieta = new TH1F("h1_ele_sigmaietaieta", "#sigma_{i#eta i#eta}", 50, 0, 0.08 );

  TH1F *ha_ele_trkProb  = new TH1F("ha_ele_trkProb", "Electron trk prob", 50, 0, 1);
  TH1F *ha_ele_trk_ndof = new TH1F("ha_ele_trk_ndof", "e trk nDOF", 45, 0, 45);
  TH1F *ha_ele_id_iso   = new TH1F("ha_ele_id_iso", "mva iso electron", 50, -1.01, 1.01);
  TH1F *ha_ele_id_noIso = new TH1F("ha_ele_id_noIso", "mva noIso electron", 50, -1.01, 1.01 ) ;
  TH1F *ha_ele_dxy_sig  = new TH1F("ha_ele_dxy_sig", "dxy/#sigma(dxy)(e)", 50, -50, 50);
  TH1F *ha_ele_dxy_sig2 = new TH1F("ha_ele_dxy_sig2", "dxy/#sigma(dxy)(e)", 50, -10, 10);
  TH1F *ha_ele_sigmaietaieta = new TH1F("ha_ele_sigmaietaieta", "#sigma_{i#eta i#eta}", 50, 0, 0.08 );

  TH2F *h2_pT1vsMemu  = new TH2F("h2_pT1vsMemu", "pT1 vs M_{e#mu}", 50, 4, 7, 50, 0, 50);
  TH2F *h2_pT2vsMemu  = new TH2F("h2_pT2vsMemu", "pT2 vs M_{e#mu}", 50, 4, 7, 50, 0, 40);
  TH2F *h2_dEtavsMemu = new TH2F("h2_dEtavsMemu", "#Delta#eta vs M_{e#mu}", 50, 4, 7, 50, 0, 2 );
  TH2F *h2_dPhivsMemu = new TH2F("h2_dPhivsMemu", "#Delta#phi vs M_{e#mu}", 50, 4, 7, 50, 0, 2 );
  
  auto f00 = TFile::Open("Bs2EMu_data_hist_gm_v2.root", "RECREATE");
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if( !hlt_path ){ cout << "HLT not fired.\n"; continue;}
    int i1, i2=0;
    TLorentzVector mu4, trg4, dif;
    mu4.SetPtEtaPhiM(mu_pt->at(0), mu_eta->at(0), mu_phi->at(0), 0.0005109989461);
    trg4.SetPtEtaPhiM(t_pt->at(0), t_eta->at(0), t_phi->at(0), 0.0005109989461);
    dif = mu4 - trg4;
    //cout <<  dif.Px() << '\t' << dif.Py() << '\t' << dif.Pt() << '\t' << mu_pt->at(0) - t_pt->at(0) << endl;
    double dpT = fabs((dif.Pt())/mu_pt->at(0));
    //double dpT = fabs((mu_pt->at(0)-t_pt->at(0))/mu_pt->at(0));
    double detax = mu_eta->at(0)-t_eta->at(0);
    double dphix = fabs(mu_phi->at(0)-t_phi->at(0));
    if(dphix > M_PI) dphix =  2*M_PI - dphix;
    double t_drx = sqrt(detax*detax+dphix*dphix);
    for(i1 = 1; i1 < mu_eta->size(); i1++){
      mu4.SetPtEtaPhiM(mu_pt->at(i1), mu_eta->at(i1), mu_phi->at(i1), 0.0005109989461);
      dif = mu4 - trg4;
      double dpT1 = fabs((dif.Pt())/mu_pt->at(i1));
      //double dpT1 = fabs((mu_pt->at(i1)-t_pt->at(0))/mu_pt->at(i1));
      if(dpT > dpT1){
	dpT = dpT1;
	cout << "DeltapT/pT can also be applied independent of dR cut.\n";
	double deta1 = mu_eta->at(i1)-t_eta->at(0);
	double dphi1 = fabs(mu_phi->at(1)-t_phi->at(0));
	if(dphi1 > M_PI) dphi1 =  2*M_PI - dphi1;
	t_drx = sqrt(deta1*deta1+dphi1*dphi1);
      }
    }
    if(t_drx == 5 ) continue;
    if( dpT < 0.05 ){
      ha_trg_dr1  -> Fill(t_drx);
      ha_trg_dr2  -> Fill(t_drx);
      ha_trg_dr3  -> Fill(t_drx);
      ha_trg_dr4  -> Fill(t_drx);
    } else {
      hb_trg_dr1  -> Fill(t_drx);
      hb_trg_dr2  -> Fill(t_drx);
      hb_trg_dr3  -> Fill(t_drx);
      hb_trg_dr4  -> Fill(t_drx);
    }
    double deta = mu_eta->at(0)-t_eta->at(0);
    double dphi = fabs(mu_phi->at(0)-t_phi->at(0));
    if(dphi > M_PI) dphi =  2*M_PI - dphi;
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

    h1_trg_dr -> Fill(t_dr);
    h1_trg_dr2 -> Fill(t_dr);
    h1_trg_dr3 -> Fill(t_dr);
    h1_trg_dr4 -> Fill(t_dr);
    h1_trg_dpT -> Fill(dpT);
    h1_trg_dpT2 -> Fill(dpT);
    h1_trg_dpT3 -> Fill(dpT);

    //if(mu_trk_ndof -> at(i2) < 8 ) continue;
    //if(mu_id_soft_mva -> at(i2) < 0) continue;
    
    double pt1 = mu_pt -> at(i2);
    double pt2 = ele_pt -> at(i2);
    if( pt1 < pt2 ) { double temp_pt = pt1; pt1 = pt2; pt2 = temp_pt; }
    h1_pt1byM -> Fill( pt1 / mass -> at(i2) );
    h1_pt2byM -> Fill( pt2 / mass -> at(i2) );
    float delEta = fabs(mu_eta->at(i2) - ele_eta->at(i2));
    h1_delEta   -> Fill(delEta);
    float delPhi = fabs(mu_phi -> at(i2) - ele_phi -> at(i2));
    while(delPhi > M_PI) delPhi = 2*M_PI - delPhi;
    h1_delPhi   -> Fill(delPhi);
    
    h1_vtx_prob      -> Fill( vtx_prob -> at(i2) );
    h1_vtx_prob2     -> Fill( vtx_prob -> at(i2) );
    h1_lxy_err       -> Fill( lxy_err -> at(i2) );
    h1_lxy_err2      -> Fill( lxy_err -> at(i2) );
    h1_lxy_sig       -> Fill( lxy_sig -> at(i2) );
    h1_cos2d         -> Fill( cos2d -> at(i2) );

    h1_mu_trkProb   -> Fill( mu_trk_prob -> at(i2) );
    h1_mu_trk_ndof  -> Fill( mu_trk_ndof -> at(i2) );
    h1_mu_soft_mva  -> Fill( mu_id_soft_mva -> at(i2) );
    h1_mu_dxy_sig   -> Fill( mu_dxy_sig -> at(i2) );
    h1_mu_dxy_sig2  -> Fill( mu_dxy_sig -> at(i2) );
    if( mu_id_loose -> at(i2) ){
      ha_mu_trkProb   -> Fill( mu_trk_prob -> at(i2) );
      ha_mu_trk_ndof  -> Fill( mu_trk_ndof -> at(i2) );
      ha_mu_soft_mva  -> Fill( mu_id_soft_mva -> at(i2) );
      ha_mu_dxy_sig   -> Fill( mu_dxy_sig -> at(i2) );
      ha_mu_dxy_sig2  -> Fill( mu_dxy_sig -> at(i2) );
    }
    h1_ele_trkProb   -> Fill( ele_trk_prob -> at(i2) );
    h1_ele_trk_ndof  -> Fill( ele_trk_ndof -> at(i2) );
    h1_ele_id_iso    -> Fill( ele_id_iso -> at(i2) );
    h1_ele_id_noIso  -> Fill( ele_id_noIso -> at(i2) );
    h1_ele_dxy_sig   -> Fill( ele_dxy_sig -> at(i2) );
    h1_ele_dxy_sig2  -> Fill( ele_dxy_sig -> at(i2) );
    h1_ele_sigmaietaieta  -> Fill( ele_sigmaietaieta -> at(i2) );
    if( ele_id_iso_wpLoose -> at(i2) ){
      ha_ele_trkProb   -> Fill( ele_trk_prob -> at(i2) );
      ha_ele_trk_ndof  -> Fill( ele_trk_ndof -> at(i2) );
      ha_ele_id_iso    -> Fill( ele_id_iso -> at(i2) );
      ha_ele_id_noIso  -> Fill( ele_id_noIso -> at(i2) );
      ha_ele_dxy_sig   -> Fill( ele_dxy_sig -> at(i2) );
      ha_ele_dxy_sig2  -> Fill( ele_dxy_sig -> at(i2) );
      ha_ele_sigmaietaieta  -> Fill( ele_sigmaietaieta -> at(i2) );
    }
    h2_pT1vsMemu   -> Fill(mass -> at(i2), pt1);
    h2_pT2vsMemu   -> Fill(mass -> at(i2), pt2);
    h2_dEtavsMemu  -> Fill(mass -> at(i2), delEta);
    h2_dPhivsMemu  -> Fill(mass -> at(i2), delPhi);

  }
  h2_dPhivsMemu -> Draw("colz");

  h1_trg_dr    -> Write();
  h1_trg_dr2   -> Write();
  h1_trg_dr3   -> Write();
  h1_trg_dr4   -> Write();
  h1_trg_dpT   -> Write();
  h1_trg_dpT2  -> Write();
  h1_trg_dpT3  -> Write();

  ha_trg_dr1   -> Write();
  ha_trg_dr2   -> Write();
  ha_trg_dr3   -> Write();
  ha_trg_dr4   -> Write();
  hb_trg_dr1   -> Write();
  hb_trg_dr2   -> Write();
  hb_trg_dr3   -> Write();
  hb_trg_dr4   -> Write();

  h1_pt1byM  -> Write();
  h1_pt2byM  -> Write();
  h1_delEta  -> Write();
  h1_delPhi  -> Write();

  h1_vtx_prob      -> Write();
  h1_vtx_prob2     -> Write();
  h1_lxy_err       -> Write();
  h1_lxy_err2      -> Write();
  h1_lxy_sig       -> Write();
  h1_cos2d         -> Write();

  h1_mu_trkProb   -> Write();
  h1_mu_trk_ndof  -> Write();
  h1_mu_soft_mva  -> Write();
  h1_mu_dxy_sig   -> Write();
  h1_mu_dxy_sig2  -> Write();
  
  ha_mu_trkProb   -> Write();
  ha_mu_trk_ndof  -> Write();
  ha_mu_soft_mva  -> Write();
  ha_mu_dxy_sig   -> Write();
  ha_mu_dxy_sig2  -> Write();

  h1_ele_trkProb   -> Write();
  h1_ele_trk_ndof  -> Write();
  h1_ele_id_iso    -> Write();
  h1_ele_id_noIso  -> Write();
  h1_ele_dxy_sig   -> Write();
  h1_ele_dxy_sig2  -> Write();
  h1_ele_sigmaietaieta  -> Write();

  ha_ele_trkProb   -> Write();
  ha_ele_trk_ndof  -> Write();
  ha_ele_id_iso    -> Write();
  ha_ele_id_noIso  -> Write();
  ha_ele_dxy_sig   -> Write();
  ha_ele_dxy_sig2  -> Write();
  ha_ele_sigmaietaieta  -> Write();

  h2_pT1vsMemu   -> Write();
  h2_pT2vsMemu   -> Write();
  h2_dEtavsMemu  -> Write();
  h2_dPhivsMemu  -> Write();
}
