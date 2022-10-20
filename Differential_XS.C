void Differential_XS(){

    TString filename[4];

    filename[0] = "/home/hanjie/Desktop/MicroBooNE/analysis/flat_trees/flat/bnb.ub.num.genie_v3_00_06.flat.root";
    filename[1] = "/home/hanjie/Desktop/MicroBooNE/analysis/flat_trees/flat/bnb.ub.num.genie_v2_12_10.mec.flat.root";
    filename[2] = "/home/hanjie/Desktop/MicroBooNE/analysis/flat_trees/flat/bnb.ub.num.neut_5_4_0_1.flat.root";
    filename[3] = "/home/hanjie/Desktop/MicroBooNE/analysis/flat_trees/flat/bnb.ub.num.nuwro_19_02_1.flat.root";

    TString xs_label[4]={"genie_v3","genie_v2","neut","nuwro"};

    double pi0_bin[12]={0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5};

    TList *pi_l = new TList();
    TH1F *hpi0[4];
    TH1F *hpi0_p[4];
    TH1F *h1pi0_p[4];
    TFile *f1[4];

    for(int iff=0; iff<4; iff++){

       f1[iff] = new TFile(filename[iff]);
       TTree *T = (TTree*) f1[iff]->Get("FlatTree_VARS");

       Int_t nentries = T->GetEntries();

       const int nmax = 100;
       int pdg[nmax];
       float px[nmax], py[nmax], pz[nmax];
       int nfsp;
       Char_t cc;
       float Enu_true;
       double fScal;
       int tgta;
       T->SetBranchAddress("nfsp", &nfsp);
       T->SetBranchAddress("cc", &cc);
       T->SetBranchAddress("tgta", &tgta);
       T->SetBranchAddress("Enu_true", &Enu_true);
       T->SetBranchAddress("pdg", &pdg);
       T->SetBranchAddress("px", &px);
       T->SetBranchAddress("py", &py);
       T->SetBranchAddress("pz", &pz);
       T->SetBranchAddress("fScaleFactor", &fScal);

       Int_t nNC=0;
       Int_t nNC_1=0;

       TString hpi_name = "hpi0_"+xs_label[iff];
       hpi0[iff] = new TH1F(hpi_name,hpi_name,15,0,15);

       TString h_name = "hpi0_p_"+xs_label[iff];
       hpi0_p[iff] = new TH1F(h_name,"pi0 momentum",11, pi0_bin);
       h_name = "h1pi0_p_"+xs_label[iff];
       h1pi0_p[iff] = new TH1F(h_name,"one pi0 momentum",11,pi0_bin);

       for(int i=0; i<nentries; i++){
           T->GetEntry(i);
          
           if(cc==1) continue; // only NC events
           if(Enu_true<0.275) continue;
 
           bool foundpi0=false;
           int npi0=0;
           double maxP=0;
           for(int j=0; j<nfsp; j++){
               if(pdg[j]==111){
                  foundpi0=true;
                  npi0++; 
                  double tmp_p = sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]);
                  if(tmp_p>maxP) maxP=tmp_p;
               }
           }
           if(foundpi0==false) continue; 

           nNC=nNC+1; 
           if(npi0==1){
              nNC_1++;
              h1pi0_p[iff]->Fill(maxP,fScal);
           }
           hpi0_p[iff]->Fill(maxP,fScal);
           hpi0[iff]->Fill(npi0);
       }

       pi_l->Add(hpi0[iff]);
       pi_l->Add(hpi0_p[iff]);
       pi_l->Add(h1pi0_p[iff]);
    }

    TString outfile = "output_differential_xs.root";
    TFile *fo = new TFile(outfile, "RECREATE"); 
    pi_l->Write("pi_hist", TObject::kSingleKey);
    fo->Write();
    fo->Close();

    double xs_data[11]={4.33e-40, 8.88e-40, 9.24e-40, 7.93e-40, 6.65e-40, 5.53e-40, 3.93e-40, 2.28e-40, 1.20e-40, 0.58e-40, 0.27e-40};
    double xs_stat[11]={1.10e-40, 2.03e-40, 1.77e-40, 1.27e-40, 0.97e-40, 0.84e-40, 0.68e-40, 0.52e-40, 0.42e-40, 0.30e-40, 0.19e-40};
    double xs_sys[11]={0.83e-40, 1.59e-40, 1.48e-40, 1.13e-40, 0.90e-40, 0.76e-40, 0.60e-40, 0.45e-40, 0.35e-40, 0.25e-40, 0.15e-40};

    double xs_err[11]={0};
    double bin_c[11]={0};
    for(int i=0; i<11; i++){
        xs_err[i] = sqrt(xs_stat[i]*xs_stat[i]+xs_sys[i]*xs_sys[i]);
        bin_c[i] = pi0_bin[i] + ( pi0_bin[i+1]-pi0_bin[i] )/2.0;
    }

    TGraphErrors *gxs_data = new TGraphErrors(11, bin_c, xs_data, 0, xs_err);

    for(int i=0; i<4; i++){
     for(int j=1; j<12; j++){
         double tmp = hpi0_p[i]->GetBinContent(j);
         double new_con = tmp/(pi0_bin[j+1]-pi0_bin[j]);
         hpi0_p[i]->SetBinContent(j,new_con);

         tmp = h1pi0_p[i]->GetBinContent(j);
         new_con = tmp/(pi0_bin[j+1]-pi0_bin[j]);
         h1pi0_p[i]->SetBinContent(j,new_con);
     }
    }


    int color[4] = {1,4,8,6};

    TCanvas *c1=new TCanvas("c1","c1",1000,1000);
    gxs_data->GetYaxis()->SetRangeUser(0.,2e-39);
    gxs_data->SetMarkerStyle(8);
    gxs_data->SetMarkerColor(2);
    gxs_data->Draw("AP");
    gxs_data->GetXaxis()->SetTitle("pi0 momentum");
    gxs_data->GetYaxis()->SetTitle("xs (cm2/GeV/nucleon)");

    for(int i=0; i<4; i++){
        hpi0_p[i]->Draw("HIST same");
        hpi0_p[i]->SetLineColor(color[i]);
        hpi0_p[i]->SetLineWidth(2);
        h1pi0_p[i]->Draw("HIST same");
        h1pi0_p[i]->SetLineColor(color[i]);
        h1pi0_p[i]->SetLineWidth(2);
        h1pi0_p[i]->SetLineStyle(9);
    }

    TLegend *leg1 = new TLegend(0.5,0.75, 0.85, 0.9);
    leg1->SetNColumns(2);
    leg1->AddEntry(gxs_data,"data", "P");
    for(int ii=0; ii<4; ii++){
       leg1->AddEntry(hpi0_p[ii], xs_label[ii], "L");
       leg1->AddEntry(h1pi0_p[ii], xs_label[ii]+" 1 pi", "L");
    }
    leg1->Draw();

// apply smearing matrix
    TFile *f2 = new TFile("../../4ch_diffPpi0_combined_WCdef_open_rw/output.root");
    TH2D *hsmear = (TH2D*)f2->Get("smear");

    TH1F *hpi0_p_smear = new TH1F("hpi0_p_smear","hpi0_p_smear", 11, pi0_bin);
    TH1F *h1pi0_p_smear = new TH1F("h1pi0_p_smear","h1pi0_p_smear", 11,pi0_bin);
    for(int i=1; i<12; i++){
       double tmpxs1=0;
       double tmpxs2=0;
       for(int j=1; j<12; j++){ 
           double fac = hsmear->GetBinContent(i, j);
           tmpxs1 = tmpxs1 + fac*hpi0_p[0]->GetBinContent(j);
           tmpxs2 = tmpxs2 + fac*h1pi0_p[0]->GetBinContent(j);
       }

       hpi0_p_smear->SetBinContent(i, tmpxs1);   // need to check i should start from 0 or 1
       h1pi0_p_smear->SetBinContent(i, tmpxs2);
    }

    hpi0_p_smear->Draw("HIST");
    hpi0_p[0]->Draw("HIST same");
}
