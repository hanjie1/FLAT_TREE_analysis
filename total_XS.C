void total_XS(){

    TString filename[4];

    filename[0] = "/home/hanjie/Desktop/MicroBooNE/analysis/flat_trees/flat/bnb.ub.num.genie_v3_00_06.flat.root";
    filename[1] = "/home/hanjie/Desktop/MicroBooNE/analysis/flat_trees/flat/bnb.ub.num.genie_v2_12_10.mec.flat.root";
    filename[2] = "/home/hanjie/Desktop/MicroBooNE/analysis/flat_trees/flat/bnb.ub.num.neut_5_4_0_1.flat.root";
    filename[3] = "/home/hanjie/Desktop/MicroBooNE/analysis/flat_trees/flat/bnb.ub.num.nuwro_19_02_1.flat.root";

    double total_xs[4]={0};
    int color[4] = {30,40,46,38};
    int ss[4] = {41,45,47,34};
    TGraph *gxs[4];

    for(int iff=0; iff<4; iff++){

       TFile *f1 = new TFile(filename[iff]);
       TTree *T = (TTree*) f1->Get("FlatTree_VARS");

       Int_t nentries = T->GetEntries();

       const int nmax = 100;
       int pdg[nmax];
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
       T->SetBranchAddress("fScaleFactor", &fScal);

       Int_t nNC=0;
       for(int i=0; i<nentries; i++){
           T->GetEntry(i);
          
           if(cc==1) continue; // only NC events
           if(Enu_true<0.275) continue;
 
           bool foundpi0=false;
           for(int j=0; j<nfsp; j++){
               if(pdg[j]==111){
                  foundpi0=true;
                  break;
               }
           }
           if(foundpi0==false) continue; 

           nNC=nNC+1; 
       }

       total_xs[iff] = nNC*fScal*tgta;
       delete T;
       delete f1;

       gxs[iff]=new TGraph();
       double np = iff+1;
       gxs[iff]->SetPoint(0, np, total_xs[iff]);
       gxs[iff]->SetMarkerColor(color[iff]);
       gxs[iff]->SetMarkerStyle(ss[iff]);
       gxs[iff]->SetMarkerSize(2);
    }

    double xs_data = 1.43e-38;
    double xs_data_staterr = 0.21e-38;
    double xs_data_syserr = 0.33e-38;
    double data_err = sqrt(xs_data_staterr*xs_data_staterr+xs_data_syserr*xs_data_syserr);
    double dataxx=2.5;
    
    TGraphErrors *gdata = new TGraphErrors();
    gdata->SetPoint(0, 2.5, xs_data);
    gdata->SetPointError(0, 0.0, data_err);
    gdata->SetMarkerColor(2);
    gdata->SetMarkerStyle(8);
    gdata->SetMarkerSize(2);

    TString xs_label[4]={"genie_v3","genie_v2","neut","nuwro"};

    TMultiGraph *mg = new TMultiGraph();
    for(int ii=0; ii<4; ii++)
       mg->Add(gxs[ii]);
    mg->Add(gdata);
    mg->Draw("AP");

    TLegend *leg = new TLegend(0.7,0.7, 0.85, 0.9);
    for(int ii=0; ii<4; ii++)
       leg->AddEntry(gxs[ii], xs_label[ii], "P");
    leg->AddEntry(gdata,"data", "P");
    leg->Draw();

}
