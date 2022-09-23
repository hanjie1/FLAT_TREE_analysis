void total_XS(){

    //TString filename = "/home/hanjie/Desktop/MicroBooNE/analysis/flat_trees/flat/bnb.ub.num.genie_v3_00_06.flat.root";
    //TString filename = "/home/hanjie/Desktop/MicroBooNE/analysis/flat_trees/flat/bnb.ub.num.genie_v2_12_10.mec.flat.root";
    //TString filename = "/home/hanjie/Desktop/MicroBooNE/analysis/flat_trees/flat/bnb.ub.num.neut_5_4_0_1.flat.root";
    TString filename = "/home/hanjie/Desktop/MicroBooNE/analysis/flat_trees/flat/bnb.ub.num.nuwro_19_02_1.flat.root";

    TFile *f1 = new TFile(filename);
    TTree *T = (TTree*) f1->Get("FlatTree_VARS");

    Int_t nentries = T->GetEntries();

    const int nmax = 100;
    int pdg[nmax];
    int nfsp;
    Char_t cc;
    float Enu_true;
    double fScal;
    T->SetBranchAddress("nfsp", &nfsp);
    T->SetBranchAddress("cc", &cc);
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

    double totalxs = nNC*fScal;
    cout<<"total XS:  "<<totalxs<<endl;

}
