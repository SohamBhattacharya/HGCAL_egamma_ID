# ifndef Analyze_elMu_H
# define Analyze_elMu_H


namespace CommonUtils
{
    std::vector <double> range(double start, double end, double incr)
    {
        int n = (int) std::ceil((end - start) / incr);
        
        std::vector <double> v = {};
        
        for(int i = 0; i <= n; i++)
        {
            v.push_back(start +  i*incr);
        }
        
        return v;
    }
    
    std::vector <double> v_pTbin = range(0, 2000, 0.25);
    std::vector <double> v_etaBin = range(1, 4, 0.01);
    
    
    std::vector <TTree*> buildEventWeights(
        std::vector <TChain*> v_tree,
        std::vector <double> v_treeWeight,
        std::string outBranchName,
        std::string classStr
    )
    {
        int nTree = v_tree.size();
        
        std::vector <double> v_outWeight(nTree, 0.0);
        std::vector <TTree*> v_outTree;
        
        for(int iTree = 0; iTree < nTree; iTree++)
        {
            TChain *inTree = v_tree.at(iTree);
            
            char name[5000];
            sprintf(name, "tree_eventWeight_%s_%d", classStr.c_str(), iTree+1);
            
            v_outTree.push_back(new TTree(name, name));
            TTree *outTree = v_outTree.back();
            
            outTree->Branch(outBranchName.c_str(), &v_outWeight.at(iTree));
            
            int nEvent = inTree->GetEntries();
            
            for(int iEvent = 0; iEvent < nEvent; iEvent++)
            {
                if(iEvent == 0 || iEvent == nEvent-1 || !((iEvent+1) % 10000))
                {
                    printf("Creating event weights for tree %d/%d: event %d/%d. \n", iTree+1, nTree, iEvent+1, nEvent);
                }
                
                v_outWeight.at(iTree) = v_treeWeight.at(iTree);
                
                outTree->Fill();
            }
            
            
            // Activate the branches
            //outTree->SetBranchStatus("*", 1);
        }
        
        return v_outTree;
    }
    
    
    // Returns new trees (containing the weight branches) for sig and bkg
    std::map <std::string, std::vector <TTree*> > buildPtEtaWeights(
        std::string branchName_pT,
        std::string branchName_eta,
        std::string cut_sig,
        std::string cut_bkg,
        std::vector <TChain*> v_tree_sig,
        std::vector <TChain*> v_tree_bkg,
        std::vector <double> v_treeWeight_sig,
        std::vector <double> v_treeWeight_bkg,
        std::string outBranchName_PtEtaWeight
    )
    {
        int nSig = v_tree_sig.size();
        int nBkg = v_tree_bkg.size();
        
        
        char drawStr[5000];
        sprintf(drawStr, "%s:%s", branchName_pT.c_str(), branchName_eta.c_str());
        
        
        
        // Signal
        char histName_pT_vs_eta_sig[5000];
        sprintf(histName_pT_vs_eta_sig, "pT_vs_eta_sig");
        
        char drawStr_sig[5000];
        sprintf(drawStr_sig, "%s:%s >> +%s", branchName_pT.c_str(), branchName_eta.c_str(), histName_pT_vs_eta_sig);
        
        TH2F *h2_pT_vs_eta_sig = new TH2F(histName_pT_vs_eta_sig, histName_pT_vs_eta_sig, v_etaBin.size()-1, &v_etaBin.at(0), v_pTbin.size()-1, &v_pTbin.at(0));
        h2_pT_vs_eta_sig->Sumw2();
        
        for(int iSig = 0; iSig < nSig; iSig++)
        {
            char weightStr[5000];
            sprintf(weightStr, "(%s) * %e", cut_sig.c_str(), v_treeWeight_sig.at(iSig));
            
            printf(
                "buildPtEtaWeights(...); Drawing: \n"
                "Expression : %s \n"
                "Weight     : %s \n"
                "\n",
                drawStr_sig,
                weightStr
            );
            
            v_tree_sig[iSig]->Draw(drawStr_sig, weightStr, "goff");
        }
        
        printf("%s: integral %f, meanX %f, meanY %f \n", histName_pT_vs_eta_sig, h2_pT_vs_eta_sig->Integral(), h2_pT_vs_eta_sig->GetMean(1), h2_pT_vs_eta_sig->GetMean(2));
        h2_pT_vs_eta_sig->Scale(1.0 / h2_pT_vs_eta_sig->Integral(), "width");
        
        
        printf("\n\n");
        
        
        // Background
        char histName_pT_vs_eta_bkg[5000];
        sprintf(histName_pT_vs_eta_bkg, "pT_vs_eta_bkg");
        
        char drawStr_bkg[5000];
        sprintf(drawStr_bkg, "%s:%s >> +%s", branchName_pT.c_str(), branchName_eta.c_str(), histName_pT_vs_eta_bkg);
        
        TH2F *h2_pT_vs_eta_bkg = new TH2F(histName_pT_vs_eta_bkg, histName_pT_vs_eta_bkg, v_etaBin.size()-1, &v_etaBin.at(0), v_pTbin.size()-1, &v_pTbin.at(0));
        h2_pT_vs_eta_bkg->Sumw2();
        
        for(int iBkg = 0; iBkg < nBkg; iBkg++)
        {
            char weightStr[5000];
            sprintf(weightStr, "(%s) * %e", cut_bkg.c_str(), v_treeWeight_bkg.at(iBkg));
            
            printf(
                "buildPtEtaWeights(...); Drawing: \n"
                "Expression : %s \n"
                "Weight     : %s \n"
                "\n",
                drawStr_bkg,
                weightStr
            );
            
            v_tree_bkg[iBkg]->Draw(drawStr_bkg, weightStr, "goff");
        }
        
        printf("%s: integral %f, meanX %f, meanY %f \n", histName_pT_vs_eta_bkg, h2_pT_vs_eta_bkg->Integral(), h2_pT_vs_eta_bkg->GetMean(1), h2_pT_vs_eta_bkg->GetMean(2));
        h2_pT_vs_eta_bkg->Scale(1.0 / h2_pT_vs_eta_bkg->Integral(), "width");
        
        
        // Weights
        TH2F *h2_weight_pT_vs_eta_sig = new TH2F("weight_pT_vs_eta_sig", "weight_pT_vs_eta_sig", v_etaBin.size()-1, &v_etaBin.at(0), v_pTbin.size()-1, &v_pTbin.at(0));
        h2_weight_pT_vs_eta_sig->Sumw2();
        
        TH2F *h2_weight_pT_vs_eta_bkg = new TH2F("weight_pT_vs_eta_bkg", "weight_pT_vs_eta_bkg", v_etaBin.size()-1, &v_etaBin.at(0), v_pTbin.size()-1, &v_pTbin.at(0));
        h2_weight_pT_vs_eta_bkg->Sumw2();
        
        TH1F *h1_weight_pT_sig = new TH1F("weight_pT_sig", "weight_pT_sig", v_pTbin.size()-1, &v_pTbin.at(0));
        h1_weight_pT_sig->Sumw2();
        
        TH1F *h1_weight_pT_bkg = new TH1F("weight_pT_bkg", "weight_pT_bkg", v_pTbin.size()-1, &v_pTbin.at(0));
        h1_weight_pT_bkg->Sumw2();
        
        int nBinX = h2_weight_pT_vs_eta_sig->GetNbinsX();
        int nBinY = h2_weight_pT_vs_eta_sig->GetNbinsY();
        
        for(int iBinX = 0; iBinX < nBinX; iBinX++)
        {
            for(int iBinY = 0; iBinY < nBinY; iBinY++)
            {
                double binContent_sig = h2_pT_vs_eta_sig->GetBinContent(iBinX+1, iBinY+1);
                double binError_sig = h2_pT_vs_eta_sig->GetBinError(iBinX+1, iBinY+1);
                
                double binContent_bkg = h2_pT_vs_eta_bkg->GetBinContent(iBinX+1, iBinY+1);
                double binError_bkg = h2_pT_vs_eta_bkg->GetBinError(iBinX+1, iBinY+1);
                
                double binContent_weight_sig = 0;
                double binError_weight_sig = 0;
                
                double binContent_weight_bkg = 0;
                double binError_weight_bkg = 0;
                
                if(binContent_sig)
                {
                    binContent_weight_sig = std::min(binContent_sig, binContent_bkg) / binContent_sig;
                }
                
                else
                {
                    binContent_weight_sig = 1;
                }
                
                if(binContent_bkg)
                {
                    binContent_weight_bkg = std::min(binContent_sig, binContent_bkg) / binContent_bkg;
                }
                
                else
                {
                    binContent_weight_bkg = 1;
                }
                
                h2_weight_pT_vs_eta_sig->SetBinContent(iBinX+1, iBinY+1, binContent_weight_sig);
                h2_weight_pT_vs_eta_bkg->SetBinContent(iBinX+1, iBinY+1, binContent_weight_bkg);
            }
        }
        
        std::vector <TTree*> v_tree_weight_sig;
        std::vector <TTree*> v_tree_weight_bkg;
        
        for(int iSig = 0; iSig < nSig; iSig++)
        {
            char name[5000];
            sprintf(name, "tree_pT_eta_weight_sig_%d", iSig+1);
            
            v_tree_weight_sig.push_back(new TTree(name, name));
        }
        
        for(int iBkg = 0; iBkg < nBkg; iBkg++)
        {
            char name[5000];
            sprintf(name, "tree_pT_eta_weight_bkg_%d", iBkg+1);
            
            v_tree_weight_bkg.push_back(new TTree(name, name));
        }
        
        
        TH1F *h1_temp_sig = new TH1F("h1_temp_sig", "h1_temp_sig", v_pTbin.size()-1, &v_pTbin.at(0));
        TH1F *h1_temp_bkg = new TH1F("h1_temp_bkg", "h1_temp_bkg", v_pTbin.size()-1, &v_pTbin.at(0));
        
        
        // Create sig weight branch
        std::vector <std::vector <double>*> vv_weight_pT_vs_eta_sig(nSig, 0);
        
        //# pragma omp parallel for
        for(int iSig = 0; iSig < nSig; iSig++)
        {
            TChain *tree_sig = v_tree_sig.at(iSig);
            TTree *tree_weight_sig = v_tree_weight_sig.at(iSig);
            
            std::vector <double> v_weight_pT_vs_eta_sig;
            vv_weight_pT_vs_eta_sig.at(iSig) = &v_weight_pT_vs_eta_sig;
            tree_weight_sig->Branch(outBranchName_PtEtaWeight.c_str(), &v_weight_pT_vs_eta_sig);
            
            int nEvent_sig = tree_sig->GetEntries();
            
            char drawStr_mod[5000];
            sprintf(drawStr_mod, "Length$(%s)", branchName_pT.c_str());
            
            int nEvent = tree_sig->Draw(drawStr_mod, "", "goff");
            std::vector <double> v_nObj(tree_sig->GetV1(), tree_sig->GetV1() + nEvent);
            
            int countPt = tree_sig->Draw(branchName_pT.c_str(), "", "goff");
            std::vector <double> v_pT(tree_sig->GetV1(), tree_sig->GetV1() + countPt);
            
            int countEta = tree_sig->Draw(branchName_eta.c_str(), "", "goff");
            std::vector <double> v_eta(tree_sig->GetV1(), tree_sig->GetV1() + countEta);
            
            int count = 0;
            
            for(int iEvent = 0; iEvent < nEvent; iEvent++)
            {
                if(iEvent == 0 || iEvent == nEvent-1 || !((iEvent+1) % 10000))
                {
                    printf("Creating pT-eta weights for sig. tree %d/%d: event %d/%d. \n", iSig+1, nSig, iEvent+1, nEvent);
                }
                
                v_weight_pT_vs_eta_sig.clear();
                
                int nObj = v_nObj.at(iEvent);
                
                if(nObj)
                {
                    for(int iObj = 0; iObj < nObj; iObj++)
                    {
                        double pT = v_pT.at(count);
                        double eta = v_eta.at(count);
                        
                        int binX = h2_weight_pT_vs_eta_sig->GetXaxis()->FindBin(eta);
                        int binY = h2_weight_pT_vs_eta_sig->GetYaxis()->FindBin(pT);
                        
                        //printf("Tree %d/%d weight: %0.4e \n", tree_sig->GetTreeNumber()+1, tree_sig->GetNtrees(), tree_sig->GetTree()->GetWeight());
                        double weight = h2_weight_pT_vs_eta_sig->GetBinContent(binX, binY);
                        
                        v_weight_pT_vs_eta_sig.push_back(weight);
                        
                        //printf("sig iObj %d pT %f, eta %f, weight %f \n", iObj, pT, eta, weight);
                        
                        h1_temp_sig->Fill(pT, weight);
                        
                        count++;
                    }
                }
                
                //printf("[%d] v_weight_pT_vs_eta_sig.size() %d \n", iEvent, (int) v_weight_pT_vs_eta_sig.size());
                
                //tree_weight_sig->GetBranch(outBranchName_PtEtaWeight.c_str())->Fill();
                tree_weight_sig->Fill();
            }
            
            
            // Activate the branches
            //tree_weight_sig->SetBranchStatus("*", 1);
        }
        
        
        // Create bkg weight branch
        std::vector <std::vector <double>*> vv_weight_pT_vs_eta_bkg(nBkg, 0);
        
        for(int iBkg = 0; iBkg < nBkg; iBkg++)
        {
            TChain *tree_bkg = v_tree_bkg.at(iBkg);
            TTree *tree_weight_bkg = v_tree_weight_bkg.at(iBkg);
            
            std::vector <double> v_weight_pT_vs_eta_bkg;
            vv_weight_pT_vs_eta_bkg.at(iBkg) = &v_weight_pT_vs_eta_bkg;
            tree_weight_bkg->Branch(outBranchName_PtEtaWeight.c_str(), &v_weight_pT_vs_eta_bkg);
            
            int nEvent_bkg = tree_bkg->GetEntries();
            
            char drawStr_mod[5000];
            sprintf(drawStr_mod, "Length$(%s)", branchName_pT.c_str());
            
            int nEvent = tree_bkg->Draw(drawStr_mod, "", "goff");
            std::vector <double> v_nObj(tree_bkg->GetV1(), tree_bkg->GetV1() + nEvent);
            
            int countPt = tree_bkg->Draw(branchName_pT.c_str(), "", "goff");
            std::vector <double> v_pT(tree_bkg->GetV1(), tree_bkg->GetV1() + countPt);
            
            int countEta = tree_bkg->Draw(branchName_eta.c_str(), "", "goff");
            std::vector <double> v_eta(tree_bkg->GetV1(), tree_bkg->GetV1() + countEta);
            
            int count = 0;
            
            for(int iEvent = 0; iEvent < nEvent; iEvent++)
            {
                if(iEvent == 0 || iEvent == nEvent-1 || !((iEvent+1) % 10000))
                {
                    printf("Creating pT-eta weights for bkg. tree %d/%d: event %d/%d. \n", iBkg+1, nBkg, iEvent+1, nEvent);
                }
                
                v_weight_pT_vs_eta_bkg.clear();
                
                int nObj = v_nObj.at(iEvent);
                
                if(nObj)
                {
                    for(int iObj = 0; iObj < nObj; iObj++)
                    {
                        double pT = v_pT.at(count);
                        double eta = v_eta.at(count);
                        
                        int binX = h2_weight_pT_vs_eta_bkg->GetXaxis()->FindBin(eta);
                        int binY = h2_weight_pT_vs_eta_bkg->GetYaxis()->FindBin(pT);
                        
                        //printf("Tree %d/%d weight: %0.4e \n", tree_bkg->GetTreeNumber()+1, tree_bkg->GetNtrees(), tree_bkg->GetTree()->GetWeight());
                        double weight = h2_weight_pT_vs_eta_bkg->GetBinContent(binX, binY);
                        
                        v_weight_pT_vs_eta_bkg.push_back(weight);
                        
                        //printf("bkg iObj %d pT %f, eta %f, weight %f \n", iObj, pT, eta, weight);
                        
                        h1_temp_bkg->Fill(pT, weight);
                        
                        count++;
                    }
                }
                
                //printf("[%d] v_weight_pT_vs_eta_bkg.size() %d \n", iEvent, (int) v_weight_pT_vs_eta_bkg.size());
                
                //tree_weight_bkg->GetBranch(outBranchName_PtEtaWeight.c_str())->Fill();
                tree_weight_bkg->Fill();
            }
            
            
            // Activate the branches
            //tree_weight_bkg->SetBranchStatus("*", 1);
        }
        
        
        std::map <std::string, std::vector <TTree*> > m_tree_weight;
        
        m_tree_weight["sig"] = v_tree_weight_sig;
        m_tree_weight["bkg"] = v_tree_weight_bkg;
        
        fflush(stdout);
        
        return m_tree_weight;
    }
    
    
    int loadTChain(std::string sourceFileName, TChain *chain, int nFile_max = -1)
    {
        printf("Loading input files from: %s \n", sourceFileName.c_str());
        
        std::ifstream sourceFileList;
        sourceFileList.open(sourceFileName);
        
        std::string inFileName;
        
        int nFile = 0;
        
        while(std::getline(sourceFileList, inFileName))
        {
            //if(inFileName.find("/eos/") == 0)
            //{
            //    inFileName = "root://cms-xrd-global.cern.ch/" + inFileName;
            //}
            
            TFile *inFile = (TFile*) TFile::Open(inFileName.c_str());
            
            if(inFile && !inFile->IsZombie())
            {
                //v_inFile.push_back(inFile);
                
                chain->Add(inFileName.c_str());
                
                inFile->Close();
                
                nFile++;
                //printf("Loaded %d input files. \n", nFile);
                
                if(nFile == nFile_max)
                {
                    break;
                }
            }
        }
        
        printf("Loaded %d input files. \n\n", nFile);
        
        return nFile;
    }
}


# endif
