import ROOT
#file = ROOT.TFile.Open("/data/jkarancs/CMSSW/ntuple/B2GTTreeNtupleExtra_MC_25ns_80X_QCD.root")
file = ROOT.TFile.Open("/data/jkarancs/CMSSW/SusyAnalysis/Ntuples/Validation/CMSSW_8_0_20/src/B2GTTreeNtupleExtra_MC_25ns_80X_QCD.root")
tree = file.Get("B2GTTreeMaker/B2GTree")

def printvars( tree, name, prefix_list ):
    # print variables
    for prefix in prefix_list:
        for branch in tree.GetListOfBranches():
            branchname = branch.GetName()
            if prefix in branchname:
                short = branchname if (prefix_list[0] not in prefix) else branchname[len(prefix):]
                print "    //stream.select(\""+branchname+"\", data."+name+"."+short+");"
    print "    "
    return

# Now print the contents for settings.h
print "  void selectVariables(itreestream& stream, DataStruct& data) {"

printvars( tree, "evt",             ["evt_", "SUSY_"] )
printvars( tree, "met",             ["met_"] )
printvars( tree, "pu",              ["pu_"] )
printvars( tree, "vtx",             ["vtx_"] )
printvars( tree, "syst_scale",      ["scale_"] )
printvars( tree, "syst_pdf",        ["pdf_"] )
printvars( tree, "syst_alphas",     ["alphas_"] )
printvars( tree, "filter",          ["Flag_"] )
printvars( tree, "hlt",             ["HLT_"] )
printvars( tree, "gen",             ["gen_"] )
printvars( tree, "ele",             ["el_"] )
printvars( tree, "mu",              ["mu_"] )
printvars( tree, "jetsAK4Puppi",    ["jetAK4Puppi_"] )
printvars( tree, "jetsAK8Puppi",    ["jetAK8Puppi_"] )
printvars( tree, "subjetsAK8Puppi", ["subjetAK8Puppi_"] )
printvars( tree, "genjetsAK8",      ["genjetAK8SD_"] )

print "  }"
