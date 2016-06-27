import ROOT
import math


#syst_file = "results/Plotter_out_2016_06_17_12h50m48.root"
#syst_file = "results/Plotter_out_2016_06_22_22h41m53.root"
#syst_file = "results/Plotter_out_2016_06_23_11h40m47.root"
#syst_file = "results/Plotter_out_2016_06_23_12h01m42.root" # 76X noPdf systematics separately
syst_file = "results/Plotter_out_2016_06_24_14h28m51.root" # all noPdf syst together

hname_ab = "Counts_vs_R/Syst_vs_R/FailAnyTau32Cut_PassDPhiCut_"
hname_cd = "Counts_vs_R/Syst_vs_R/PassTau32Cuts_PassDPhiCut_"

R_CUT_LOW = 0.2
R_CUT     = 0.4

low_stat = False

f = ROOT.TFile.Open(syst_file)

Backgrounds = [
    "Background",
    #"TTJetsMGHT",
    #"TTJetsMG",
    #"TTJetsNLOFXFX",
    #"TTNLO",
    #"TTNLOHerwig",
    "TTPowheg",
    #"TTPowhegmpiOFF",
    #"TTPowhegnoCR",
    #"TTPowhegHerwig",
    "QCD",
    "ZJets",
    "DYJetsToLL",
    "WJets",
    "TTX",
    "Diboson",
    "Top",
]

Systematics = [
    "Lumi       ",
    "PU         ",
    "Trigger    ",
    "JEC        ",
    "HadTopTagSF",
    "HT         ",
    "AlphaS     ",
    "Scale1     ",
    "Scale2     ",
    "Scale3     ",
]

bin_rmin = -1
bin_rmax = -1

sample_d = []

for bkg in Backgrounds:
    # histos
    ab = f.Get(hname_ab+bkg)
    cd = f.Get(hname_cd+bkg)
    # bins
    acmin = ab.GetXaxis().FindBin(R_CUT_LOW)
    acmax = ab.GetXaxis().FindBin(R_CUT)-(ab.GetXaxis().GetBinLowEdge(ab.GetXaxis().FindBin(R_CUT))-R_CUT<1e-10)
    bdmin = ab.GetXaxis().FindBin(R_CUT)
    bdmax = ab.GetNbinsX()
    Nsyst = ab.GetNbinsY()
    print bkg
    # Loop on systematics
    A0 = 0
    B0 = 0
    C0 = 0
    D0 = 0
    Dpred0 = 0
    d = []
    for nsyst in range(1, Nsyst+1):
        # ABCD counts
        A = ab.Integral(acmin, acmax, nsyst, nsyst)
        B = ab.Integral(bdmin, bdmax, nsyst, nsyst)
        C = cd.Integral(acmin, acmax, nsyst, nsyst)
        D = cd.Integral(bdmin, bdmax, nsyst, nsyst)
        d.append(D)
        # Bin errors
        Aerr = 0
        Berr = 0
        Cerr = 0
        Derr = 0
        for bin in range(acmin, acmax+1):
            Aerr += ab.GetBinError(bin, nsyst)*ab.GetBinError(bin, nsyst)
            Cerr += cd.GetBinError(bin, nsyst)*cd.GetBinError(bin, nsyst)
        for bin in range(bdmin, bdmax+1):
            Berr += ab.GetBinError(bin, nsyst)*ab.GetBinError(bin, nsyst)
            Derr += cd.GetBinError(bin, nsyst)*cd.GetBinError(bin, nsyst)
        Dpred = 0
        if A>0:
            Dpred = B * C / A
            CperAerr = math.sqrt((C*C*Aerr*Aerr + A*A*Cerr*Cerr)/(A*A*A*A))
            Dprederr = math.sqrt(B*B*CperAerr*CperAerr + (C/A)*(C/A)*Berr*Berr)
        if (nsyst == 1):
            A0 = A
            B0 = B
            C0 = C
            D0 = D
            D0pred = Dpred
        if low_stat:
            #print "#syst: "+str(nsyst-1)+" Asyst/A: "+str(A/A0)+" Bsyst/B: "+str(B/B0)+" Csyst/C: "+str(C/C0)+" Dpredsyst/Dpres: "+str(Dpred/D0pred)
            if (nsyst > 1):
                if (Nsyst==21): print Systematics[(nsyst-2)/2],
                print ("(+1 sigma): " if nsyst%2==0 else "(-1 sigma): "),
                if (D0pred != 0):
                    print "{0:+.03f}%".format((Dpred/D0pred-1)*100)
                else:
                    print "NaN (D0pred = 0)"
        else:
            #print "#syst: "+str(nsyst-1)+" Asyst: "+str(A)+" Bsyst: "+str(B)+" Csyst: "+str(C)+" Dsyst: "+str(D)
            #print "#syst: "+str(nsyst-1)+" Asyst/A: "+str(A/A0)+" Bsyst/B: "+str(B/B0)+" Csyst/C: "+str(C/C0)+" Dsyst/D: "+str(D/D0)
            if (nsyst == 1):
                print "#syst: "+str(nsyst-1)+" pred (D): "+str(Dpred)+" +- "+str(Dprederr)+" obs. (D):"+str(D)+" +- "+str(Derr)
            else:
                #print "#syst: "+str(nsyst-1)+" Asyst/A: "+str(A/A0)
                if (Nsyst==21): print Systematics[(nsyst-2)/2],
                print ("(+1 sigma): " if nsyst%2==0 else "(-1 sigma): "),
                if (D0 != 0):
                    print "{0:+.03f}%".format((D/D0-1)*100)
                else:
                    print "(D0=0), count = "+str(D)
    sample_d.append(d)

sum = 0
for i in range(1, len(Backgrounds)): sum += sample_d[i][0]
print
print "D (sum, bkg): "+str(sum)+" = "+str(sample_d[0][0])
