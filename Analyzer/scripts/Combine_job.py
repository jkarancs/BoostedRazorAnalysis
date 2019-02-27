import subprocess, os, sys, time
import ROOT as rt
from optparse import OptionParser

# --------------------- Options -------------------------
# Read options from command line
usage = "Usage: python %prog filelists [options]"
parser = OptionParser(usage=usage)
parser.add_option('-M','--method',      dest="method",      type="string",       default="AsymptoticLimits",    help="Model (default=AsymptoticLimits)")
parser.add_option('-d','--datacard',    dest="datacard",    type="string",       default="",                    help='Location of the data card')
parser.add_option('--signif',           dest="signif",      action="store_true", default=False,                 help='Calculate significance')
parser.add_option('--blind',            dest="blind",       action="store_true", default=False,                 help='Do not use data for limits/significance')
(opt,args) = parser.parse_args()

# Additional arguments are the boxes
boxes = args

# Automatically derive these from the data card name
logfile   = str(opt.datacard).replace("/cards/","/combine/").replace(".txt",".log")
directory = str(opt.datacard).split("/cards/")[0]
model     = str(opt.datacard).split("/cards/")[1].split("_")[1].replace("SMS-","")
mg        = str(opt.datacard).split("/cards/")[1].split("_")[2]
mchi      = str(opt.datacard).split("/cards/")[1].split("_")[3]
box       ="_".join(opt.datacard.split("/cards/")[1][:-4].split("_")[4:])

energy = 13
if "HL" in directory:
    energy = 14
elif "HE" in directory:
    energy = 27

# -------------------- Functions ------------------------
icommand=0
def special_call(cmd, verbose=1):
    global icommand, opt
    if verbose:
        print("[%d]" % icommand),
        for i in xrange(len(cmd)): print cmd[i],
        print ""
    ntry = 0
    while True:
        try:
            if subprocess.call(cmd):
                print "ERROR: Problem executing command:"
                print("[%d]" % icommand)
                for i in xrange(len(cmd)): print cmd[i],
                print ""
                print "exiting."
                sys.exit()
        except:
            print "Could not excecute command: "
            print("[%d]" % icommand)
            for i in xrange(len(cmd)): print cmd[i],
            print ""
            print "Wait 10s and continue"
            time.sleep(10)
            ntry += 1
            if ntry == 20: sys.exit()
            continue
        break
    if verbose: print ""
    sys.stdout.flush()
    icommand+=1

# Run command with stdout/stderr saved to logfile
def logged_call(cmd, logfile, append=False):
    global opt
    dirname = os.path.dirname(logfile)
    if dirname != "" and not os.path.exists(dirname):
        special_call(["mkdir", "-p", os.path.dirname(logfile)], 0)
    ntry = 0
    while True:
        try:
            with open(logfile, "a" if append else "w") as log:
                proc = subprocess.Popen(cmd, stdout=log, stderr=log, close_fds=True)
                proc.wait()
        except:
            print "Could not write to disk (IOError), wait 10s and continue"
            time.sleep(10)
            ntry += 1
            if ntry == 20: sys.exit()
            continue
        break

def writeTree(box, model, directory, mg, mchi, xsecULObs, xsecULExpPlus2, xsecULExpPlus, xsecULExp, xsecULExpMinus, xsecULExpMinus2, signif):
    tmpFileName = "%s/%s_xsecUL_mg_%s_mchi_%s_%s.root" %("/tmp", model, mg, mchi, box)
    fileOut = rt.TFile.Open(tmpFileName, "recreate")

    if opt.signif:
        signifTree = rt.TTree("signifTree", "signifTree")
        try:
            from ROOT import MyStruct
        except ImportError:
            myStructCmd = "struct MyStruct{Double_t mg;Double_t mchi; Double_t x; Double_t y; Double_t signif; }"
            rt.gROOT.ProcessLine(myStructCmd)
            from ROOT import MyStruct
        
        s = MyStruct()
        signifTree.Branch("mg", rt.AddressOf(s,"mg"),'mg/D')
        signifTree.Branch("mchi", rt.AddressOf(s,"mchi"),'mchi/D')
        signifTree.Branch("x", rt.AddressOf(s,"x"),'x/D')
        signifTree.Branch("y", rt.AddressOf(s,"y"),'y/D')
        signifTree.Branch("signif_%s"%box, rt.AddressOf(s,"signif"),'signif/D')
    else:
        xsecTree = rt.TTree("xsecTree", "xsecTree")
        try:
            from ROOT import MyStruct
        except ImportError:
            myStructCmd = "struct MyStruct{Double_t mg;Double_t mchi; Double_t x; Double_t y;"
            ixsecUL = 0
            myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+0)
            myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+1)
            myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+2)
            myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+3)
            myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+4)
            myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+5)
            ixsecUL+=6
            myStructCmd += "}"
            rt.gROOT.ProcessLine(myStructCmd)
            from ROOT import MyStruct
        
        s = MyStruct()
        xsecTree.Branch("mg", rt.AddressOf(s,"mg"),'mg/D')
        xsecTree.Branch("mchi", rt.AddressOf(s,"mchi"),'mchi/D')
        xsecTree.Branch("x", rt.AddressOf(s,"x"),'x/D')
        xsecTree.Branch("y", rt.AddressOf(s,"y"),'y/D')    
    
    s.mg = mg
    s.mchi = mchi
    if 'T1x' in model:
        s.x = float(model[model.find('x')+1:model.find('y')].replace('p','.'))
        s.y = float(model[model.find('y')+1:].replace('p','.'))
    elif model == 'T1bbbb':
        s.x = 1
        s.y = 0
    elif model == 'T1tttt':
        s.x = 0
        s.y = 1
    else:
        s.x = -1
        s.y = -1

    if opt.signif:
        s.signif = signif
        
        signifTree.Fill()
        
        fileOut.cd()
        signifTree.Write()
        
        fileOut.Close()
        
        outputFileName = "%s/%s_signif_mg_%s_mchi_%s_%s.root" %(directory, model, mg, mchi, box)
        print "INFO: significance values being written to %s"%outputFileName
        
        special_call(["mv",tmpFileName,outputFileName],0)
    else:
        ixsecUL = 0
        xsecTree.Branch("xsecULObs_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+0)),'xsecUL%i/D'%(ixsecUL+0))
        xsecTree.Branch("xsecULExpPlus2_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+1)),'xsecUL%i/D'%(ixsecUL+1))
        xsecTree.Branch("xsecULExpPlus_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+2)),'xsecUL%i/D'%(ixsecUL+2))
        xsecTree.Branch("xsecULExp_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+3)),'xsecUL%i/D'%(ixsecUL+3))
        xsecTree.Branch("xsecULExpMinus_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+4)),'xsecUL%i/D'%(ixsecUL+4))
        xsecTree.Branch("xsecULExpMinus2_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+5)),'xsecUL%i/D'%(ixsecUL+5))
        exec 's.xsecUL%i = xsecULObs[ixsecUL]'%(ixsecUL+0)
        exec 's.xsecUL%i = xsecULExpPlus2[ixsecUL]'%(ixsecUL+1)
        exec 's.xsecUL%i = xsecULExpPlus[ixsecUL]'%(ixsecUL+2)
        exec 's.xsecUL%i = xsecULExp[ixsecUL]'%(ixsecUL+3)
        exec 's.xsecUL%i = xsecULExpMinus[ixsecUL]'%(ixsecUL+4)
        exec 's.xsecUL%i = xsecULExpMinus2[ixsecUL]'%(ixsecUL+5)
        ixsecUL += 4
        
        xsecTree.Fill()
        
        fileOut.cd()
        xsecTree.Write()
        
        fileOut.Close()
        
        outputFileName = "%s/%s_xsecUL_mg_%s_mchi_%s_%s.root" %(directory, model, mg, mchi, box)
        print "INFO: xsec UL values being written to %s"%outputFileName
        
        special_call(["mv",tmpFileName,outputFileName],0)
        
    return outputFileName
# -------------------------------------------------------

# Combine datacards for boxes if needed
if len(boxes)>1:
    cmd = ["combineCards.py"]
    for i in range(len(boxes)):
        cmd.append(boxes[i]+"="+directory+"/cards/RazorBoost_SMS-"+model+"_"+mg+"_"+mchi+"_"+boxes[i]+".txt")
    logged_call(cmd, opt.datacard)

# Fix data card if needed (i.e. input root file is not in directory+"/cards/")
with open(opt.datacard) as card:
    sed_args = ""
    for line in card:
        if ".root" in line:
            root_file = line.split()[3]
            if os.path.dirname(root_file) != directory+"/cards":
                new_name = directory+"/cards/"+os.path.basename(root_file)
                print "Changing data card input filename:"
                print root_file
                print "--> "+new_name
                sed_args += "s;"+root_file+";"+new_name+";;"
                root_file = new_name
            if not os.path.exists(root_file):
                "\nERROR: input root file for combine does not exists: "+root_file
                sys.exit()
    if len(sed_args):
        special_call(["sed","-i",sed_args,opt.datacard],0)
        print ""

print "INFO: combine output being written to "+logfile
#if opt.blind:
#    if opt.method == "Significance":
#        logged_call(["combine", "-M", opt.method, "-t", "-1", "--expectSignal=1", "-d", opt.datacard], logfile)
#    elif opt.method == "AsymptoticLimits":
#        logged_call(["combine", "-M", opt.method, "--run", "blind", "-d", opt.datacard], logfile)
#else:
#    logged_call(["combine", "-M", opt.method, "-d", opt.datacard], logfile)

#get theory cross sections and errors
refXsec = 0
for line in open('./data/stop'+str(energy)+'TeV.txt' if ('T2' in model) else './data/gluino'+str(energy)+'TeV.txt','r'):
    line = line.replace('\n','')
    if mg==line.split(',')[0]:
        refXsec = float(line.split(',')[1]) #pb

# write results to a tree
xsecULObs       = 0.
xsecULExpMinus2 = 0.
xsecULExpMinus  = 0.
xsecULExp       = 0.
xsecULExpPlus   = 0.
xsecULExpPlus2  = 0.
signif          = 0.
with open(logfile) as log:
    for line in log:
        if "Observed Limit:" in line:
            xsecULObs = refXsec*float(line.split()[4])
        elif "Expected  2.5%:" in line:
            xsecULExpMinus2 = refXsec*float(line.split()[4])
        elif "Expected 16.0%:" in line:
            xsecULExpMinus = refXsec*float(line.split()[4])
        elif "Expected 50.0%:" in line:
            xsecULExp = refXsec*float(line.split()[4])
        elif "Expected 84.0%:" in line:
            xsecULExpPlus = refXsec*float(line.split()[4])
        elif "Expected 97.5%:" in line:
            xsecULExpPlus2 = refXsec*float(line.split()[4])
        elif "Significance:" in line:
            signif = float(line.split()[1])

print ('./data/stop'+str(energy)+'TeV.txt' if ('T2' in model) else './data/gluino'+str(energy)+'TeV.txt')
print xsecULExp
print refXsec

if opt.signif:
    writeTree(box, model, directory, float(mg), float(mchi), [0], [0], [0], [0], [0], [0], signif)
else:
    writeTree(box, model, directory, float(mg), float(mchi), [xsecULObs], [xsecULExpPlus2], [xsecULExpPlus], [xsecULExp], [xsecULExpMinus], [xsecULExpMinus2], 0)
