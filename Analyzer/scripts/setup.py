import os, re, sys, glob, socket, subprocess

# 74X
#LATEST_NTUPLE_EOS="Skim_Feb22_1AK8JetPt300"
#LATEST_NTUPLE_GRID18="Skim_Feb22_1AK8JetPt300"

# 76X
#LATEST_NTUPLE_EOS="Skim_Apr28_1AK8JetPt300"
#LATEST_NTUPLE_EOS="Skim_May21_1AK8JetPt300"
#LATEST_NTUPLE_EOS="Skim_Jul07_SignalRegion"
#LATEST_NTUPLE_EOS="Skim_Jul13_SignalRegion"
#LATEST_NTUPLE_EOS="Skim_Jul19_SignalRegion"
#LATEST_NTUPLE_GRID18="Apr13_edm_Apr01"
#LATEST_NTUPLE_GRID18="Skim_Apr28_1AK8JetPt300"
#LATEST_NTUPLE_GRID18="Jun08"
#LATEST_NTUPLE_GRID18="Skim_May21_1AK8JetPt300"
#LATEST_NTUPLE_GRID18="Skim_Jul19_SignalRegion"

# 80X
#LATEST_NTUPLE_EOS="Skim_Aug30_1AK8JetPt300"
#LATEST_NTUPLE_EOS="Skim_Oct31_2Jet_1JetAK8"
#LATEST_NTUPLE_EOS="Skim_Feb26_1JetAK8_0p04R2"
#LATEST_NTUPLE_EOS="Skim_Mar09_1JetAK8_0p04R2"
#LATEST_NTUPLE_EOS="Skim_May19_1JetAK8_0p04R2"
#LATEST_NTUPLE_EOS="Skim_Sep26"
#LATEST_NTUPLE_EOS="Skim_Oct11"
LATEST_NTUPLE_EOS="Skim_Dec12"
#LATEST_NTUPLE_GRID18="Aug17"
#LATEST_NTUPLE_GRID18="Skim_Aug30_1AK8JetPt300"
#LATEST_NTUPLE_GRID18="Oct24"
#LATEST_NTUPLE_GRID18="Skim_Oct31_2Jet_1JetAK8"
#LATEST_NTUPLE_GRID18="Jan12"
#LATEST_NTUPLE_GRID18="Skim_Feb26_1JetAK8_0p04R2"
#LATEST_NTUPLE_GRID18="Skim_Mar09_1JetAK8_0p04R2"
#LATEST_NTUPLE_GRID18="May10"
#LATEST_NTUPLE_GRID18="May10_part7"
#LATEST_NTUPLE_GRID18="Skim_May19_1JetAK8_0p04R2"
#LATEST_NTUPLE_GRID18="Sep26_part1"
#LATEST_NTUPLE_GRID18="Sep26_part2"
#LATEST_NTUPLE_GRID18="Skim_Oct11"
#LATEST_NTUPLE_GRID18="Nov30_part3"
LATEST_NTUPLE_GRID18="Skim_Dec12"

ANA_BASE = os.environ['CMSSW_BASE']+'/src/BoostedRazorAnalysis/Analyzer'
if 'grid18.kfki.hu' in socket.gethostname(): ANA_BASE='/data/jkarancs/CMSSW/Analyzer'
DIR = ANA_BASE+'/ntuple/Latest'

if 'lxplus' in socket.gethostname():
    print 'Running on lxplus'
    if os.path.lexists(ANA_BASE+'/ntuple/Latest'):
        if os.path.islink(ANA_BASE+'/ntuple/Latest'):
            print 'Removing soft-link directory (used by previous version): '+ANA_BASE+'/ntuple/Latest'
            os.remove(ANA_BASE+'/ntuple/Latest')
        elif os.path.isdir(ANA_BASE+'/ntuple/Latest'):
            print 'Removing previous soft-links (if any) within: '+ANA_BASE+'/ntuple/Latest'
            for subdir in os.listdir(ANA_BASE+'/ntuple/Latest'):
                if os.path.islink(ANA_BASE+'/ntuple/Latest/'+subdir):
                    os.remove(ANA_BASE+'/ntuple/Latest/'+subdir)
    if not os.path.exists(ANA_BASE+'/ntuple/Latest'):
        print "Making directory for latest ntuple symlinks: "+ANA_BASE+'/ntuple/Latest'
        os.makedirs(ANA_BASE+'/ntuple/Latest')
    if os.path.exists(ANA_BASE+'/ntuple/eos_viktor/'+LATEST_NTUPLE_EOS):
        print 'Creating symlinks to the latest ntuples in Viktor\'s EOS folder: ntuple/eos_viktor/'+LATEST_NTUPLE_EOS+'/ ... ',
        for viktor_subdir in os.listdir(ANA_BASE+'/ntuple/eos_viktor/'+LATEST_NTUPLE_EOS):
            if os.path.isdir(ANA_BASE+'/ntuple/eos_viktor/'+LATEST_NTUPLE_EOS+'/'+viktor_subdir):
                source = os.path.realpath(ANA_BASE+'/ntuple/eos_viktor/'+LATEST_NTUPLE_EOS+'/'+viktor_subdir)
                target = ANA_BASE+'/ntuple/Latest/'+source.split("/")[-1]
                os.symlink(source, target)
        print 'Done.'
    if os.path.exists(ANA_BASE+'/ntuple/eos_janos/'+LATEST_NTUPLE_EOS):
        print 'Creating symlinks to the latest ntuples in Janos\' EOS folder: ntuple/eos_janos/'+LATEST_NTUPLE_EOS+'/ ... ',
        for janos_subdir in os.listdir(ANA_BASE+'/ntuple/eos_janos/'+LATEST_NTUPLE_EOS):
            if os.path.isdir(ANA_BASE+'/ntuple/eos_janos/'+LATEST_NTUPLE_EOS+'/'+janos_subdir):
                source = os.path.realpath(ANA_BASE+'/ntuple/eos_janos/'+LATEST_NTUPLE_EOS+'/'+janos_subdir)
                target = ANA_BASE+'/ntuple/Latest/'+source.split("/")[-1]
                if os.path.islink(target): target += "_2"
                os.symlink(source, target)
        print 'Done.'
    if os.path.exists(ANA_BASE+'/ntuple/eosuser_janos/'+LATEST_NTUPLE_EOS):
        print 'Creating symlinks to the latest ntuples in Janos\' EOS USER folder: ntuple/eosuser_janos/'+LATEST_NTUPLE_EOS+'/ ... ',
        for janos_subdir in os.listdir(ANA_BASE+'/ntuple/eosuser_janos/'+LATEST_NTUPLE_EOS):
            if os.path.isdir(ANA_BASE+'/ntuple/eosuser_janos/'+LATEST_NTUPLE_EOS+'/'+janos_subdir):
                source = os.path.realpath(ANA_BASE+'/ntuple/eosuser_janos/'+LATEST_NTUPLE_EOS+'/'+janos_subdir)
                target = ANA_BASE+'/ntuple/Latest/'+source.split("/")[-1]
                if os.path.islink(target): target += "_2"
                os.symlink(source, target)
        print 'Done.'
elif 'grid18.kfki.hu' in socket.gethostname():
    print 'Running on grid18 (Budapest)'

    if os.path.lexists(ANA_BASE+'/ntuple/Latest'):
        if os.path.islink(ANA_BASE+'/ntuple/Latest'):
            print 'Removing soft-link directory (used by previous version): '+ANA_BASE+'/ntuple/Latest'
            os.remove(ANA_BASE+'/ntuple/Latest')
        elif os.path.isdir(ANA_BASE+'/ntuple/Latest'):
            print 'Removing previous soft-links (if any) within: '+ANA_BASE+'/ntuple/Latest'
            for subdir in os.listdir(ANA_BASE+'/ntuple/Latest'):
                if os.path.islink(ANA_BASE+'/ntuple/Latest/'+subdir):
                    os.remove(ANA_BASE+'/ntuple/Latest/'+subdir)
    if not os.path.exists(ANA_BASE+'/ntuple/Latest'):
        print "Making directory for latest ntuple symlinks: "+ANA_BASE+'/ntuple/Latest'
        os.makedirs(ANA_BASE+'/ntuple/Latest')
    
    #print 'Creating symlinks to the latest ntuples in ntuple/grid18/'+LATEST_NTUPLE_GRID18+'/ ... ',
    #for grid18_subdir in os.listdir(ANA_BASE+'/ntuple/grid18/'+LATEST_NTUPLE_GRID18):
    #    if os.path.isdir(ANA_BASE+'/ntuple/grid18/'+LATEST_NTUPLE_GRID18+'/'+grid18_subdir):
    #        source = os.path.realpath(ANA_BASE+'/ntuple/grid18/'+LATEST_NTUPLE_GRID18+'/'+grid18_subdir)
    #        target = ANA_BASE+'/ntuple/Latest/'+source.split("/")[-1]
    #        os.symlink(source, target)
    #    print 'Done.'
    if os.path.exists(ANA_BASE+'/ntuple/grid18_data/'+LATEST_NTUPLE_GRID18):
        print 'Creating symlinks to the latest ntuples in /data drive: ntuple/grid18_data/'+LATEST_NTUPLE_GRID18+'/ ... ',
        for data_subdir in os.listdir(ANA_BASE+'/ntuple/grid18_data/'+LATEST_NTUPLE_GRID18):
            if os.path.isdir(ANA_BASE+'/ntuple/grid18_data/'+LATEST_NTUPLE_GRID18+'/'+data_subdir):
                source = os.path.realpath(ANA_BASE+'/ntuple/grid18_data/'+LATEST_NTUPLE_GRID18+'/'+data_subdir)
                target = ANA_BASE+'/ntuple/Latest/'+source.split("/")[-1]
                os.symlink(source, target)
        print 'Done.'
    if os.path.exists(ANA_BASE+'/ntuple/grid18_data_6tb/'+LATEST_NTUPLE_GRID18):
        print 'Creating symlinks to the latest ntuples in /data_6tb drive: ntuple/grid18_data_6tb/'+LATEST_NTUPLE_GRID18+'/ ... ',
        for data_6tb_subdir in os.listdir(ANA_BASE+'/ntuple/grid18_data_6tb/'+LATEST_NTUPLE_GRID18):
            if os.path.isdir(ANA_BASE+'/ntuple/grid18_data_6tb/'+LATEST_NTUPLE_GRID18+'/'+data_6tb_subdir):
                if not os.path.islink(ANA_BASE+'/ntuple/grid18_data_6tb/'+LATEST_NTUPLE_GRID18+'/'+data_6tb_subdir):
                    source = os.path.realpath(ANA_BASE+'/ntuple/grid18_data_6tb/'+LATEST_NTUPLE_GRID18+'/'+data_6tb_subdir)
                    target = ANA_BASE+'/ntuple/Latest/'+source.split("/")[-1]
                    os.symlink(source, target)
        print 'Done.'
else:
    print "Error: not on lxplus or grid18 (Budapest)"
    sys.exit()

print "Creating file lists ... ",
if not os.path.exists(ANA_BASE+'/filelists/data'): os.makedirs(ANA_BASE+'/filelists/data')
if not os.path.exists(ANA_BASE+'/filelists/signals'): os.makedirs(ANA_BASE+'/filelists/signals')
if not os.path.exists(ANA_BASE+'/filelists/backgrounds'): os.makedirs(ANA_BASE+'/filelists/backgrounds')
for txtfile in glob.glob('filelists/*/*.txt'):
    os.remove(txtfile)

for directory in os.listdir(DIR):
    if os.path.isdir(DIR+'/'+directory):
        if os.listdir(DIR+'/'+directory):
            # Data
            if re.compile('.*20[1-2][0-9][A-J].*').match(directory):
                txtname = directory
                if txtname.endswith('_recovery'): txtname = txtname[:-9]
                if txtname.endswith('_2'): txtname = txtname[:-2]
                flist = open(ANA_BASE+'/filelists/data/'+txtname+'.txt', 'a')
                for files in os.listdir(DIR+'/'+directory):
                    filename = os.path.realpath(DIR+'/'+directory+'/'+files)
                    if 'lxplus' in socket.gethostname():
                        if filename.startswith("/eos/user"):
                            filename = "root://eosuser/"+filename
                        else:
                            filename = "root://eoscms/"+filename
                    print>>flist, filename
            # Signals
            elif re.compile('.*T[1-9][t,b,c,q][t,b,c,q].*').match(directory):
                txtname = directory
                if txtname.endswith('_recovery'): txtname = txtname[:-9]
                if txtname.endswith('_2'): txtname = txtname[:-2]
                flist = open(ANA_BASE+'/filelists/signals/'+txtname+'.txt', 'a')
                for files in os.listdir(DIR+'/'+directory):
                    filename = os.path.realpath(DIR+'/'+directory+'/'+files)
                    if 'lxplus' in socket.gethostname():
                        if filename.startswith("/eos/user"):
                            filename = "root://eosuser/"+filename
                        else:
                            filename = "root://eoscms/"+filename
                    print>>flist, filename
            # Backgrounds
            else:
                txtname = directory
                # Please synchronize the lines below with get_xsec_totweight_from_txt_file in AnalysisBase.h
                if txtname.endswith('_recovery'): txtname = txtname[:-9]
                if txtname.endswith('_2'):      txtname = txtname[:-2]
                if txtname.endswith('_ext1'):   txtname = txtname[:-5]
                if txtname.endswith('_ext2'):   txtname = txtname[:-5]
                if txtname.endswith('_ext3'):   txtname = txtname[:-5]
                if txtname.endswith('_ext4'):   txtname = txtname[:-5]
                if txtname.endswith('_backup'): txtname = txtname[:-7]
                flist = open(ANA_BASE+'/filelists/backgrounds/'+txtname+'.txt', 'a')
                for files in os.listdir(DIR+'/'+directory):
                    filename = os.path.realpath(DIR+'/'+directory+'/'+files)
                    if 'lxplus' in socket.gethostname():
                        if filename.startswith("/eos/user"):
                            filename = "root://eosuser/"+filename
                        else:
                            filename = "root://eoscms/"+filename
                    print>>flist, filename

print "Creating temp file list directories (for batch and split jobs) ... ",
if not os.path.exists(ANA_BASE+'/filelists_tmp/data'): os.makedirs(ANA_BASE+'/filelists_tmp/data')
if not os.path.exists(ANA_BASE+'/filelists_tmp/signals'): os.makedirs(ANA_BASE+'/filelists_tmp/signals')
if not os.path.exists(ANA_BASE+'/filelists_tmp/backgrounds'): os.makedirs(ANA_BASE+'/filelists_tmp/backgrounds')

print 'Done.'
