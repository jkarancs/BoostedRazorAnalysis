import random as rd
import sys
import datetime

ngaus      =  18 # Number of gaussians to generate (8 sigma for syst, 10 sigma for SF)
nmax_scale =   3 # Number of scale variations to consider
nmax_pdf   = 100 # Number of pdfs to consider

def generate_line():
    # Here we generate one line in the systematics file
    # it consists of n gaussians, and then a uniform number for the pdfs

    sysline = ""
    for i in range(ngaus):
        number = rd.gauss(0,1)
        sysline = sysline + str(number) + " "

    # randrange returns integers
    scalenumber = rd.randrange(1,nmax_scale+1)
    pdfnumber   = rd.randrange(1,nmax_pdf+1)
    sysline = sysline + str(scalenumber) + " " + str(pdfnumber)

    return sysline

if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print "Usage: python %s nsyst" % (sys.argv[0])
        exit()

    nsyst = int(sys.argv[1])
    today = datetime.date.today()
    fname = "systematics/%s.txt" % (str(today).replace("-","_"))

    systfile = open(fname,'w')

    for i in range(nsyst):
        l = generate_line()
        systfile.write(l+"\n")

    systfile.close()
