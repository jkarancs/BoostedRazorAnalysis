import ROOT
import os

def set_default_style_():
    ROOT.gStyle.SetPaperSize(20.,20.);
    ROOT.gStyle.SetTitleFont(42,"xyz");
    ROOT.gStyle.SetCanvasBorderMode(0);
    ROOT.gStyle.SetCanvasColor(0);
    ROOT.gStyle.SetErrorX(0);
    ROOT.gStyle.SetFrameBorderMode(0);
    #ROOT.gStyle.SetFrameFillColor(0);
    #ROOT.gStyle.SetFrameFillStyle(0);
    ROOT.gStyle.SetFrameLineWidth(2);
    #ROOT.gStyle.SetLineWidth(2);
    ROOT.gStyle.SetOptStat(0);
    ROOT.gStyle.SetOptTitle(0);
    ROOT.gStyle.SetPadBorderMode(0);
    ROOT.gStyle.SetPadColor(0);
    ROOT.gStyle.SetPadTickX(1);
    ROOT.gStyle.SetPadTickY(1);
    ROOT.gStyle.SetPalette(1);
    ROOT.gStyle.SetTitleBorderSize(0);
    ROOT.gStyle.SetTitleFillColor(0);
    ROOT.gStyle.SetTitleStyle(0);
    ROOT.gStyle.SetTitleX(1);
    ROOT.gStyle.SetTitleY(1);
    ROOT.gStyle.SetTitleAlign(33);

def custom_can(h, canname, gx = 1, gy = 1,
               histosize_x = 500, histosize_y = 500,
               mar_left = 90, mar_right = 20, mar_top = 20, mar_bottom = 60,
               title_align = 33, title_x = 0.99, title_y = 0.99):
    #set_default_style_()
    if len(h.GetTitle()): mar_top += 25
    maxlabelsize = 0
    for binx in range(1, h.GetNbinsX()+1):
        label = h.GetXaxis().GetBinLabel(binx)
        if len(label)>maxlabelsize: maxlabelsize = len(label)
    if maxlabelsize>6: mar_bottom = maxlabelsize*10
    titlefontsize = 32
    labelfontsize = 20
    yoffset_x = mar_left - titlefontsize - 4
    xoffset_y = mar_bottom - titlefontsize - 4
    zoffset_x = mar_right - titlefontsize - 4
    padsize_x = histosize_x + mar_left + mar_right
    padsize_y = histosize_y + mar_top + mar_bottom
    padsize = min(padsize_x, padsize_y)
    padratio_yx = 1.0 if padsize_y/padsize_x else float(padsize_y)/padsize_x
    padratio_xy = 1.0 if padsize_x/padsize_y else float(padsize_x)/padsize_y
    xoffset = (float(xoffset_y)/titlefontsize+0.5) * padratio_xy /1.6
    yoffset = (float(yoffset_x)/titlefontsize+0.5) * padratio_yx /1.6
    zoffset = (float(zoffset_x)/titlefontsize+0.5) * padratio_yx /1.6
    titlesize = float(titlefontsize)/padsize
    labelsize = float(labelfontsize)/padsize
    if maxlabelsize>0: h.GetXaxis().SetLabelSize(1.5*labelsize)
    if (len(h.GetTitle())):
        ROOT.gStyle.SetOptTitle(1)
        ROOT.gStyle.SetTitleH(titlefontsize/padsize)
        ROOT.gStyle.SetTitleFontSize(titlesize*0.8)
        ROOT.gStyle.SetTitleBorderSize(0)
        ROOT.gStyle.SetTitleAlign(title_align)
        ROOT.gStyle.SetTitleX(title_x)
        ROOT.gStyle.SetTitleY(title_y)
    h.SetTitleFont(42,"xyz")
    h.SetLabelFont(42,"xyz")
    h.SetTitleSize(titlesize,"xyz")
    h.SetLabelSize(labelsize,"xyz")
    h.GetXaxis().SetTitleOffset(xoffset)
    h.GetYaxis().SetTitleOffset(yoffset)
    h.GetZaxis().SetTitleOffset(zoffset)
    h.GetYaxis().SetDecimals(1)
    h.GetZaxis().SetDecimals(1)
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetTitleH(titlefontsize/padsize)
    ROOT.gStyle.SetTitleFontSize(titlesize)
    canvas = ROOT.TCanvas(canname, h.GetTitle(), padsize_x + 4, padsize_y + 26)
    pad = canvas.cd(1)
    pad.SetLeftMargin(float(mar_left)/padsize_x)
    pad.SetRightMargin(float(mar_right)/padsize_x)
    pad.SetTopMargin(float(mar_top)/padsize_y)
    pad.SetBottomMargin(float(mar_bottom)/padsize_y)
    canvas.SetGrid(gx,gy)
    return canvas

def draw_mr_bins(h, ymin, ymax, combine_bins, keep,
                 mrbins  = [ 0.8, 1.0, 1.2, 1.6, 2.0, 4.0 ],
                 r2bins  = [ 0.08, 0.12, 0.16, 0.24, 0.4, 1.5  ]):
    for i in range(len(mrbins)-1):
        bins = range(i*5+1, i*5+6)
        if combine_bins and i==3: bins = range(16,20)
        if combine_bins and i==4: bins = range(20,23)
        maxcont = -9999
        for binx in bins:
            if h.GetBinContent(binx)>maxcont:
                maxcont = h.GetBinContent(binx)
        y2 = [maxcont+(ymax-ymin)*0.1, maxcont+(ymax-ymin)*0.175]
        if ymin!=0: y2 = [maxcont*((ymax/ymin)**0.1), maxcont*((ymax/ymin)**0.175)]
        x = i*5
        if i==4 and combine_bins: x = x-1
        if i!=0:
            line = ROOT.TLine(x,ymin,x,y2[0])
            line.SetLineStyle(2)
            line.Draw()
            keep.append(line)
        x = 2.5 + i*5
        if i==3 and combine_bins: x = x-0.5
        if i==4 and combine_bins: x = x-2
        if i==0:
            bin_lat = ROOT.TLatex(x, y2[1], "M_{R} (TeV)")
            bin_lat.SetTextAlign(22)
            bin_lat.SetTextSize(0.04)
            bin_lat.Draw()
            keep.append(bin_lat)
        lat = ROOT.TLatex(x, y2[0], ("[%1.1f, %1.1f]" % (mrbins[i], mrbins[i+1])))
        lat.SetTextAlign(22)
        lat.SetTextSize(0.04)
        lat.Draw()
        keep.append(lat)

def add_stack_ratio_plot(c, xmin, xmax, keep, add_labels=True, combine_bins=True,
                         mrbins  = [ 800, 1000, 1200, 1600, 2000, 4000 ],
                         r2bins  = [ 0.08, 0.12, 0.16, 0.24, 0.4, 1.5  ],
                         remove=False, debug = 0):
    # Canvas division sizes
    mar_top    = 45.
    y1         = 365.
    mid2       = 10.
    y2         = 115.
    mar_bottom = 100.
    mar_left   = 90.
    x          = 500.
    mar_right  = 20.
    x_can    = mar_left+x+mar_right
    y_can    = mar_top+y1+mid2*2+y2+mar_bottom
    padsize1 = min(mar_top+y1+mid2,    x_can)
    padsize2 = min(mar_bottom+y2+mid2, x_can)
    labelfontsize = 20.
    titlefontsize = 32.
    leg_y2 = 0.9 # not used values, read from orig
    ok = False
    if debug: print "Start debugging: "+c.GetName()
    if debug: print "ok"
    if (c.GetListOfPrimitives().GetEntries()>2):
        # Histos
        if debug: print "ok1"
        Data = c.GetListOfPrimitives().At(1)
        if debug: print "ok1"
        MCstack = c.GetListOfPrimitives().At(2)
        if debug: print "ok1"
        syst_err = c.GetListOfPrimitives().At(3)
        if debug: print "ok1"
        leg = c.GetListOfPrimitives().At(c.GetListOfPrimitives().GetEntries()-1)
        if debug: print "ok1"
        if not MCstack.GetTitle()=="0":
            ratio = Data.Clone(Data.GetName()+"_num")
            ratio.SetDirectory(0)
            if debug: print "ok2"
            mc_sum = MCstack.GetHists().At(0).Clone(Data.GetName()+"_den")
            mc_sum.SetDirectory(0)
            mc_sum_syst = 0
            if debug: print "ok2"
            for iStack in range(1, MCstack.GetHists().GetEntries()):
                h = MCstack.GetHists().At(iStack)
                mc_sum.Add(h.Clone())
            if debug: print "ok2"
            den_stat_err     = mc_sum.Clone("den_stat_err")
            if debug: print "ok2"
            den_total_err = syst_err.Clone("den_total_err")
            if debug: print "ok2"
            # Instead of Divide(), scale the error of num, and plot error of den around 1
            ratio.Divide(mc_sum)
            if debug: print "ok2"
            for bin in range(1, ratio.GetNbinsX()+1):
                if (mc_sum.GetBinContent(bin)!=0):
                    ratio  .SetBinContent(bin, Data.GetBinContent(bin)/mc_sum.GetBinContent(bin))
                    ratio  .SetBinError  (bin, Data.GetBinError(bin)  /mc_sum.GetBinContent(bin))
                    den_stat_err.SetBinContent(bin, 1)
                    den_stat_err.SetBinError  (bin, mc_sum.GetBinError(bin)  /mc_sum.GetBinContent(bin))
                    den_total_err.SetBinContent(bin, den_total_err.GetBinContent(bin)/mc_sum.GetBinContent(bin))
                    den_total_err.SetBinError  (bin, den_total_err.GetBinError(bin)  /mc_sum.GetBinContent(bin))
                else:
                    ratio  .SetBinContent(bin, 0)
                    ratio  .SetBinError  (bin, 0)
                    den_stat_err.SetBinContent(bin, 1)
                    den_stat_err.SetBinError  (bin, 0)
                    den_total_err.SetBinContent(bin, 1)
                    den_total_err.SetBinError  (bin, 0)
            if debug: print "ok2"
            # Legend
            # Remove Non-Data non-stack plots (eg. signal)
            # indices:
            # 0: Data, 1: stack, 2: Data again, 3+: (signals), 3+nsig: Legend
            if debug: print "ok2"
            # Styles
            heightratio1 = float(padsize1)/y_can
            Data .SetTitleSize  (Data.GetYaxis().GetTitleSize()  /heightratio1,"y")
            Data .SetTitleOffset(Data.GetYaxis().GetTitleOffset()*heightratio1,"y")
            Data .SetLabelSize(labelfontsize/padsize1,"xyz")
            ratio.SetLabelSize(labelfontsize/padsize2,"xyz")
            ratio.SetTitleSize(titlefontsize/padsize2,"xyz")
            ratio.GetYaxis().SetRangeUser(0,3)
            ratio.GetYaxis().SetNdivisions(305)
            ratio.GetYaxis().SetTitle("#frac{Data}{Estimate}")
            if debug: print "ok2"
            heightratio2 = float(padsize2)/y_can
            #ratio.SetTitleOffset(ratio.GetYaxis().GetTitleOffset()*heightratio2,"y")
            ratio.GetYaxis().SetTitleOffset(0.33)
            ratio.SetTitle("")
            ratio.SetMarkerStyle(20)
            ratio.SetMarkerColor(1)
            ratio.SetLineColor(1)
            if debug: print "ok2"
            # New Canvas
            left_mar = c.GetLeftMargin()
            right_mar = c.GetRightMargin()
            logScale = c.GetLogy()
            c = ROOT.TCanvas(c.GetName()+"_Ratio", c.GetTitle(), int(x_can+4), int(y_can+26)) # 600, 600
            c.Divide(1,2)
            if debug: print "ok2"
            # Pad 1 (x: 90+500+20 x y: 45+350+10)
            p = c.cd(1)
            p.SetGrid(c.GetGridx(),c.GetGridy())
            p.SetPad(0,float(padsize2)/y_can,1,1)
            if debug: print "ok2"
            p.SetTopMargin(mar_top/(mar_top+y1+mid2))
            p.SetBottomMargin(0)
            p.SetLeftMargin(left_mar)
            p.SetRightMargin(right_mar)
            if debug: print "ok2"
            if (logScale): p.SetLogy(1)
            Data.Draw("PE1")
            MCstack.Draw("SAMEHIST")
            syst_err.Draw("SAME E2")
            leg.Draw("SAME")
            if debug: print "ok2"
            Data.Draw("SAMEPE1")
            if debug: print "ok3"
            draw_mr_bins(Data, Data.GetMinimum(),Data.GetMaximum(), combine_bins, keep, mrbins)
            if debug: print "ok3"
            ROOT.gPad.Update()
            if debug: print "ok3"
            # Pad 2 (x: 90+500+20 x y: 60+150+10)
            p2 = c.cd(2)
            p2.SetPad(0,0,1,float(padsize2)/y_can)
            p2.SetLogy(0)
            p2.SetGrid(0,1)
            p2.SetTopMargin(float(mid2)/padsize2)
            p2.SetBottomMargin(float(mar_bottom)/padsize2)
            p2.SetLeftMargin(left_mar)
            p2.SetRightMargin(right_mar)
            if debug: print "ok3"
            den_stat_err.SetFillColor(1)
            den_stat_err.SetFillStyle(3004)
            den_stat_err.SetMarkerStyle(0)
            den_stat_err.SetMarkerColor(0)
            den_total_err.SetFillColor(ROOT.kGray) # 920
            den_total_err.SetFillStyle(1001)
            den_total_err.SetMarkerStyle(0)
            den_total_err.SetMarkerColor(0)
            ratio.Draw("PE1")
            den_total_err.Draw("SAME E2")
            den_stat_err.Draw("SAME E2")
            ratio.Draw("SAME PE1")
            if debug: print "ok3"
            if (xmin==xmax):
                xmin = ratio.GetYaxis().GetXmin()
                xmax = ratio.GetYaxis().GetXmax()
            if debug: print "ok3"
            l = ROOT.TLine(xmin, 1, xmax, 1)
            l.SetLineWidth(2)
            #l.SetLineColor(2)
            l.SetLineStyle(2)
            l.Draw()
            ROOT.gPad.Update()
            if add_labels:
                Razor_labels = []
                binx = 0
                if combine_bins:
                    for i in range(len(mrbins)-1):
                        for j in range(len(r2bins)-1):
                            binx += 1
                            if binx<=18:
                                Razor_labels.append("[%.2f, %.2f]" % (r2bins[j], r2bins[j+1]) )
                            elif binx==19:
                                Razor_labels.append("[%.2f, %.2f]" % (r2bins[j], r2bins[j+2]) )
                            elif binx>=21 and binx<=22:
                                Razor_labels.append("[%.2f, %.2f]" % (r2bins[j], r2bins[j+1]) )
                            elif binx==23:
                                Razor_labels.append("[%.2f, %.2f]" % (r2bins[j], r2bins[j+3]) )
                else:
                    for i in range(len(mrbins)-1):
                        for j in range(len(r2bins)-1):
                            binx += 1
                            Razor_labels.append("[%.2f, %.2f]" % r2bins[j], r2bins[j+1])
                ratio.GetXaxis().SetLabelColor(0)
                labelsize = ratio.GetXaxis().GetLabelSize()
                ymin = ratio.GetMinimum()
                ymax = ratio.GetMaximum()
                offset = (ymax-ymin) * ratio.GetXaxis().GetLabelOffset() * 5
                for i in range(len(Razor_labels)):
                    lat = ROOT.TLatex(0.5+i, ymin-offset, Razor_labels[i])
                    lat.SetTextAlign(32)
                    lat.SetTextAngle(90)
                    lat.SetTextFont(ratio.GetXaxis().GetLabelFont())
                    lat.SetTextSize(labelsize)
                    lat.Draw("SAME")
                    ROOT.SetOwnership(lat, False)
            c.Write()
            if debug: print "ok3"
            ok = 1

def add_ratio_plot(c, xmin, xmax, keep, add_labels=True, combine_bins=True,
                   mrbins  = [ 800, 1000, 1200, 1600, 2000, 4000 ],
                   r2bins  = [ 0.08, 0.12, 0.16, 0.24, 0.4, 1.5  ],
                   yratio = 1.0, remove=False, debug = 0):
    # Canvas division sizes
    mar_top    = 45.
    y1         = 365.
    mid2       = 10.
    y2         = 115.
    mar_bottom = 100.
    mar_left   = 90.
    x          = 500.
    mar_right  = 20.
    x_can    = mar_left+x+mar_right
    y_can    = mar_top+y1+mid2*2+y2+mar_bottom
    padsize1 = min(mar_top+y1+mid2,    x_can)
    padsize2 = min(mar_bottom+y2+mid2, x_can)
    labelfontsize = 20.
    titlefontsize = 32.
    leg_y2 = 0.9 # not used values, read from orig
    ok = False
    if debug: print "Start debugging: "+c.GetName()
    if debug: print "ok"
    # Histos
    if debug: print "ok1"
    num = c.GetListOfPrimitives().At(0)
    if debug: print "ok1"
    den = c.GetListOfPrimitives().At(1)
    if debug: print "ok1"
    den2 = c.GetListOfPrimitives().At(2)
    if debug: print "ok1"
    leg = c.GetListOfPrimitives().At(c.GetListOfPrimitives().GetEntries()-1)
    if debug: print "ok1"
    ratio = num.Clone(num.GetName()+"_ratio")
    ratio.SetDirectory(0)
    if debug: print "ok2"
    den_stat_err     = den.Clone("den_stat_err")
    if debug: print "ok2"
    # Instead of Divide(), scale the error of num, and plot error of den around 1
    #ratio.Divide(den)
    if debug: print "ok2"
    ratiomax = 0
    for binx in range(1, ratio.GetNbinsX()+1):
        if (den.GetBinContent(binx)!=0):
            if num.GetBinContent(binx)/den.GetBinContent(binx)>ratiomax:
                ratiomax = num.GetBinContent(binx)/den.GetBinContent(binx)
            ratio  .SetBinContent(binx, num.GetBinContent(binx) /den.GetBinContent(binx))
            ratio  .SetBinError  (binx, num.GetBinError(binx)   /den.GetBinContent(binx))
            den_stat_err.SetBinContent(binx, yratio)
            den_stat_err.SetBinError  (binx, yratio*den.GetBinError(binx)/den.GetBinContent(binx))
        else:
            ratio  .SetBinContent(binx, 0)
            ratio  .SetBinError  (binx, 0)
            den_stat_err.SetBinContent(binx, 0)
            den_stat_err.SetBinError  (binx, 0)
    if debug: print "ok2"
    # Legend
    # indices:
    # 0: num, 1: mc, 2: Legend
    if debug: print "ok2"
    # Styles
    heightratio1 = float(padsize1)/y_can
    num .SetTitleSize  (num.GetYaxis().GetTitleSize()  /heightratio1,"y")
    num .SetTitleOffset(num.GetYaxis().GetTitleOffset()*heightratio1,"y")
    num .SetLabelSize(labelfontsize/padsize1,"xyz")
    ratio.SetLabelSize(labelfontsize/padsize2,"xyz")
    ratio.SetTitleSize(titlefontsize/padsize2,"xyz")
    #ratio.GetYaxis().SetRangeUser(0,3*yratio)
    if yratio == 1.0:
        ratio.GetYaxis().SetRangeUser(0,3)
    else:
        ratio.GetYaxis().SetRangeUser(0,(min(4,int(ratiomax/yratio))+1)*yratio)
    ratio.GetYaxis().SetNdivisions(305)
    ratio.GetYaxis().SetTitle("Ratio")
    if debug: print "ok2"
    heightratio2 = float(padsize2)/y_can
    #ratio.SetTitleOffset(ratio.GetYaxis().GetTitleOffset()*heightratio2,"y")
    ratio.GetYaxis().SetTitleOffset(0.5)
    ratio.SetTitle("")
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerColor(1)
    ratio.SetLineColor(1)
    if debug: print "ok2"
    # New Canvas
    left_mar = c.GetLeftMargin()
    right_mar = c.GetRightMargin()
    logScale = c.GetLogy()
    can2 = ROOT.TCanvas(c.GetName()+"_Ratio", c.GetTitle(), int(x_can+4), int(y_can+26)) # 600, 600
    can2.Divide(1,2)
    if debug: print "ok2"
    # Pad 1 (x: 90+500+20 x y: 45+350+10)
    p = can2.cd(1)
    p.SetGrid(c.GetGridx(),c.GetGridy())
    p.SetPad(0,float(padsize2)/y_can,1,1)
    if debug: print "ok2"
    p.SetTopMargin(mar_top/(mar_top+y1+mid2))
    p.SetBottomMargin(0)
    p.SetLeftMargin(left_mar)
    p.SetRightMargin(right_mar)
    if debug: print "ok2"
    if (logScale): p.SetLogy(1)
    num.Draw("PE1")
    den.Draw("SAME HIST")
    den2.Draw("SAME PE1")
    leg.Draw("SAME")
    #if debug: print "ok2"
    #num.Draw("SAME PE1")
    draw_mr_bins(num, num.GetMinimum(),num.GetMaximum(), combine_bins, keep, mrbins)
    if debug: print "ok3"
    ROOT.gPad.Update()
    if debug: print "ok3"
    # Pad 2 (x: 90+500+20 x y: 60+150+10)
    p2 = can2.cd(2)
    p2.SetPad(0,0,1,float(padsize2)/y_can)
    p2.SetLogy(0)
    p2.SetGrid(0,(yratio==1.0))
    p2.SetTopMargin(float(mid2)/padsize2)
    p2.SetBottomMargin(float(mar_bottom)/padsize2)
    p2.SetLeftMargin(left_mar)
    p2.SetRightMargin(right_mar)
    if debug: print "ok3"
    ratio.Draw("PE1")
    den_stat_err.SetFillColor(1)
    den_stat_err.SetFillStyle(3004)
    den_stat_err.SetMarkerStyle(0)
    den_stat_err.Draw("SAME E2")
    ratio.Draw("SAME PE1")
    if debug: print "ok3"
    if (xmin==xmax):
        xmin = ratio.GetYaxis().GetXmin()
        xmax = ratio.GetYaxis().GetXmax()
    if debug: print "ok3"
    l = ROOT.TLine(xmin, yratio, xmax, yratio)
    l.SetLineWidth(2)
    #l.SetLineColor(2)
    l.SetLineStyle(2)
    l.Draw()
    ROOT.gPad.Update()
    if add_labels:
        add_r2_labels(ratio, combine_bins, mrbins, r2bins)
        #Razor_labels = []
        #binx = 0
        #if combine_bins:
        #    for i in range(len(mrbins)-1):
        #        for j in range(len(r2bins)-1):
        #            binx += 1
        #            if binx<=18:
        #                Razor_labels.append("[%.2f, %.2f]" % (r2bins[j], r2bins[j+1]) )
        #            elif binx==19:
        #                Razor_labels.append("[%.2f, %.2f]" % (r2bins[j], r2bins[j+2]) )
        #            elif binx>=21 and binx<=22:
        #                Razor_labels.append("[%.2f, %.2f]" % (r2bins[j], r2bins[j+1]) )
        #            elif binx==23:
        #                Razor_labels.append("[%.2f, %.2f]" % (r2bins[j], r2bins[j+3]) )
        #else:
        #    for i in range(len(mrbins)-1):
        #        for j in range(len(r2bins)-1):
        #            binx += 1
        #            Razor_labels.append("[%.2f, %.2f]" % r2bins[j], r2bins[j+1])
        #ratio.GetXaxis().SetLabelColor(0)
        #labelsize = ratio.GetXaxis().GetLabelSize()
        #ymin = ratio.GetMinimum()
        #ymax = ratio.GetMaximum()
        #offset = (ymax-ymin) * ratio.GetXaxis().GetLabelOffset() * 5
        #for i in range(len(Razor_labels)):
        #    lat = ROOT.TLatex(0.5+i, ymin-offset, Razor_labels[i])
        #    lat.SetTextAlign(32)
        #    lat.SetTextAngle(90)
        #    lat.SetTextFont(ratio.GetXaxis().GetLabelFont())
        #    lat.SetTextSize(labelsize)
        #    lat.Draw("SAME")
        #    ROOT.SetOwnership(lat, False)
    can2.Write()
    if debug: print "ok3"
    ok = 1
    return can2

def add_r2_labels(h, combine_bins, 
                  mrbins  = [ 800, 1000, 1200, 1600, 2000, 4000 ],
                  r2bins  = [ 0.08, 0.12, 0.16, 0.24, 0.4, 1.5  ]):
    Razor_labels = []
    binx = 0
    if combine_bins:
        for i in range(len(mrbins)-1):
            for j in range(len(r2bins)-1):
                binx += 1
                if binx<=18:
                    Razor_labels.append("[%.2f, %.2f]" % (r2bins[j], r2bins[j+1]) )
                elif binx==19:
                    Razor_labels.append("[%.2f, %.2f]" % (r2bins[j], r2bins[j+2]) )
                elif binx>=21 and binx<=22:
                    Razor_labels.append("[%.2f, %.2f]" % (r2bins[j], r2bins[j+1]) )
                elif binx==23:
                    Razor_labels.append("[%.2f, %.2f]" % (r2bins[j], r2bins[j+3]) )
    else:
        for i in range(len(mrbins)-1):
            for j in range(len(r2bins)-1):
                binx += 1
                Razor_labels.append("[%.2f, %.2f]" % r2bins[j], r2bins[j+1])
    h.GetXaxis().SetLabelColor(0)
    labelsize = h.GetXaxis().GetLabelSize()
    ymin = h.GetMinimum()
    ymax = h.GetMaximum()
    offset = (ymax-ymin) * h.GetXaxis().GetLabelOffset() * 5
    for i in range(len(Razor_labels)):
        #print str(i)+" "+str(ymin-offset)
        lat = ROOT.TLatex(0.5+i, ymin-offset, Razor_labels[i])
        lat.SetTextAlign(32)
        lat.SetTextAngle(90)
        lat.SetTextFont(h.GetXaxis().GetLabelFont())
        lat.SetTextSize(labelsize)
        lat.Draw("SAME")

def add_bin_labels(h, combine_bins,
                   mrbins  = [ 800, 1000, 1200, 1600, 2000, 4000 ],
                   r2bins  = [ 0.08, 0.12, 0.16, 0.24, 0.4, 1.5  ]):
    binx = 0
    if combine_bins:
        for i in range(len(mrbins)-1):
            for j in range(len(r2bins)-1):
                binx += 1
                if binx<=18:
                    h.GetXaxis().SetBinLabel(binx,"[%.2f, %.2f]" % (r2bins[j], r2bins[j+1]))
                elif binx==19:
                    h.GetXaxis().SetBinLabel(binx,"[%.2f, %.2f]" % (r2bins[j], r2bins[j+2]))
                elif binx>=21 and binx<=22:
                    h.GetXaxis().SetBinLabel(binx-1,"[%.2f, %.2f]" % (r2bins[j], r2bins[j+1]))
                elif binx==23:
                    h.GetXaxis().SetBinLabel(binx-1,"[%.2f, %.2f]" % (r2bins[j], r2bins[j+3]))
    else:
        for i in range(len(mrbins)-1):
            for j in range(len(r2bins)-1):
                binx += 1
                h.GetXaxis().SetBinLabel(binx,"[%.2f, %.2f]" % (r2bins[j], r2bins[j+1]))
    h.GetXaxis().LabelsOption("v")
    h.GetXaxis().SetTitle("")

# Silence stdout/stderr
class suppress_stdout_stderr(object):
    def __init__(self):
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        self.save_fds = [os.dup(1), os.dup(2)]
    def __enter__(self):
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)
    def __exit__(self, *_):
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        for fd in self.null_fds + self.save_fds:
            os.close(fd)

def save_plot(can, name, plotname, write=True):
    # Check if the directory exists (if not create it first)
    dirname = os.path.dirname(plotname)
    if dirname != "" and not os.path.exists(dirname):
        os.makedirs(dirname)
    with suppress_stdout_stderr():
        can.SaveAs(plotname+".png")
    if write:
        if name != "":
            can.Write(name)
        else:
            can.Write()
