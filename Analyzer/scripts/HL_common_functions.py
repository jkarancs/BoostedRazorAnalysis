import ROOT
import os, sys, subprocess, time, math

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

def draw_mr_bins(vh, ymin, ymax, combine_bins, keep,
                 mrbins  = [ 0.8, 1.0, 1.2, 1.6, 2.0, 4.0 ],
                 r2bins  = [ 0.08, 0.12, 0.16, 0.24, 0.4, 1.5  ],
                 exceptions = {},
                 xoffset = 0, yoffset=0, textsize=0.04):
    for i in range(len(mrbins)-1):
        bins = range(i*5+1, i*5+6)
        if combine_bins and i==3: bins = range(16,20)
        if combine_bins and i==4: bins = range(20,23)
        maxcont = -9999
        for binx in bins:
            for h in vh:
                if h.InheritsFrom("TH1"):
                    if h.GetBinContent(binx)>maxcont:
                        maxcont = h.GetBinContent(binx)
                elif h.InheritsFrom("THStack"):
                    sum = 0;
                    for j in range(h.GetHists().GetEntries()):
                        sum = sum + max(0.,h.GetHists().At(j).GetBinContent(binx))
                    if sum>maxcont: maxcont = sum
        y2 = maxcont+(ymax-ymin)*0.1 + yoffset
        if ymin!=0: y2 = maxcont*((ymax/ymin)**0.1 + yoffset)
        for iexc, y2exc in exceptions.iteritems():
            if i==iexc: y2 = y2exc
        y3 = y2 +(ymax-ymin)*0.1
        if ymin!=0: y3 = y2*((ymax/ymin)**0.075)
        x = i*5
        if i==4 and combine_bins: x = x-1
        if i!=0:
            line = ROOT.TLine(x,ymin,x,y2)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()
            keep.append(line)
        x = 2.5 + i*5
        if i==3 and combine_bins: x = x-0.5
        if i==4 and combine_bins: x = x-2
        if i==0:
            bin_lat = ROOT.TLatex(x+xoffset, y3, "M_{R} (TeV)")
            bin_lat.SetTextAlign(22)
            bin_lat.SetTextSize(textsize)
            bin_lat.Draw()
            keep.append(bin_lat)
        num1 = "%d" if mrbins[i]==float(int(mrbins[i])) else "%1.1f"
        num2 = "%d" if mrbins[i+1]==float(int(mrbins[i+1])) else "%1.1f"
        lat = ROOT.TLatex(x+xoffset if i<4 else x, y2, ("["+num1+", "+num2+"]") % (mrbins[i], mrbins[i+1]))
        lat.SetTextAlign(22)
        lat.SetTextSize(textsize)
        lat.Draw()
        keep.append(lat)

def add_stack_ratio_plot(c, xmin, xmax, keep, add_labels=True, combine_bins=True,
                         mrbins  = [ 800, 1000, 1200, 1600, 2000, 4000 ],
                         r2bins  = [ 0.08, 0.12, 0.16, 0.24, 0.4, 1.5  ],
                         exceptions = {},
                         remove=False, legx1 = 0.16, debug = 0):
    # Canvas division sizes
    mar_top    = 45.
    y1         = 365.
    mid2       = 10.
    y2         = 115.
    mar_bottom = 132.
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
        stat_err = c.GetListOfPrimitives().At(4)
        if debug: print "ok1"
        for i in range(c.GetListOfPrimitives().GetEntries()):
            prim = c.GetListOfPrimitives().At(i)
            if prim.GetTitle().startswith("Legend"):
                leg = prim
                break
        # Add also signals if they are there
        vh_signals = []
        if "BkgEstimate" in c.GetName():
            for i in range(c.GetListOfPrimitives().GetEntries()):
                if "signal" in c.GetListOfPrimitives().At(i).GetName():
                    vh_signals.append(c.GetListOfPrimitives().At(i))
        if debug: print "ok1"
        if not MCstack.GetTitle()=="0":
            ratio = Data.Clone(Data.GetName()+"_num")
            keep.append(ratio)
            if debug: print "ok2"
            mc_sum = MCstack.GetHists().At(0).Clone(Data.GetName()+"_den")
            keep.append(mc_sum)
            mc_sum_syst = 0
            if debug: print "ok2"
            for iStack in range(1, MCstack.GetHists().GetEntries()):
                h = MCstack.GetHists().At(iStack)
                mc_sum.Add(h.Clone())
            if debug: print "ok2"
            den_stat_err     = mc_sum.Clone("den_stat_err")
            keep.append(den_stat_err)
            if debug: print "ok2"
            den_total_err = syst_err.Clone("den_total_err")
            keep.append(den_total_err)
            if debug: print "ok2"
            # Instead of Divide(), scale the error of num, and plot error of den around 1
            ratio.Divide(mc_sum)
            if debug: print "ok2"
            for bin in range(1, ratio.GetNbinsX()+1):
                if (mc_sum.GetBinContent(bin)!=0):
                    ratio  .SetBinContent(bin, Data.GetBinContent(bin)/mc_sum.GetBinContent(bin))
                    if Data.GetBinError(bin)==0:
                        #Draw the Garwood confidence interval for 0 counts [0,1.83]
                        ratio.SetBinError(bin, 1.83 /mc_sum.GetBinContent(bin))
                    else:
                        ratio.SetBinError(bin, Data.GetBinError(bin)/mc_sum.GetBinContent(bin))
                    #den_stat_err.SetBinContent(bin, 1)
                    den_stat_err.SetBinContent(bin, 0)
                    den_stat_err.SetBinError  (bin, mc_sum.GetBinError(bin)  /mc_sum.GetBinContent(bin))
                    #den_total_err.SetBinContent(bin, den_total_err.GetBinContent(bin)/mc_sum.GetBinContent(bin))
                    den_total_err.SetBinContent(bin, den_total_err.GetBinContent(bin)/mc_sum.GetBinContent(bin) - 1)
                    den_total_err.SetBinError  (bin, den_total_err.GetBinError(bin)  /mc_sum.GetBinContent(bin))
                else:
                    ratio  .SetBinContent(bin, 0)
                    ratio  .SetBinError  (bin, 0)
                    den_stat_err.SetBinContent(bin, 1)
                    den_stat_err.SetBinError  (bin, 0)
                    #den_total_err.SetBinContent(bin, 1)
                    den_total_err.SetBinContent(bin, 0)
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
            max_range = int(math.ceil(ratio.GetBinContent(ratio.GetMaximumBin())))
            ratio.GetYaxis().SetRangeUser(-0.5,0.5)
            ratio.GetYaxis().SetNdivisions(504)
            #ratio.GetYaxis().SetRangeUser(0,max_range)
            #ratio.GetYaxis().SetNdivisions(501+max_range)
            #ratio.GetYaxis().SetTitle("#frac{Data}{Estimate}")
            ratio.GetYaxis().SetTitle("Rel. unc.")
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
            c = ROOT.TCanvas(c.GetName()+"_Ratio", c.GetTitle(), int(x_can+4), int(y_can+26)) # 600, 600
            keep.append(c)
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
            Data.Draw("AXIS")
            MCstack.Draw("SAME HIST")
            #syst_err.Draw("SAME E2")
            #stat_err.Draw("SAME E2")
            for h_signal in vh_signals: h_signal.Draw("SAME HIST")
            leg.Draw("SAME")
            if debug: print "ok2"
            #Data.Draw("SAMEPE0")
            if debug: print "ok3"            
            # Draw also Garwood intervals for 0 counts [0,1.83]
            zero = Data.Clone(Data.GetName()+"_zeroes")
            keep.append(zero)
            ymin = zero.GetMinimum()
            if ymin>0: zero.SetMarkerStyle(1)
            for binx in range(1,zero.GetNbinsX()+1):
                if zero.GetBinContent(binx)>0:
                    zero.SetBinContent(binx,0)
                    zero.SetBinError  (binx,0)
                else:
                    if ymin>0: zero.SetBinContent(binx, ymin*1.000001)
                    zero.SetBinError(binx, 1.83-ymin*1.000001)
            #zero.Draw("SAME PE")
            if debug: print "ok3"
            draw_mr_bins([Data, MCstack], Data.GetMinimum(),Data.GetMaximum(), combine_bins, keep, mrbins, r2bins, exceptions)
            if debug: print "ok3"
            ROOT.gPad.RedrawAxis()
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
            #den_stat_err.SetFillColor(1)
            #den_stat_err.SetFillStyle(3004)
            #den_stat_err.SetMarkerStyle(0)
            #den_stat_err.SetMarkerColor(0)
            #den_total_err.SetFillColor(ROOT.kGray) # 920
            #den_total_err.SetFillStyle(1001)
            #den_total_err.SetMarkerStyle(0)
            #den_total_err.SetMarkerColor(0)
            #den_total_err.SetFillColor(13)
            #den_total_err.SetFillStyle(3001)
            den_stat_err.SetLineColor(1)
            den_stat_err.SetMarkerStyle(0)
            den_stat_err.SetMarkerColor(0)
            den_stat_err.SetFillColor(1)
            den_stat_err.SetFillStyle(3004)
            den_total_err.SetLineColor(1)
            den_total_err.SetMarkerStyle(0)
            den_total_err.SetMarkerColor(0)
            den_total_err.SetFillColor(ROOT.kGray)
            den_total_err.SetFillStyle(1001)
            ratio.Draw("AXIS")
            den_total_err.Draw("SAME E2")
            den_stat_err.Draw("SAME E2")
            #ratio.Draw("SAME PE0")
            if debug: print "ok3"
            if (xmin==xmax):
                xmin = ratio.GetXaxis().GetXmin()
                xmax = ratio.GetXaxis().GetXmax()
            if debug: print "ok3"
            l = ROOT.TLine(xmin, 1, xmax, 1)
            l.SetLineWidth(2)
            #l.SetLineColor(2)
            l.SetLineStyle(2)
            l.Draw()
            keep.append(l)
            if debug: print "ok3"
            # Add legend to indicate stat/total error
            legx1 = 0.16
            legx2 = legx1 + 0.34
            legy2 = 1.0   - (15.0 / (mid2+y2+mar_bottom))
            legy1 = legy2 - (25.0 / (mid2+y2+mar_bottom))
            leg2 = ROOT.TLegend(legx1, legy1, legx2, legy2, "")
            leg2.SetNColumns(2)
            leg2.SetFillColor(0)
            leg2.SetFillStyle(0)
            leg2.SetBorderSize(0)
            leg2.SetTextSize(12.5/(mid2+y2+mar_bottom))
            leg2.AddEntry(den_stat_err,  "Stat. unc.",         "f")
            leg2.AddEntry(den_total_err, "Stat. + syst. unc.", "f")
            leg2.Draw("SAME")
            keep.append(leg2)
            if debug: print "ok3"
            ROOT.gPad.Update()
            if debug: print "ok3"
            ROOT.gPad.RedrawAxis()
            if add_labels: add_r2_labels(ratio, combine_bins, keep, mrbins, r2bins)
            #c.Write()
            if debug: print "ok3"
            ok = 1
    return c

def add_ratio_plot(c, xmin, xmax, keep, add_labels=True, combine_bins=True,
                   mrbins  = [ 800, 1000, 1200, 1600, 2000, 4000 ],
                   r2bins  = [ 0.08, 0.12, 0.16, 0.24, 0.4, 1.5  ],
                   yratio = 1.0, remove=False, debug = 0):
    # Canvas division sizes
    mar_top    = 45.
    y1         = 365.
    mid2       = 10.
    y2         = 115.
    mar_bottom = 132.
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
    for i in range(c.GetListOfPrimitives().GetEntries()):
        prim = c.GetListOfPrimitives().At(i)
        if prim.GetTitle().startswith("Legend"):
            leg = prim
            break
    keep.append(leg)
    if debug: print "ok1"
    ratio = num.Clone(num.GetName()+"_ratio")
    keep.append(ratio)
    if debug: print "ok2"
    #den_stat_err     = den.Clone("den_stat_err")
    #keep.append(den_stat_err)
    #if debug: print "ok2"
    # Instead of Divide(), scale the error of num, and plot error of den around 1
    #ratio.Divide(den)
    if debug: print "ok2"
    ratiomax = 0
    for binx in range(1, ratio.GetNbinsX()+1):
        if (den.GetBinContent(binx)!=0):
            ymax = (num.GetBinContent(binx)+num.GetBinError(binx))/den.GetBinContent(binx)
            if ymax>ratiomax: ratiomax = ymax
            ratio  .SetBinContent(binx, num.GetBinContent(binx) /den.GetBinContent(binx))
            ratio  .SetBinError  (binx, num.GetBinError(binx)   /den.GetBinContent(binx))
            #den_stat_err.SetBinContent(binx, yratio)
            #den_stat_err.SetBinError  (binx, yratio*den.GetBinError(binx)/den.GetBinContent(binx))
        else:
            ratio  .SetBinContent(binx, 0)
            ratio  .SetBinError  (binx, 0)
            #den_stat_err.SetBinContent(binx, 0)
            #den_stat_err.SetBinError  (binx, 0)
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
        ratio.GetYaxis().SetRangeUser(0,(int(ratiomax/yratio)+1)*yratio)
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
    can2.SetGrid(c.GetGridx(),c.GetGridy())
    keep.append(can2)
    can2.Divide(1,2)
    if debug: print "ok2"
    # Pad 1 (x: 90+500+20 x y: 45+350+10)
    p = can2.cd(1)
    p.SetPad(0,float(padsize2)/y_can,1,1)
    if debug: print "ok2"
    p.SetTopMargin(mar_top/(mar_top+y1+mid2))
    p.SetBottomMargin(0)
    p.SetLeftMargin(left_mar)
    p.SetRightMargin(right_mar)
    if debug: print "ok2"
    if (logScale): p.SetLogy(1)
    num.Draw("PE0")
    den.Draw("SAME HIST")
    den2.Draw("SAME PE0")
    leg.Draw("SAME")
    #if debug: print "ok2"
    #num.Draw("SAME PE0")
    draw_mr_bins([num], num.GetMinimum(),num.GetMaximum(), combine_bins, keep, mrbins, r2bins)
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
    ratio.Draw("P")
    #den_stat_err.SetFillColor(1)
    #den_stat_err.SetFillStyle(3004)
    #den_stat_err.SetMarkerStyle(0)
    #den_stat_err.Draw("SAME E2")
    #ratio.Draw("SAME PE0")
    if debug: print "ok3"
    if (xmin==xmax):
        xmin = ratio.GetYaxis().GetXmin()
        xmax = ratio.GetYaxis().GetXmax()
    if debug: print "ok3"
    l = ROOT.TLine(xmin, yratio, xmax, yratio)
    l.SetLineWidth(2)
    #l.SetLineColor(2)
    l.SetLineStyle(2)
    #l.Draw()
    ROOT.gPad.Update()
    if add_labels: add_r2_labels(ratio, combine_bins, keep, mrbins, r2bins)
    #can2.Write()
    if debug: print "ok3"
    ok = 1
    return can2

def add_r2_labels(h, combine_bins, keep, 
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
        keep.append(lat)
    h.GetXaxis().SetTitle("R^{2}")
    h.GetXaxis().SetTitleOffset(2.2)

def add_bin_labels(h, combine_bins,
                   mrbins  = [ 800, 1000, 1200, 1600, 2000, 4000 ],
                   r2bins  = [ 0.08, 0.12, 0.16, 0.24, 0.4, 1.5  ],
                   offset = 2.2):
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
    h.GetXaxis().SetTitle("R^{2}")
    h.GetXaxis().SetTitleOffset(offset)

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
    ROOT.gPad.Update()
    # Check if the directory exists (if not create it first)
    dirname = os.path.dirname(plotname)
    if dirname != "" and not os.path.exists(dirname):
        os.makedirs(dirname)
    with suppress_stdout_stderr():
        can.SaveAs(plotname+".png")
        can.SaveAs(plotname+".pdf")
        can.SaveAs(plotname+".C")
    if write:
        if name != "":
            can.Write(name)
        else:
            can.Write()

# Show and run command with stdout on screen
icommand=0
def special_call(cmd, run=1, verbose=1):
    global icommand
    if verbose:
        if run:
            print("[%d]" % icommand),
        else:
            print("[dry]"),
        for i in xrange(len(cmd)): print cmd[i],
        print ""
    if run:
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
def logged_call(cmd, logfile, run=1):
    dirname = os.path.dirname(logfile)
    if dirname != "" and not os.path.exists(dirname):
        special_call(["mkdir", "-p", os.path.dirname(logfile)], run, 0)
    if run:
        ntry = 0
        while True:
            try:
                with open(logfile, "a") as log:
                    proc = subprocess.Popen(cmd, stdout=log, stderr=log, close_fds=True)
                    proc.wait()
            except:
                print "Could not write to disk (IOError), wait 10s and continue"
                time.sleep(10)
                ntry += 1
                if ntry == 20: sys.exit()
                continue
            break
    else:
        proc = subprocess.call(["echo", "[dry]"]+cmd+[">", logfile])

def add_cms_era(plot, approval, keep, energy=14, intlumi=3000, prefix = "", twoframe=False):
    ROOT.gPad.Update() # Forces to create frames
    if plot.InheritsFrom("TCanvas"):
        # Check if its a two frame plot
        if plot.GetListOfPrimitives().GetEntries()==2:
            if plot.GetListOfPrimitives().At(0).GetName() != "TFrame":
                f1 = plot.GetListOfPrimitives().At(0).GetListOfPrimitives().At(0)
                f2 = plot.GetListOfPrimitives().At(1).GetListOfPrimitives().At(0)
                if f1.GetName()=="TFrame" and f2.GetName()=="TFrame": twoframe = True
        if twoframe:
            # cd to first frame and find histo
            plot.GetListOfPrimitives().At(0).cd()
            h = plot.GetListOfPrimitives().At(0).GetListOfPrimitives().At(1)
            if h.InheritsFrom("TH1"): add_cms_era(h, approval, keep, energy, intlumi, prefix, twoframe)
        else:
            xmin = plot.GetListOfPrimitives().At(0).GetX1()
            xmax = plot.GetListOfPrimitives().At(0).GetX2()
            ymin = plot.GetListOfPrimitives().At(0).GetY1()
            ymax = plot.GetListOfPrimitives().At(0).GetY2()
            era_and_prelim_lat_(approval, xmin, xmax, ymin, ymax, keep, energy, intlumi, prefix, twoframe, False, plot)
    elif plot.InheritsFrom("TH2"):
        xmin = plot.GetXaxis().GetBinLowEdge(plot.GetXaxis().GetFirst())
        xmax = plot.GetXaxis().GetBinUpEdge(plot.GetXaxis().GetLast())
        ymin = plot.GetYaxis().GetBinLowEdge(plot.GetYaxis().GetFirst())
        ymax = plot.GetYaxis().GetBinUpEdge(plot.GetYaxis().GetLast())
        era_and_prelim_lat_(approval, xmin, xmax, ymin, ymax, keep, energy, intlumi, prefix, twoframe)
    elif plot.InheritsFrom("TH1"):
        xmin = plot.GetXaxis().GetBinLowEdge(plot.GetXaxis().GetFirst())
        xmax = plot.GetXaxis().GetBinUpEdge(plot.GetXaxis().GetLast())
        ymin = plot.GetMinimum()
        ymax = plot.GetMaximum()
        era_and_prelim_lat_(approval, xmin, xmax, ymin, ymax, keep, energy, intlumi, prefix, twoframe)

def era_and_prelim_lat_(approval, xmin, xmax, ymin, ymax, keep, energy=14, intlumi=35.9, prefix = "", twoframe=False, inside=False, can=0):
    app = approval/10
    scale = 1.0
    if twoframe: scale = 10.0/7.0
    if app:
        # Latex example: #font[22]{Times bold} and #font[12]{Times Italic}
        text = ""
        if app==1: text = "CMS #scale[0.7]{#font[52]{Work in progress}}"
        if app==2: text = "CMS #scale[0.7]{#font[52]{Preliminary}}"
        if app==3: text = "CMS"
        if app==4: text = "#scale[0.8]{CMS Simulation }#scale[0.6]{#font[52]{Work in progress}}"
        if app==5: text = "#scale[0.8]{CMS Simulation }#scale[0.6]{#font[52]{Preliminary}}"
        if app==6:
            text = "CMS Phase-2 #font[52]{Projection}"
            scale = scale * 0.8
        if app==7:
            text = "CMS HE-LHC #font[52]{Projection}"
            scale = scale * 0.8
        #if app==7: text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}"
        if app==8: text = "CMS #scale[0.7]{#font[52]{Preliminary 2016}}"
        if app==9: text = "CMS #scale[0.7]{#font[52]{Preliminary 2018}}"
        x = xmin+(xmax-xmin)/20.0 if inside else xmin
        y = ymax-(ymax-ymin)/10.0 if inside else ymax+(ymax-ymin)/40.0
        if ymin>0:
            # Assume log scale
            if ymax/ymin>1000:
                if inside:
                    y = math.exp(math.log(ymax)-(math.log(ymax)-math.log(ymin))/10.0)
                else:
                    y = math.exp(math.log(ymax)+(math.log(ymax)-math.log(ymin))/40.0)
        cms_lat = ROOT.TLatex(x, y, text) 
        cms_lat.SetLineWidth(2)
        cms_lat.SetTextSize(cms_lat.GetTextSize()*scale)
        cms_lat.Draw()
        keep.append(cms_lat)
    era = approval%10
    if era:
        text = ""
        if intlumi<1000:
            totintlumi = ("%.1f" % intlumi)+" fb^{-1}"
        else:
            totintlumi = ("%.0f" % (intlumi/1000.0))+" ab^{-1}"
        if era==1: text = "#sqrt{s}=7 TeV"
        if era==2: text = "#sqrt{s}=8 TeV"
        if era==3: text = "#sqrt{s}=13 TeV"
        if era==4: text = "Run 2, #sqrt{s}=13 TeV"
        if era==5: text = "#scale[0.9]{"+prefix+totintlumi+" ("+str(energy)+" TeV)}"
        if era==6: text = "#scale[0.9]{"+prefix+totintlumi+" ("+str(energy)+" TeV)}"
        if era==7: text = "#scale[0.9]{"+prefix+totintlumi+" ("+str(energy)+" TeV)}"
        #if era==6: text = "#scale[0.65]{W ana, 35.9 fb^{-1} (13 TeV)}"
        #if era==7: text = "#scale[0.65]{Top ana, 35.9 fb^{-1} (13 TeV)}"
        y = ymax+(ymax-ymin)/25.0
        if ymin>0:
            # Assume log scale
            if ymax/ymin>1000:
                y = math.exp(math.log(ymax)+(math.log(ymax)-math.log(ymin))/25.0)
        era_lat = ROOT.TLatex(xmax, y, text)
        era_lat.SetTextAlign(32)
        era_lat.SetTextSize(0.04*scale)
        era_lat.SetTextFont(42)
        era_lat.SetLineWidth(2)
        era_lat.Draw()
        keep.append(era_lat)
        # Add W/Top box to the legend label
        if can and (era==6 or era==7):
            # find legend
            legend = 0
            for i in range(can.GetListOfPrimitives().GetEntries()):
                prim = can.GetListOfPrimitives().At(i)
                if prim.GetTitle().startswith("Legend"):
                    legend = prim
                    break
            if legend:
                header = legend.GetHeader()
                if not "category" in header:
                    if len(header): header += ", "
                    header += "W category" if era==6 else "Top category"
                    legend.SetHeader(header)
    ROOT.gPad.Update()
