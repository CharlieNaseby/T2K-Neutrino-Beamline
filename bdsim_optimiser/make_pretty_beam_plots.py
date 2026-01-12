import ROOT
import sys
import pandas as pd


def strip_whitespace(df):
    df.columns = df.columns.str.strip()
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    return df


def draw_quadrupole(ax, s_start, length, name):
    ax.cd()
#    print("drawing quadrupole:", name, "at", s_start, "length", length)
    rect = ROOT.TBox(s_start, -1, s_start + length, 1)
    ROOT.SetOwnership(rect, False)
    rect.SetFillColor(ROOT.kRed - 7)
#    rect.SetLineColor(ROOT.kBlue + 2)
    rect.Draw("same")
    label = ROOT.TLatex()
    label.SetTextAngle(-45)
    label.SetTextSize(0.04)
    label.DrawLatex(s_start + length / 2, -1.2, name)
    #ax.Update()

def draw_bend(ax, s_start, length, name):
    ax.cd()
#    print("drawing bend:", name, "at", s_start, "length", length)
    rect = ROOT.TBox(s_start, -1, s_start + length, 1)

    ROOT.SetOwnership(rect, False)
    rect.SetFillColor(ROOT.kBlue - 7)
#    rect.SetLineColor(ROOT.kRed + 2)
    rect.Draw("same")
    label = ROOT.TLatex()
    label.SetTextAngle(-45)
    label.SetTextSize(0.04)
    label.DrawLatex(s_start + length / 2, -1.2, name)
    #ax.Update()

def draw_ssem(ax, s_position, name):
    ax.cd()
#    print("drawing ssem:", name, "at", s_position)
    line = ROOT.TLine(s_position, -1, s_position, 1)

    ROOT.SetOwnership(line, False)
    line.SetLineColor(ROOT.kGreen + 2)
    line.SetLineWidth(2)
    line.Draw("same")
    label = ROOT.TLatex()
    #text at 45 degree angle
    label.SetTextAngle(-45)
    label.SetTextSize(0.04)
    label.DrawLatex(s_position, -1.2, name)
    #ax.Update()


def draw_graph(pad, chisq, title, runnum, xrange, graph, data=None):
    pad.cd()
    graph.SetLineColor(ROOT.kRed)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetRangeUser(xrange[0], xrange[1])
    graph.Draw()
    graph.SetTitle("MR Run "+str(runnum))
    if data:
        data.SetMarkerStyle(20)
        data.SetMarkerSize(0.8)
        data.Draw("same P*")
    chisq_text = ROOT.TPaveText(0.6, 0.7, 0.9, 0.9, "NDC")
    chisq_text.SetShadowColor(0)
    chisq_text.SetBorderSize(1)

    chisq_text.SetFillColor(0)
    chisq_text.SetFillStyle(0)
    
    chisq_text.AddText("Total Fit #chi^{2} = " + str(chisq))
    ROOT.SetOwnership(chisq_text, False)
    chisq_text.Draw("same")
    graph.GetYaxis().SetTitle(title)
    #pad.update()

def draw_diagram(pad, xrange, line):
    pad.cd()
    hist = ROOT.TH1F("hist_diag", "", 1000, xrange[0], xrange[1])
    hist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
    hist.GetXaxis().SetTitle("Position along beamline (m)")
    hist.Draw()
    hist.SetMinimum(-2)
    hist.SetMaximum(2)
    for element in line.itertuples():
        if element.type == 'quadrupole':
            draw_quadrupole(pad, element.s_start + float(element.mark) - element.polelength / 2, element.polelength, element.element)
        elif element.type == 'rbend' or element.type == 'sbend':
            draw_bend(pad, element.s_start + float(element.mark) - element.polelength / 2, element.polelength, element.element)
        elif element.type == 'ssem':
            draw_ssem(pad, element.s_start + float(element.mark), element.element)
    #pad.Update()
    

def draw(line, chisq, name, title, runnum, graph, data=None):
    canv = ROOT.TCanvas("canv_" + name, "Beam Optics", 800, 600)
    canv.cd()
    ROOT.gStyle.SetOptStat(0)
    graphpad = ROOT.TPad("graphpad_" + name, "Graph Pad", 0.0, 0.4, 1.0, 1.0)
    diagrampad = ROOT.TPad("diagrampad_" + name, "Diagram Pad", 0.0, 0.0, 1.0, 0.4)
    graphpad.SetBottomMargin(0.0)
    diagrampad.SetTopMargin(0.0)
    canv.cd()
    graphpad.Draw()
    graphpad.cd()


    xrange = [graph.GetXaxis().GetXmin(), graph.GetXaxis().GetXmax()]

    draw_graph(graphpad, chisq, title, runnum, xrange, graph, data)
    canv.cd()
    diagrampad.Draw()
    diagrampad.cd()
    draw_diagram(diagrampad, xrange, line)
    canv.Update()
    canv.SaveAs("pretty_plots/pretty_" + str(runnum) + "_" + name + ".pdf")
    canv.SaveAs("pretty_plots/pretty_" + str(runnum) + "_" + name + ".png")




if __name__ == "__main__":

    line = strip_whitespace(pd.read_csv("../fujii-san.csv", header=0, skipinitialspace=True)) #the csv containing beampipe properties
    line['length'] = line['length'].astype(float)
    #replace 'NA' in 'mark' column with 0
    line['mark'] = line['mark'].replace('NA', 0)
    line['mark'] = line['mark'].astype(float)
    line['polelength'] = line['polelength'].astype(float)
    line['length'] = line['length']*0.001  #convert from mm to m
    line['mark'] = line['mark']*0.001  #convert from mm to m
    line['polelength'] = line['polelength']*0.001  #convert from mm to m
    line['s_start'] = line['length'].shift().cumsum()
    beamline_length = line['length'].sum()




    if(len(sys.argv) != 2):
        print("Usage: python make_pretty_beam_plots.py <run_number>")
        sys.exit(1)

    fit_result = ROOT.TFile("fit_results_" + sys.argv[1] + ".root", "READ")
    chisq = fit_result.Get("chisq")
    chisq = chisq[0]


    fitted_beam = ROOT.TFile("fit_results_" + sys.argv[1] + "_plots.root", "READ")

    sim = ["opt_x", "opt_y", "opt_width_x", "opt_width_y"]
    data = ["data_xbds", "data_ybds", "data_xwbds", "data_ywbds"]
    names = ["x", "y", "width_x", "width_y"]
    titles = ["Beam Position X (mm)", "Beam Position Y (mm)", "Beam Width X (mm)", "Beam Width Y (mm)"]
    for i in range(len(sim)):
        graph = fitted_beam.Get(sim[i])
        datapoints = fitted_beam.Get(data[i])
        draw(line, chisq, names[i], titles[i], sys.argv[1], graph, datapoints)

    sim = ["opt_alpha_x", "opt_alpha_y", "opt_beta_x", "opt_beta_y"]
    names = ["alpha_x", "alpha_y", "beta_x", "beta_y"]
    titles = ["Alpha X", "Alpha Y", "Beta X (m)", "Beta Y (m)"]

    for i in range(len(sim)):
        graph = fitted_beam.Get(sim[i])
        draw(line, chisq, names[i], titles[i], sys.argv[1], graph)

