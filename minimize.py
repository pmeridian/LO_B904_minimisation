import ROOT as R

rgen=R.TRandom3(100)

def smeared_scaled(tree,hFit,norm,scale,smear):
    hOut=hFit.Clone('hOut')
    hOut.Reset()
    for j in range(5):
        for e in tree:
            hOut.Fill(e.charge_sig[0]*scale+rgen.Gaus(0,smear),norm)
#    hOut.Sumw2()
    hOut.Scale(1/5.)
    hOut.Smooth()
    return hOut

def chi2_func(norm,scale,smear):
    chi2=0
    hOut=smeared_scaled(fft_digi_ori,h_2,norm,scale,smear)
    for i in range(min_bin,max_bin):
        chi2+=(h_2.GetBinContent(i)-hOut.GetBinContent(i))**2/(h_2.GetBinError(i)**2+(hOut.GetBinError(i)*5)**2)
    return chi2

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--run_ori', default='410_HPK_CH08_OV1.0_T18_run125')
parser.add_argument('--run_target', default='479_HPK_CH08_OV1.0_Tm40_run166')
parser.add_argument('--norm', type=float, default=5.8)
parser.add_argument('--scale', type=float, default=0.765)
parser.add_argument('--smear', type=float, default=230)

args = parser.parse_args()

run_ori=args.run_ori
f1=R.TFile('data/h4Reco_%s.root'%run_ori)

fft_digi_ori=f1.Get("fft_digi")
h_ori=R.TH1F('h_ori','h_ori',200,0,10000)

for e in fft_digi_ori:
    h_ori.Fill(e.charge_sig[0])

ps_ori=R.TSpectrum(6)
n_peaks_ori=ps_ori.Search(h_ori,3,"goff",0.8)
if(n_peaks_ori>0):
    print(n_peaks_ori,ps_ori.GetPositionX()[n_peaks_ori-1])

run_target=args.run_target
f2=R.TFile('data/h4Reco_%s.root'%run_target)

fft_digi_2=f2.Get("fft_digi")
fft_digi_2.SetBranchStatus("*",0)
fft_digi_2.SetBranchStatus("charge_sig",1)

h_2=R.TH1F('h_2','h_2',200,0,10000)
for e in fft_digi_2:
    h_2.Fill(e.charge_sig[0])

ps=R.TSpectrum(6)
n_peaks=ps.Search(h_2,3,"goff",0.8)
if(n_peaks>0):
    print(n_peaks,ps.GetPositionX()[n_peaks-1])

import numpy as np

ranges={
    'norm':{'ini':args.norm,'step':0.1},
    'scale':{'ini':args.scale,'step':0.005},
    'smear':{'ini':args.smear,'step':5}
}

minG={}
minX=max([.7*ps.GetPositionX()[n_peaks-1],1000])
maxX=1.7*ps.GetPositionX()[n_peaks-1]
min_bin=h_2.FindBin(minX)
max_bin=h_2.FindBin(maxX)
init_values={}

#init
print("Init values")
for ipar,(par,par_range) in enumerate(ranges.items()):
    minG[par]=R.TGraph()
    init_values[par]=par_range['ini']
print(init_values)

#scan
tval={}
min_values=[99999.]
tolerance=1
min_diff=99999.
iter=0

print("Start minimisation with tolerance %f"%tolerance)

watch=R.TStopwatch()
while (abs(min_diff)>tolerance):
    watch.Start()
    for ipar,(par,par_range) in enumerate(ranges.items()):
        #reinit
        minG[par].Clear()
        for ini_par in ranges:
            tval[ini_par]=init_values[ini_par]
        for ival,a in enumerate(np.arange(par_range['ini']-5*par_range['step'],par_range['ini']+5*par_range['step'],par_range['step'])):
            tval[par]=a
            val=chi2_func(tval['norm'],tval['scale'],tval['smear'])
            print(iter,par,ival,tval,val)
            minG[par].SetPoint(ival,a,val)
        minG[par].Fit('pol2')
    for par in ranges:
        init_values[par]=minG[par].GetFunction("pol2").GetMinimumX()
    min_values.append(chi2_func(init_values['norm'],init_values['scale'],init_values['smear']))
    min_diff=min_values[-1]-min_values[-2]
    watch.Stop()
    print("=====> ",iter,min_values[-1],min_diff,init_values)
    watch.Print()
    iter+=1

for ipar,(par,par_range) in enumerate(ranges.items()):
    #minG[par].Draw("A*")
    #minG[par].Fit("pol2")
    min_x=minG[par].GetFunction("pol2").GetMinimumX()
    min_val=minG[par].GetFunction("pol2").GetMinimum()
    delta_plus=minG[par].GetFunction("pol2").GetX(min_val+1,min_x,min_x+5*ranges[par]['step'])-min_val
    delta_minus=min_val-minG[par].GetFunction("pol2").GetX(min_val+1,min_x-5*ranges[par]['step'],min_x)
    delta=(delta_plus+delta_minus)*0.5
    ranges[par]['min']=min_x
    ranges[par]['err']=delta

h_fit_min=smeared_scaled(fft_digi_ori,h_2,ranges['norm']['min'],ranges['scale']['min'],ranges['smear']['min'])

c1=R.TCanvas('c1','c1',800,600)

#Draw fit result
R.gStyle.SetOptTitle(0)
R.gStyle.SetOptStat(0)

l=R.TLegend(0.1,0.7,0.4,0.86)
l.SetFillColorAlpha(0,0)
l.SetBorderSize(0)
l.SetTextSize(0.03)

h_2.SetStats(0)
h_2.Draw('PE')
h_2.SetMarkerStyle(20)
h_2.SetMarkerSize(1.2)
h_2.SetMarkerColor(R.kBlack)
h_2.SetLineColor(R.kBlack)
h_2.Draw('PE')
h_2.SetMaximum(h_2.GetMaximum()*1.6)
h_2.GetXaxis().SetRangeUser(minX,maxX)
h_2.GetXaxis().SetTitle('ADC')
l.AddEntry(h_2,'Target','P')

h_fit_min.Draw('PSAME')
h_fit_min.SetMarkerColor(R.kRed)
h_fit_min.SetMarkerStyle(24)
h_fit_min.SetMarkerSize(1.2)
#h_fit_min.SetLineColor(R.kRed)
l.AddEntry(h_fit_min,'Fit (Scaled + DCR Noise)','P')


#h_2.SetMaximum(4000)
h_2.GetXaxis().SetRangeUser(0.3*ps.GetPositionX()[n_peaks-1],2.2*ps.GetPositionX()[n_peaks-1])
#h_2.GetXaxis().SetTitle('ADC')
h_ori_c=h_ori.Clone('h_ori_c')
h_ori_c.Scale(h_2.GetEntries()/h_ori.GetEntries())
h_ori_c.Draw('SAME')
l.AddEntry(h_ori,'Original','PL')

t=R.TLatex()
t.SetTextSize(0.035)

results=R.TGraphErrors()
for ipar,(par,r) in enumerate(ranges.items()):
    t.DrawLatexNDC(0.6,0.825-ipar*0.04,"%s = %.3f #pm %.3f"%(par,r['min'],r['err']))
    results.SetPoint(ipar,ipar,r['min'])
    results.SetPointError(ipar,0,r['err'])
t.DrawLatexNDC(0.12,0.91,"%s"%(run_target))    
l.Draw()
c1.Draw()
for ext in ['.png','.pdf']:
    c1.SaveAs('results/fit_%s%s'%(run_target,ext))

fResults=R.TFile('results/fit_%s.root'%run_target,'RECREATE')
h_2.Write()
h_ori.Write()
h_fit_min.Write()
results.Write('fitResults')
fResults.Close()
