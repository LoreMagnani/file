#cerca file chimati output_<indice>.root nella stessa cartella (e file lhe nella cartella padre se è attiva la rispettiva parte di codice)
#crea file chiamati histo_<nome variabile>.root contenenti un istogramma per ogni peso definito + un istogramma histo_Data vuoto per combine
#in base alle parti di codice attive crea istogrammi mirati su una certa variabile chiedendo xmin e xmax, oppure crea istogrammi in modo automatizzato 
#una porzione di codice è utile per creare istogrammi velocemente e scegliere i valori su cui ottimizzare senza perdere molto tempo
import sys
import ROOT
import argparse
import os
import awkward as ak
import pylhe
import vector
vector.register_awkward()
import numpy as np
import matplotlib.pyplot as plt
import hist
sys.path.insert(0, '/..')
ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True) 


#ottiene i percorsi assoluti dei file root e li salva in un array
folder_path = os.getcwd()
prefix = "output"
extension = ".root"
unsorted_files = [file for file in os.listdir(folder_path) if file.startswith(prefix) and file.endswith(extension)]
unsorted_files = np.array(unsorted_files)
def custom_sort(file_name):
    return int(file_name.split('_')[1].split('.')[0])
files = sorted(unsorted_files, key=custom_sort)
files = list(map(lambda k: folder_path + "/" + k, files))

n_files = len(files)

#creo dataframe
df = ROOT.RDataFrame("tree", files)

#definizione dei pesi
weights = {
        "w_sm": f"(weights[0])",
        "w_lin_cWtil": f"(0.5*(weights[2] - weights[1]))",
        "w_quad_cWtil": f"(-w_sm + 0.5*(weights[2] + weights[1]))",
        "w_sm_lin_quad_cWtil": f"(weights[2])",

        "w_lin_cHGtil": f"(0.5*(weights[4] - weights[3]))",
        "w_quad_cHGtil": f"(-w_sm + 0.5*(weights[4] + weights[3]))",
        "w_sm_lin_quad_cHGtil": f"(weights[4])",

        "w_lin_cHWtil": f"(0.5*(weights[6] - weights[5]))",
        "w_quad_cHWtil": f"(-w_sm + 0.5*(weights[6] + weights[5]))",
        "w_sm_lin_quad_cHWtil": f"(weights[6])",

        "w_lin_cHBtil": f"((0.5*(weights[8] - weights[7])))",
        "w_quad_cHBtil": f"((-w_sm + 0.5*(weights[8] + weights[7])))",
        "w_sm_lin_quad_cHBtil": f"((weights[8]))",

        "w_lin_cHWBtil": f"((0.5*(weights[10] - weights[9])))",
        "w_quad_cHWBtil": f"((-w_sm + 0.5*(weights[10] + weights[9])))",
        "w_sm_lin_quad_cHWBtil": f"((weights[10]))",

        "w_lin_cbWIm": f"((0.5*(weights[12] - weights[11])))",
        "w_quad_cbWIm": f"((-w_sm + 0.5*(weights[12] + weights[11])))",
        "w_sm_lin_quad_cbWIm": f"((weights[12]))",

        "w_lin_cbBIm": f"((0.5*(weights[14] - weights[13])))",
        "w_quad_cbBIm": f"((-w_sm + 0.5*(weights[14] + weights[13])))",
        "w_sm_lin_quad_cbBIm": f"((weights[14]))",

        "w_lin_clebQIm": f"((0.5*(weights[16] - weights[15])))",
        "w_quad_clebQIm": f"((-w_sm + 0.5*(weights[16] + weights[15])))",
        "w_sm_lin_quad_clebQIm": f"((weights[16]))"
        }    

#questo riempe il dataframe con array di pesi (non so a che servano)
for key in weights.keys():
    df = df.Define(key, weights[key])
    print(df.AsNumpy(columns=[key])[key])


histos = []

print(df.AsNumpy(columns=["nom_weights"])["nom_weights"])


#trovo somma pesi nominali. necessita di file lhe in ../
lhe_path = os.getcwd() + "/.."
prelhe = "1"
extlhe = ".lhe"
unsorted_fis = [fi for fi in os.listdir(lhe_path) if fi.startswith(prelhe) and fi.endswith(extlhe)]
unsorted_fis = np.array(unsorted_fis)
def custom_sort_lhe(fi_name):
    return int(fi_name.split('_')[1].split('.')[0])
fis = sorted(unsorted_fis, key=custom_sort_lhe)
fis = list(map(lambda k: lhe_path + "/" + k, fis))
tot_nom_weights = 0
for index in range(len(fis)):
    fi = fis[index]
    print(fi , index)
    xs = pylhe.read_lhe_init(fi)['procInfo'][0]['xSection']
    weight_names = ['rwgt_' + str(i) for i in range(1, 46)]
    arr = pylhe.to_awkward(pylhe.read_lhe_with_attributes(fi) , weight_names)
    nom_weights = (arr.eventinfo.weight)
    tot_nom_weights += ak.sum(nom_weights)

#crea istogrammi in modo automatizzato
vars = ["mll", "deltaPhill", "deltaEtall", "ptl1", "ptl2", "mJJ", "deltaPhiJJ", "deltaEtaJJ", "ptJ1", "ptJ2"]
xmin = [70 , -np.pi , -4 , 0 , 0 , 0  , -np.pi , -9 , 0 , 0]
xmax = [110 , np.pi , 4 , 900 , 400 , 5000 , np.pi , 9 , 1100 , 1000]
lum = 100
for w in weights.keys():
    #counts = df.Count().GetValue()
    counts = tot_nom_weights
    weight = f'100 * ({w} * nom_weights / {counts})'
    #weight = f'{lum} * 1000.0 * ' + w
    df = df.Redefine(w, weight)
for j in range(len(vars)):
    print(vars[j] + '\n')
    for w in weights.keys():
        print('\t' + w + '\n')
        print(df.AsNumpy(columns=[w])[w])
        histos.append(df.Histo1D(('histo_' + '_'.join(w.split('_')[1:]), "", 40, xmin[j], xmax[j]), vars[j], w))

    f = ROOT.TFile('hist_' + vars[j] + '.root' , 'RECREATE')
    f.cd()  

    savedData = False
    for histo in histos:
        histo.Write()
        if not savedData:
            h = histo.Clone()
            for i in range(0, h.GetNbinsX() +1):
                h.SetBinContent(i, 0)
            h.SetName('histo_Data')
            h.Write()
            savedData = True

    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)

    for i in range(25):
        canvas.Clear()
        histos[-i].Draw()
        canvas.SaveAs('histo_' + vars[j] + '_' + '_'.join(histos[-i].GetName().split('_')[1:]) + ".png")

    f.Close()
       
    print(vars[j] + ' fatto')

"""
#qua uso parser per riconoscere argomenti per fare solo alcuni istogrammi
def defaultParser():

    parser = argparse.ArgumentParser()

    parser.add_argument("-v",
                        type = str,
                        )

    parser.add_argument("-xmin",
                        type = int,
                        )

    parser.add_argument("-xmax",
                        type = int,
                        )
    return parser

parser = defaultParser()
args = parser.parse_args()

variable = args.v
lum = 100
xmax = args.xmax
xmin = args.xmin

if xmax == 3:
    xmax = np.pi
if xmin == -3:
    xmin = -np.pi

#questo crea istogrammi mirati sulla variabile 
for w in weights.keys():
        #weight = f'{lum} * 1000.0 * ' + w
        weight = f'{lum} * 1000.0 * ({w} * nom_weights)'
        df = df.Redefine(w, weight)
        histos.append(df.Histo1D(('histo_' + '_'.join(w.split('_')[1:]), "", 40, int(xmin), int(xmax)), variable, w))

f = ROOT.TFile('hist_' + variable + '.root' , 'RECREATE')
f.cd()   

savedData = False
for histo in histos:
    histo.Write()
    if not savedData:
        h = histo.Clone()
        for i in range(0, h.GetNbinsX() +1):
            h.SetBinContent(i, 0)
        h.SetName('histo_Data')
        h.Write()
        savedData = True

canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)

for histo in histos:
    canvas.Clear()
    histo.Draw()
    canvas.SaveAs(histo.GetName() + '_' + variable + ".png")

f.Close()
"""

"""
#righe utili per produrre velocemente istogrammi e cercare velocemente valori di xmin e xmax per automatizzare la produzione
pi = np.pi
for w in weights.keys():
    weight = f'{lum} * 1000.0 * ' + w
    df = df.Redefine(w, weight)
    histos.append(df.Histo1D(('histo_mll' + '_' + '_'.join(w.split('_')[1:]), "", 40, 70, 110), "mll", w))
    histos.append(df.Histo1D(('histo_deltaPhill' + '_' + '_'.join(w.split('_')[1:]), "", 40, -pi, pi), "deltaPhill", w))
    histos.append(df.Histo1D(('histo_deltaEtall' + '_' + '_'.join(w.split('_')[1:]), "", 40, -4, 4), "deltaEtall", w))
    histos.append(df.Histo1D(('histo_ptl1' + '_' + '_'.join(w.split('_')[1:]), "", 40, 0, 900), "ptl1", w))
    #histos.append(df.Histo1D(('histo_phil1' + '_' + '_'.join(w.split('_')[1:]), "", 40, -pi, pi), "phil1", w))
    #histos.append(df.Histo1D(('histo_etal1' + '_' + '_'.join(w.split('_')[1:]), "", 40, -pi, pi), "etal1", w))
    histos.append(df.Histo1D(('histo_ptl2' + '_' + '_'.join(w.split('_')[1:]), "", 40, 0, 400), "ptl2", w))
    #histos.append(df.Histo1D(('histo_phil2' + '_' + '_'.join(w.split('_')[1:]), "", 40, -pi, pi), "phil2", w))
    #histos.append(df.Histo1D(('histo_etal2' + '_' + '_'.join(w.split('_')[1:]), "", 40, -pi, pi), "etal2", w))
    histos.append(df.Histo1D(('histo_mJJ' + '_' + '_'.join(w.split('_')[1:]), "", 40, 0, 5000), "mJJ", w))
    histos.append(df.Histo1D(('histo_deltaPhiJJ' + '_' + '_'.join(w.split('_')[1:]), "", 40, -pi, pi), "deltaPhiJJ", w))
    histos.append(df.Histo1D(('histo_deltaEtaJJ' + '_' + '_'.join(w.split('_')[1:]), "", 40, -9, 9), "deltaEtaJJ", w))
    histos.append(df.Histo1D(('histo_ptJ1' + '_' + '_'.join(w.split('_')[1:]), "", 40, 0, 1100), "ptJ1", w))
    #histos.append(df.Histo1D(('histo_phiJ1' + '_' + '_'.join(w.split('_')[1:]), "", 40, -pi, pi), "phiJ1", w))
    #histos.append(df.Histo1D(('histo_etaJ1' + '_' + '_'.join(w.split('_')[1:]), "", 40, -5, 5), "etaJ1", w))
    histos.append(df.Histo1D(('histo_ptJ2' + '_' + '_'.join(w.split('_')[1:]), "", 40, 0, 1000), "ptJ2", w))
    #histos.append(df.Histo1D(('histo_phiJ2' + '_' + '_'.join(w.split('_')[1:]), "", 40, -pi, pi), "phiJ2", w))
    #histos.append(df.Histo1D(('histo_etaJ2' + '_' + '_'.join(w.split('_')[1:]), "", 40, -5, 5), "etaJ2", w))
"""





