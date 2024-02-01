import sys
sys.path.insert(0, '/..')
import ROOT
import pylhe
import awkward as ak
import numpy as np
import os
import uproot
import vector
vector.register_awkward()

def main():
    #ottiene i percorsi assoluti dei file lhe e li salva in un array
    folder_path = os.getcwd()
    prefix = sys.argv[1]
    extension = ".lhe"
    unsorted_files = [file for file in os.listdir(folder_path) if file.startswith(prefix) and file.endswith(extension)]
    unsorted_files = np.array(unsorted_files)
    def custom_sort(file_name):
        return int(file_name.split('_')[1].split('.')[0])
    files = sorted(unsorted_files, key=custom_sort)
    files = list(map(lambda k: folder_path + "/" + k, files))

    #questo array contiene le sezioni d'urto di ogni file
    xs = []

    #itero su ogni file contenuto in files
    for index in range(len(files)):
        file = files[index]
 
        print(file, index)

        #salvo la sezione d'urto
        xs.append(pylhe.read_lhe_init(file)['procInfo'][0]['xSection'])

        #trova nomi dei pesi
        weight_names = ['rwgt_' + str(i) for i in range(1, 46)]

        # Lettura dei dati LHE
        arr = pylhe.to_awkward(pylhe.read_lhe_with_attributes(file) , weight_names)

        mll = (arr.particles.vector[:, -4]+arr.particles.vector[:, -3]).M
        deltaPhill = (arr.particles.vector[:, -4]).deltaphi(arr.particles.vector[:, -3])
        deltaEtall = (arr.particles.vector[:, -4]).deltaeta(arr.particles.vector[:, -3])

        firstIsFirst_l = np.array(arr.particles.vector[:, -4].pt > arr.particles.vector[:, -3].pt)

        ptl1  = np.array(arr.particles.vector[:, -4].pt)
        phil1 = np.array(arr.particles.vector[:, -4].phi)
        etal1 = np.array(arr.particles.vector[:, -4].eta)
        ptl1 [~firstIsFirst_l] = np.array(arr.particles.vector[:, -3].pt) [~firstIsFirst_l]
        phil1[~firstIsFirst_l] = np.array(arr.particles.vector[:, -3].phi)[~firstIsFirst_l]
        etal1[~firstIsFirst_l] = np.array(arr.particles.vector[:, -3].eta)[~firstIsFirst_l]

        ptl2  = np.array(arr.particles.vector[:, -3].pt)
        phil2 = np.array(arr.particles.vector[:, -3].phi)
        etal2 = np.array(arr.particles.vector[:, -3].eta)
        ptl2 [~firstIsFirst_l] = np.array(arr.particles.vector[:, -4].pt) [~firstIsFirst_l]
        phil2[~firstIsFirst_l] = np.array(arr.particles.vector[:, -4].phi)[~firstIsFirst_l]
        etal2[~firstIsFirst_l] = np.array(arr.particles.vector[:, -4].eta)[~firstIsFirst_l]

        #-------------------Jet variables-------------------------------

        mJJ = (arr.particles.vector[:, -2]+arr.particles.vector[:, -1]).M

        firstIsFirst_J = np.array(arr.particles.vector[:, -2].pt > arr.particles.vector[:, -1].pt)

        deltaPhiJJ = (arr.particles.vector[:, -2]).deltaphi(arr.particles.vector[:, -1])
        deltaEtaJJ = (arr.particles.vector[:, -2]).deltaeta(arr.particles.vector[:, -1])

        ptJ1 = np.array( arr.particles.vector[:, -2].pt)
        phiJ1  =np.array( arr.particles.vector[:, -2].phi)
        etaJ1 =np.array( arr.particles.vector[:, -2].eta)
        ptJ1 [~firstIsFirst_J] =np.array( arr.particles.vector[:, -1].pt) [~firstIsFirst_J]
        phiJ1 [~firstIsFirst_J] = np.array(arr.particles.vector[:, -1].phi) [~firstIsFirst_J]
        etaJ1 [~firstIsFirst_J] = np.array(arr.particles.vector[:, -1].eta) [~firstIsFirst_J]

        ptJ2 = np.array( arr.particles.vector[:, -1].pt)
        phiJ2  =np.array( arr.particles.vector[:, -1].phi)
        etaJ2 = np.array(arr.particles.vector[:, -1].eta)
        ptJ2 [~firstIsFirst_J] =np.array( arr.particles.vector[:, -2].pt) [~firstIsFirst_J]
        phiJ2 [~firstIsFirst_J] = np.array(arr.particles.vector[:, -2].phi) [~firstIsFirst_J]
        etaJ2 [~firstIsFirst_J] = np.array(arr.particles.vector[:, -2].eta) [~firstIsFirst_J]

        weights = arr.weights.values
        print(weights)

        weights = np.array(arr.weights.values.tolist())  # Converti la lista di liste in un array numpy
        weights = np.pad(weights, ((0, 0), (0, 45 - len(weights[0]))), 'constant', constant_values=0)
        print(weights)

        #calcolo peso nominale per tutti gli eventi e controllo la somma restituisca la sezione d'urto
        nom_weights = arr.eventinfo.weight * xs[index] * 1000
        print("cross-section del file: " + str(xs[index]))
        print("somma pesi nominali: " + str(ak.sum(nom_weights / ak.sum(arr.eventinfo.weight))))

        cuts = ((ptl1 > 25) & (ptl2 > 15)) \
    & (((arr.particles[:, -3].id*arr.particles[:, -4].id) == -11*11) \
    | ((arr.particles[:, -3].id*arr.particles[:, -4].id) == -13*13)) \
    & ak.all(abs(arr.particles[:, [-3, -4]].vector.eta) < 2.5, axis=1) \
    & (abs(mll - 91) < 15) \
    & (mJJ>200) \
    & (ptJ1>30) & (ptJ2>30) \
    & ak.all(abs(arr.particles[:, [-2, -1]].vector.eta) < 4.7, axis=1) 


        rootfile = uproot.recreate('output_' + str(index) + '.root')

        rootfile['tree'] = {
            'mll': ak.fill_none(mll[cuts], 0),
            'deltaPhill': ak.fill_none(deltaPhill[cuts], 0),
            'deltaEtall': ak.fill_none(deltaEtall[cuts], 0),
            'ptl1': ak.fill_none(ptl1[cuts], 0),
            'phil1': ak.fill_none(phil1[cuts], 0),
            'etal1': ak.fill_none(etal1[cuts], 0),
            'ptl2': ak.fill_none(ptl2[cuts], 0),
            'phil2': ak.fill_none(phil2[cuts], 0),
            'etal1': ak.fill_none(etal2[cuts], 0),

            'mJJ': ak.fill_none(mJJ[cuts], 0),
            'deltaPhiJJ': ak.fill_none(deltaPhiJJ[cuts], 0),
            'deltaEtaJJ': ak.fill_none(deltaEtaJJ[cuts], 0),
            'ptJ1': ak.fill_none(ptJ1[cuts], 0),
            'phiJ1': ak.fill_none(phiJ1[cuts], 0),
            'etaJ2': ak.fill_none(etaJ1[cuts], 0),
            'ptJ2': ak.fill_none(ptJ2[cuts], 0),
            'phiJ2': ak.fill_none(phiJ2[cuts], 0),
            'etaJ1': ak.fill_none(etaJ2[cuts], 0), 
            'weights': weights[cuts],
            'nom_weights': ak.fill_none(nom_weights[cuts], 0),
        }

        rootfile.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Funzionamento: Python3 an.py -nome_dei_file-")
        sys.exit(1)

    main()
