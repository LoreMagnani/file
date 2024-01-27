import uproot
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import mplhep as hep
import os

d = hep.style.CMS

output_path = os.getcwd()

# might change d 
# d['xtick.major.size'] = 20

plt.style.use([d, hep.style.firamath])


#f = uproot.open('output_hist.root')
effes = []
for nome_file in os.listdir(output_path):
    if nome_file.endswith(".root") and nome_file.startswith("hist_"):
        percorso_file = os.path.join(output_path, nome_file)
        if os.path.isfile(percorso_file):
            f = uproot.open(nome_file)
            effes.append(f)
        else:
            print(percorso_file , "non esiste")

        print("File: " , nome_file)

        opes = ["cWtil" , "cHGtil" , "cHWtil" , "cHBtil" , "cHWBtil" , "cbWIm" , "cbBIm" , "clebQIm"]
        for ope in opes:

            def make_error_band(bins, err):
                _x = np.array(list(zip(bins-1e-8, bins))).flatten()[1:-1]
                _err = np.array(list(zip(err, err))).flatten()
                return _x, _err

            def ratioFun(num, den, _min=1e-8, _max=10.0, tolerance=1e-12):
                ratio=np.ones(len(num))
                # check where the denominator is almost zero and num is not -> set ratio to _max
                mask1 = (den < tolerance) & (num > tolerance)
                # check where both the den and num are almost zero -> set ratio to _min (discard points)
                mask2 = (den < tolerance) & (num < tolerance)
                # check where the den is greater than zero -> set ratio to num/den
                mask3 = (den > tolerance)
                
                ratio[mask1] = _max
                ratio[mask2] = _min
                ratio[mask3] = num[mask3] / den[mask3]
                
                # check if ratio is out of bounds for min and max
                mask_min = ratio < _min
                ratio[mask_min] = _min

                mask_max = ratio > _max
                ratio[mask_max] = _max
                
                return ratio



            def get_darker_color(color):
                rgb = list(matplotlib.colors.to_rgba(color)[:-1])
                darker_factor = 4/5
                rgb[0] = (rgb[0] * darker_factor)
                rgb[1] = (rgb[1] * darker_factor)
                rgb[2] = (rgb[2] * darker_factor)
                return tuple(rgb)

            def get_shade_color(color):
                rgb = list(matplotlib.colors.to_rgba(color)[:-1])
                darker_factor = 1/5
                rgb[0] = ((1.0 - rgb[0]) * darker_factor) + rgb[0]
                rgb[1] = ((1.0 - rgb[1]) * darker_factor) + rgb[1]
                rgb[2] = ((1.0 - rgb[2]) * darker_factor) + rgb[2]
                return tuple(rgb)

            def nome_Latex(variableName):
                new_nome = ""

                # Aggiungi la parte iniziale
                if 'delta' in variableName:
                    new_nome += r'\Delta'
                if 'Eta' in variableName:
                    new_nome += r'\eta'
                elif 'Phi' in variableName:
                    new_nome += r'\phi'
                elif 'm' in variableName:
                    new_nome += r'm'
                elif 'pt' in variableName:
                    new_nome += r'p_{t'

                # Aggiungi la parte finale
                if 'JJ' in variableName:
                    new_nome += '_{jj}'
                elif 'll' in variableName:
                    new_nome += '_{ll}'
                elif 'J1' in variableName:
                    new_nome += '_{j1}}'
                elif 'J2' in variableName:
                    new_nome += '_{j2}}'
                elif 'l1' in variableName:
                    new_nome += '_{l1}}'
                elif 'l2' in variableName:
                    new_nome += '_{l2}}'

                return new_nome

            groupPlot = {}

            groupPlot['SM']  = {
                            'nameHR' : 'SM',
                            'isSignal' : 0,
                            'color': 'lightgray',    # kGray + 1
                            'samples'  : ['sm']
            }

            operator1 = ope
            operator1Value = 1
            operator2 = 'per ora no'
            operator2Value = 0.01
            groupPlot['Lin1']  = {
                            'nameHR' : 'Lin '+ operator1 + ' = ' + str(operator1Value),
                            'isSignal' : 0,
                            'color': 'purple',    # kGray + 1
                            'samples'  : ['lin_'+operator1]
            }

            groupPlot['Quad1']  = {
                            'nameHR' : 'Quad '+operator1 + ' = ' + str(operator1Value),
                            'isSignal' : 0,
                            'color': 'red',    # kGray + 1
                            'samples'  : ['quad_' + operator1]
            }

            groupPlot['Lin2']  = {
                            'nameHR' : 'Lin '+ operator2 + ' = ' + str(operator2Value),
                            'isSignal' : 0,
                            'color': 'pink',    # kGray + 1
                            'samples'  : ['lin_'+operator2]
            }

            groupPlot['Quad2']  = {
                            'nameHR' : 'Quad '+operator2 + ' = ' + str(operator2Value),
                            'isSignal' : 0,
                            'color': 'orange',    # kGray + 1
                            'samples'  : ['quad_' + operator2]
            }

            groupPlot['Mix']  = {
                            'nameHR' : 'Mix ' + operator1 + ' = ' + str(operator1Value) + ', ' + operator2 + ' = ' + str(operator2Value),
                            'isSignal' : 0,
                            'color': 'lightgreen',    # kGray + 1
                            'samples'  : ['mixed_' + operator1 + '_' + operator2]
            }

            """
            variableName = 'mll'
            variableNameLatex = r'm_{ll}'
            variableNameUnit = 'GeV'
            """
            variableName = nome_file[nome_file.find("_") + 1:nome_file.find(".")]
            print(variableName)
            variableNameLatex = nome_Latex(variableName)
            print(variableNameLatex)
            variableNameUnit = 'GeV'

            samples = {}
            bins = -1
            c_data = 0

            flow = False

            for histo in f.keys():
                if type(bins) == type(-1):
                    bins = f[histo].to_numpy(flow)[1]

                h = f[histo].to_numpy(flow)[0]
                e = np.power(f[histo].errors(flow), 2)
                found = False
                histoName = f[histo].name[len('histo_'):]

                for sampleName in groupPlot.keys():
                    if histoName in groupPlot[sampleName]['samples']:
                        if 'Lin1' in sampleName:
                            h = operator1Value * h
                            e = operator1Value**2 * e
                        elif 'Quad1' in sampleName:
                            h = operator1Value**2 * h
                            e = operator1Value**4 * e
                        elif 'Lin2' in sampleName:
                            h = operator2Value * h
                            e = operator2Value**2 * e
                        elif 'Quad2' in sampleName:
                            h = operator2Value**2 * h
                            e = operator2Value**4 * e
                        elif 'Mix' in sampleName:
                            h = operator1Value * operator2Value * h
                            e = operator1Value**2 * operator2Value**2 * e
                        if sampleName in samples.keys():
                            samples[sampleName]['h'] += h
                            samples[sampleName]['e2']+= e
                        else:
                            samples[sampleName] = {
                                    'h': h,
                                    'e2': e,
                            }
                        found = True
                        break
                if not found:
                    print('Could not find', histoName, 'in groupPlot')



            #print(list(groupPlot.keys()))

            x0 = bins[0]
            x1 = bins[-1]
            print('\n\nx0:', x0, '\tx1:', x1, '\n\n')
            y0 = 1.0
            y1 = 10**6
            ratio0=0.1
            ratio1=3.0
            tolerance = 1e-8


            # nuis_err2_up = np.zeros(len(list(samples.values())[0]['h']))
            # nuis_err2_do = np.zeros(len(list(samples.values())[0]['h']))
            # sig_nuis_err2_up = np.zeros(len(list(samples.values())[0]['h']))
            # sig_nuis_err2_do = np.zeros(len(list(samples.values())[0]['h']))

            #thstack = np.zeros(len(list(samples.values())[0]['h']))
            thstack = []
            thstack_nuis_err2_up = []
            thstack_nuis_err2_do = []
            c_mcs_nominal = []

            #nuis_err2_do = np.zeros(len(c_mcs_nominal[0][1]))
            for sample in samples.keys():
                    thstack_nuis_err2_up.append((sample, samples[sample]['e2']))
                    thstack_nuis_err2_do.append((sample, samples[sample]['e2']))
                    thstack.append((sample, samples[sample]['h']))
                    c_mcs_nominal.append((sample, samples[sample]['h']))

            # stat_err     = np.sqrt(nuis_err2_up)
            # sig_stat_err = np.sqrt(sig_nuis_err2_up)

            # # FIXME fill nuisances
            # nuis_err_up = np.sqrt(nuis_err2_up)
            # nuis_err_do = np.sqrt(nuis_err2_do)
            # sig_nuis_err_up = np.sqrt(sig_nuis_err2_up)
            # sig_nuis_err_do = np.sqrt(sig_nuis_err2_do)

            #print(c_mcs_nominal[0][1])

            #c_mcs_nominal = sorted(c_mcs_nominal, key=lambda k: np.sum(k[1]))
            #print(list(groupPlot.keys()))

            c_mcs_nominal = sorted(c_mcs_nominal, key=lambda k: list(groupPlot.keys()).index(k[0]), reverse=True)
            thstack = list(map(lambda k: k[1], sorted(thstack, key=lambda k: list(groupPlot.keys()).index(k[0]), reverse=True)))
            thstack_nuis_err2_up = list(map(lambda k: k[1], sorted(thstack_nuis_err2_up, key=lambda k: list(groupPlot.keys()).index(k[0]), reverse=True)))
            thstack_nuis_err2_do = list(map(lambda k: k[1], sorted(thstack_nuis_err2_do, key=lambda k: list(groupPlot.keys()).index(k[0]), reverse=True)))

            # _names = sorted(names, key=lambda k: list(groupPlot.keys()).index(k), reverse=True)
            # sortingMap = [names.index(_name) for _name in _names]
            # print(sortingMap)

            # thstack = sorted(thstack, key=sortingMap.index, reverse=True)
            # thstack_nuis_err2_up = sorted(thstack_nuis_err2_up, key=lambda k: list(groupPlot.keys()).index(k[0]), reverse=True)
            # thstack_nuis_err2_do = sorted(thstack_nuis_err2_do, key=lambda k: list(groupPlot.keys()).index(k[0]), reverse=True)
            # print(c_mcs_nominal)

                    
            centers = (bins[:-1]+bins[1:])/2
            #err_bar_x = (centers[1]-centers[0])/2
            err_bar_x = 0


            fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3,1]}, figsize=(10,10), dpi=100)#dpi=100
            fig.tight_layout(pad=-0.5)

            hep.cms.label('Work in progress', data=False, ax=ax[0])


            thstack_base = np.zeros(len(thstack[0]))

            for i in range(len(c_mcs_nominal)-1, -1, -1):
                name = groupPlot[c_mcs_nominal[i][0]]['nameHR']
                color = groupPlot[c_mcs_nominal[i][0]]['color']
                thstack_base += c_mcs_nominal[i][1]
                integral = str(round(np.sum(c_mcs_nominal[i][1])*1e+2,2))+'1e-2'
                ax[0].hist(centers, weights=thstack_base, bins=bins, label=name + f' [{integral}]', color=color,
                        edgecolor=get_darker_color(color), histtype='step', fill=True, linewidth=2, zorder=(-100+i))


            for i in range(len(c_mcs_nominal)-1, -1, -1):
                name = groupPlot[c_mcs_nominal[i][0]]['nameHR']
                color = groupPlot[c_mcs_nominal[i][0]]['color']
                #thstack_base += c_mcs_nominal[i][1]
                print(ax[0].hist(centers, weights=c_mcs_nominal[i][1], bins=bins, color=color,
                        edgecolor=get_darker_color(color), histtype='step', fill=False, linewidth=2, zorder=1))
                top_err_band_up = c_mcs_nominal[i][1] + np.sqrt(thstack_nuis_err2_up[i])
                top_err_band_do = c_mcs_nominal[i][1] - np.sqrt(thstack_nuis_err2_do[i])
                x, top_err_band_up = make_error_band(bins, top_err_band_up)
                _, top_err_band_do = make_error_band(bins, top_err_band_do)
                ax[0].fill_between(x , top_err_band_do, top_err_band_up,  hatch='///', color=color, zorder=2, alpha=1, facecolor='none', linewidth=0.0)



            # # integral = str(round(np.sum(c_data),2))
            # # ax[0].errorbar(centers, c_data, err_bar_data, err_bar_x, fmt='ko', markersize=6, ecolor='black', label='DATA'+ f' [{integral}]', zorder=10)

            # # signalName = ''
            # # for h in groupPlot.keys():
            # #     if groupPlot[h].get('isSignal', 0) == 1:
            # #         signalName = h
            # #         break
            # # print(signalName)
            # # print(signalName == '')

            # # h_signal = list(filter(lambda k: k[0] == signalName, c_mcs_nominal))[0][1]
            # # includeSigInTopBand = False
            # # label = 'Stat + Syst'
            # # if includeSigInTopBand:
            # #     mc = sum(thstack) + sum(sig_thstack)
            # #     top_err_band_up = mc.copy() + np.sqrt( np.power(nuis_err_up, 2) + np.power(sig_nuis_err_up, 2) )
            # #     top_err_band_do = mc.copy() - np.sqrt( np.power(nuis_err_do, 2) + np.power(sig_nuis_err_do, 2) )
            # # else:
            # #     mc = sum(thstack) 
            # #     top_err_band_up = mc.copy() + nuis_err_up
            # #     top_err_band_do = mc.copy() - nuis_err_do
            # #     label += ' on Bkg.'
            # sm_lin_quad = sum(thstack)
            # top_err_band_up = sm_lin_quad.copy() + np.sqrt( sum(thstack_nuis_err2_up) )
            # top_err_band_do = sm_lin_quad.copy() - np.sqrt( sum(thstack_nuis_err2_do) )
            # # top_err_band_do = sm_lin_quad.copy() - np.sqrt( np.power(nuis_err_do, 2) + np.power(sig_nuis_err_do, 2) )
                
            # # print(top_err_band_do, mc-top_err_band_do)

            # label = 'Stat'
            # x, top_err_band_up = make_error_band(bins, top_err_band_up)
            # _, top_err_band_do = make_error_band(bins, top_err_band_do)
            # ax[0].fill_between(x , top_err_band_do, top_err_band_up,  hatch='///', color='darkgrey', zorder=9, alpha=1, facecolor='none', linewidth=0.0, label=label)


            # # mc = sum(thstack)

            # # # ratio=np.ones(len(mc))
            # # # mask = mc>tolerance
            # # # ratio[mask] = c_data[mask]/mc[mask]

            # # # mask = (mc<tolerance) & (c_data>tolerance)
            # # # ratio[mask] = 10

            # # ratio = ratioFun(c_data, mc)

            # # # # If outside the plot put on limit
            # # # mask = ratio > ratio1
            # # # ratio[mask] =  ratio1
            # # # mask = ratio < ratio0
            # # # ratio[mask] =  ratio0


            # # # Data points in ratio error

            # # # err_bar = np.ones(len(mc))
            # # # mask = mc > tolerance
            # # # err_bar[mask] = err_bar_data[mask] / mc[mask]

            # # err_bar = ratioFun(err_bar_data, mc)


            # # # err_band_up = np.ones(len(mc))
            # # # mask = mc>tolerance
            # # # err_band_up[mask] = nuis_err_up[mask]/mc[mask]
            # # err_band_up = ratioFun(nuis_err_up, mc)

            # # # err_band_do = np.ones(len(mc))
            # # # mask = mc>tolerance
            # # # err_band_do[mask] = nuis_err_do[mask]/mc[mask]
            # # err_band_do = ratioFun(nuis_err_do, mc)

            # # x, err_band_up = make_error_band(bins, err_band_up)
            # # _, err_band_do = make_error_band(bins, err_band_do)

            # # ax[1].fill_between(x , 1 - err_band_do, 1 + err_band_up, alpha=0.2, hatch='///', color='grey', label='Stat + Syst')

            # # ax[1].errorbar(centers, ratio, err_bar, err_bar_x, fmt='ko', markersize=6, ecolor='black')

            # # ax[1].plot(bins, np.ones(len(bins)), '-', color='black')


            # # # stat err band

            # # # err_band_up = np.ones(len(mc))
            # # # mask = mc>tolerance
            # # # err_band_up[mask] = stat_err[mask]/mc[mask]

            # # err_band_up = ratioFun(stat_err, mc)
            # # err_band_do = err_band_up.copy()

            # # x, err_band_up = make_error_band(bins, err_band_up)
            # # _, err_band_do = make_error_band(bins, err_band_do)

            # # ax[1].fill_between(x , 1 - err_band_do, 1 + err_band_up, hatch='\\\\', color='red', facecolor='none', label='Stat only')


            # # Data / Pred.
            # # mc = sum(thstack)  + sum(sig_thstack)
            sm = c_mcs_nominal[2][1]
            sm_lin = sum([c_mcs_nominal[2][1], c_mcs_nominal[1][1]])
            sm_lin_quad = sum([c_mcs_nominal[2][1] , c_mcs_nominal[1][1] , c_mcs_nominal[0][1]])
            # print(c_mcs_nominal[2][0] )


            # # err_band_up = np.ones(len(mc))

            # # sm = thstack[2]
            # # sm_lin = sum(thstack[1:3])
            # # sm_lin_quad = sum(thstack)
            # # ratio=np.ones(len(mc))
            # # mask = mc>tolerance
            # # ratio[mask] = c_data[mask]/mc[mask]

            # # mask = (mc<tolerance) & (c_data>tolerance)
            # # ratio[mask] = 10

            ratio = ratioFun(sm_lin, sm)
            ratio2 = ratioFun(sm_lin_quad, sm)

            # print(ratio)
            # print(ratio2)

            # # mask = ratio > ratio1
            # # ratio[mask] =  ratio1
            # # mask = ratio < ratio0
            # # ratio[mask] =  ratio0


            # # Data points in ratio error
            # # err_bar = ratioFun(sm, sm_lin)
            # # err_bar2 = ratioFun(sm, sm_lin_quad)
            err_bar = ratioFun(np.sqrt(sum(thstack_nuis_err2_do[1:3])),  sm)
            # print(err_bar)
            err_bar2 = ratioFun(np.sqrt(sum(thstack_nuis_err2_do[0:3])),  sm)
            # # err_bar2 = np.zeros(len(ratio))


            # # err_band_up = ratioFun(np.sqrt(sum(thstack_nuis_err2_up[1:3])), sm_lin)
            # # err_band_do = ratioFun(np.sqrt(sum(thstack_nuis_err2_do[1:3])), sm_lin)

            # # x, err_band_up = make_error_band(bins, err_band_up)
            # # _, err_band_do = make_error_band(bins, err_band_do)

            err_band_up2 = ratioFun(np.sqrt(sum(thstack_nuis_err2_up[2:3])), sm)
            err_band_do2 = ratioFun(np.sqrt(sum(thstack_nuis_err2_do[2:3])), sm)

            _, err_band_up2 = make_error_band(bins, err_band_up2)
            _, err_band_do2 = make_error_band(bins, err_band_do2)

            ax[1].errorbar(centers, ratio, err_bar, err_bar_x, fmt='ko', markersize=6, markerfacecolor='purple', label='(SM + LIN) / SM')
            ax[1].errorbar(centers, ratio2, err_bar2, err_bar_x, fmt='ko', markersize=6, markerfacecolor='red', label='BSM / SM')
            ax[1].fill_between(x , 1 - err_band_do2, 1 + err_band_up2, alpha=0.4, hatch='///', color='lightgray', label='Stat SM')
            ax[1].plot(bins, np.ones(len(bins)), '-', color='black')

            # ax[1].legend(ncol=3, fontsize='x-small', loc='lower left')
            # # nuis_sym = (np.sqrt(np.power(nuis_err_up, 2) + np.power(sig_nuis_err_up, 2)) + np.sqrt(np.power(nuis_err_do, 2) + np.power(sig_nuis_err_do, 2)))/2
            # # chi2 = np.sum(np.power(ratio-1, 2) / ( np.power(nuis_sym/mc, 2) + np.power(err_bar, 2)))
            # # #chi2_v = (np.power(ratio-1, 2) / ( np.power(nuis_sym/mc, 2) + np.power(err_bar, 2)))
            # # #print(chi2_v)
            # # ax[2].plot([], [],' ', label=r'$ \chi ^2 _0 $=' + str(round(chi2/len(ratio),2)))
            # # print('\n\nChi2:', chi2/len(ratio))


            # # # stat only
            # # # err_band_up = np.ones(len(mc))
            # # # mask = mc>tolerance
            # # # err_band_up[mask] =  / mc[mask]
            # # err_band_up = ratioFun(np.sqrt(np.power(stat_err, 2) + np.power(sig_stat_err, 2)), mc)
            # # err_band_do = err_band_up.copy()


            # x, err_band_up = make_error_band(bins, err_band_up)
            # _, err_band_do = make_error_band(bins, err_band_do)
            # # ax[2].fill_between(x , 1 - err_band_do, 1 + err_band_up, hatch='\\\\', color='red', facecolor='none', label='Stat only')
            # ax[2].fill_between(x , 1 - err_band_do, 1 + err_band_up, hatch='\\\\', color='red', facecolor='none')

            # #ax[2].legend()
            handles, labels = plt.gca().get_legend_handles_labels()
            # print(labels)
            order = [1,2,0]
            ax[1].legend([handles[idx] for idx in order],[labels[idx] for idx in order], ncol=3, fontsize='x-small')
            #ax[1].legend(ncol=3, fontsize='x-small')

            y0 = np.min(thstack[0])
            #y1 =np.max(sm_lin_quad)

            eps = (y1-y0)/2
            y0 -= eps
            # y0 = max(y0, 1)
            y1 += eps*5

            #y0 = max(np.min(thstack), 1)
            #y0 = 1e+1
            print('\n\nMin on the y:', y0)
            #y1 = np.max(sm_lin_quad) *1e+2

            # eps = (y1-y0)/20
            # y0 -= eps
            # y0 = max(y0, 1)
            # y1 += eps*5

            #y0 = 1
            #y1 = 1e+9
            ax[0].set_yscale('log')
            ax[0].set_ylim(y0, y1)
            # ax[0].set_ylim(1e-12, 1e-2)

            ax[1].set_ylim(ratio0,ratio1)
            ax[1].set_xlim(x0,x1)

            # plt.draw()


            # # check if upper limit in ratio on y-axis has a label (this will interfere with the upper pad)
            # cond1 = ax[1].get_ylim()[1] in list(map(lambda k: k.get_position()[1], ax[1].get_yticklabels()))

            # ax0_ticklabels_pos = list(map(lambda k: k.get_position()[1], ax[0].get_yticklabels()))
            # # check if the lower limit in the upper pad on y-axis has a label
            # cond2 = ax[0].get_ylim()[0] in ax0_ticklabels_pos

            # if cond1 and cond2:
            #     index = ax0_ticklabels_pos.index(ax[0].get_ylim()[0])
            #     plt.setp(ax[0].get_yticklabels()[index], visible=False)





            ax[0].set_ylabel(f'$d \sigma / d {variableNameLatex}$ [pb/{variableNameUnit}]', fontsize=30)
            # #ax[1].set_ylabel('Data / Pred.', loc='center', fontsize=20)
            ax[1].set_ylabel('BSM / SM', loc='center', fontsize=20)
            ax[1].set_xlabel(f'${variableNameLatex}$ [{variableNameUnit}]', fontsize=30)
            ax[1].tick_params(axis='y', labelsize='x-small', labelleft=True)


            ax[0].legend(ncol=2, fontsize='x-small', loc='upper center')

            def prettyPrintNumber(n):
                return str(round(n, 2)).replace('.', 'p')

            fig.savefig(f"{output_path}/plot_{variableName}_{operator1}_{prettyPrintNumber(operator1Value)}.png", facecolor='white', pad_inches=0.1, bbox_inches='tight')
            #fig.savefig(f"{output_path}test_plot_{variableName}_{operator1}_{prettyPrintNumber(operator1Value)}_{operator2}_{prettyPrintNumber(operator2Value)}.png", facecolor='white', pad_inches=0.1, bbox_inches='tight')