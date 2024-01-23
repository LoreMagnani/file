#!/usr/bin/env python

import os
import sys
#import ConfigParser
import subprocess

if __name__ == "__main__":

    #switchOn = ['2','4','5','7','9','21','22','24','25','29','30','31','32','33','34']
    #switchOn = ['2', '4', '7', '21', '24', '25', '31', '32', '33', '34'] # production succeeded
    #switchOn = [str(i) for i in range(16, 61)]
    switchOn = [str(i) for i in range(1, 54)]
    #switchOn = ['2']
    params = [( '1' , 'cGtil'),
            ('2' , 'cWtil'),
            ('3' , 'cHGtil'),
            ('4' , 'cHWtil'),
            ('5' , 'cHBtil'),
            ('6' , 'cHWBtil'),
            ('7' , 'cuGIm'),
            ('8' , 'ctGIm'),
            ('9' , 'cuWIm'),
            ('10' , 'ctWIm'),
            ('11' , 'cuBIm'),
            ('12' , 'ctBIm'),
            ('13' , 'cdGIm'),
            ('14' , 'cbGIm'),
            ('15' , 'cdWIm'),
            ('16' , 'cbWIm'),
            ('17' , 'cdBIm'),
            ('18' , 'cbBIm'),
            ('19' , 'cuHIm'),
            ('20' , 'ctHIm'),
            ('21' , 'cdHIm'),
            ('22' , 'cbHIm'),
            ('23' , 'cHudIm'),
            ('24' , 'cHtbIm'),
            ('25' , 'cutbd1Im'),
            ('26' , 'cutbd8Im'),
            ('27' , 'cjQtu1Im'),
            ('28' , 'cjQtu8Im'),
            ('29' , 'cjQbd1Im'),
            ('30' , 'cjQbd8Im'),
            ('31' , 'cjujd1Im'),
            ('32' , 'cjujd8Im'),
            ('33' , 'cjujd11Im'),
            ('34' , 'cjujd81Im'),
            ('35' , 'cQtjd1Im'),
            ('36' , 'cQtjd8Im'),
            ('37' , 'cjuQb1Im'),
            ('38' , 'cjuQb8Im'),
            ('39' , 'cQujb1Im'),
            ('40' , 'cQujb8Im'),
            ('41' , 'cjtQd1Im'),
            ('42' , 'cjtQd8Im'),
            ('43' , 'cQtQb1Im'),
            ('44' , 'cQtQb8Im'),
            ('45' , 'ceHIm'),
            ('46' , 'ceWIm'),
            ('47' , 'ceBIm'),
            ('48' , 'cledjIm'),
            ('49' , 'clebQIm'),
            ('50' , 'cleju1Im'),
            ('51' , 'cleju3Im'),
            ('52' , 'cleQt1Im'),
            ('53' , 'cleQt3Im')]
    non_null_operators = []
    for param in params:
        if param[0] not in switchOn : continue
        if 'launch_Zjj_' + param[1] + '_QU.txt' in os.listdir('.'): continue

        f_launchfile = open ('launch_Zjj_' + param[1] + '_QU.txt', 'w')
        f_launchfile.write ('import model SMEFTsim_topU3l_MwScheme_UFO-' + param[1] + '_massless\n')
        f_launchfile.write('define p = g u c d s b u~ c~ d~ s~ b~\n')
        f_launchfile.write('define j = p\n')
        f_launchfile.write('define l+ = e+ mu+ \n')
        f_launchfile.write('define l- = e- mu- \n')
        f_launchfile.write ('generate p p > l+ l- j j SMHLOOP=0 QED=99 QCD=0 NP=1 NP^2==2\n')
        f_launchfile.write ('output Zjj_' + param[1] + '_QU')
        f_launchfile.close ()
        p = 'Zjj_'+param[1]+'_QU'
        
        cmd = subprocess.call("../MG5_aMC_v2_6_5/bin/mg5_aMC launch_"+p+".txt", shell=True)
        if not p in os.listdir('.'):
            print("folder for LI was not created by madgraph, exit")
            continue
        else:
            non_null_operators.append(param[0]+"QU")

        if 'launch_Zjj_' + param[1] + '_LI.txt' in os.listdir('.'): continue

        f_launchfile = open ('launch_Zjj_' + param[1] + '_LI.txt', 'w')
        f_launchfile.write ('import model SMEFTsim_topU3l_MwScheme_UFO-' + param[1] + '_massless\n')
        f_launchfile.write('define p = g u c d s b u~ c~ d~ s~ b~\n')
        f_launchfile.write('define j = p\n')
        f_launchfile.write('define l+ = e+ mu+ \n')
        f_launchfile.write('define l- = e- mu- \n')
        f_launchfile.write ('generate p p > l+ l- j j SMHLOOP=0 QED=99 QCD=0 NP=1 NP^2==1\n')
        f_launchfile.write ('output Zjj_' + param[1] + '_LI')
        f_launchfile.close ()
        p = 'Zjj_'+param[1]+'_LI'
        cmd = subprocess.call("../MG5_aMC_v2_6_5/bin/mg5_aMC launch_"+p+".txt", shell=True)
        if not p in os.listdir('.'):
            print("folder for LI was not created by madgraph, exit")
        else:
            non_null_operators.append(param[0]+"LI")

    print(non_null_operators)
