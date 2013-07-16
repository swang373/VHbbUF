#!/usr/bin/env python
import os, sys, re
from math import sqrt
import decimal

wzfile = open("logWZ.log","r");
htmlfile = open("index.html","w");
htmlfile.write('<?xml version="1.0" encoding="iso-8859-1"?>')
htmlfile.write('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">')
htmlfile.write('<html>')
htmlfile.write('<head>')
htmlfile.write('<title>LumiMonitorFromEWKMu </title>')
htmlfile.write('<meta name="author" content="Michele de Gruttola">')
htmlfile.write('<link rel="stylesheet" href="mystyle.css" type="text/css" media="screen" charset=iso-8859-1>')
htmlfile.write('<style type="text/css"></style>')
htmlfile.write('<meta http-equiv="content-type" content="text/html; charset=UTF-8">')
htmlfile.write('</head>')
htmlfile.write('<body>')
htmlfile.write('<table class=main>')
htmlfile.write('<tr><td>')
htmlfile.write('<h2> Luminosity expectation counting W and Z into muons ( using dqm minimum bias(may27rereco) till runs 1357xx + run2010/Mu prompt-reco ), assuming B/S = 10% for W and 1% for Z; disclaimer:all numbers are very very preliminary </h2>')
htmlfile.write('</tr> </td>')
###lumi table
htmlfile.write('<TABLE border="1"')
htmlfile.write('<TR><TH>run<TH>LumiFromWmn(nW*0.24) 1 / nb<TH>LumiFromZmm(nZ*2.8) 1 / nb')

run = ''
nW = ''
nZ1HLT = ''
nZ2HLT= ''
nZNotIso= '' 
nZMuSta =''
nZMuTrk= ''
totLumiW=0
totLumiZ=0
    
for i in wzfile.readlines():
    line = i.split()
    if len(line)>2 :
           if (line[1][0:35] == 'NofW[50-200]fromEWKMuLumiMonitorDQM'):
             run= i.split()[0]
             nW= i.split()[2]
             print run, 'nW', nW
           if (line[1][0:8] == 'NofZ1HLT'):
            # run =  i.split()[0]
             nZ1HLT= i.split()[2]
             if run=='135445':
                 nZ1HLT=0  
             print run,  'nZ1hlt', nZ1HLT
           if (line[1][0:8] == 'NofZ2HLT'):
             run =  i.split()[0]
             nZ2HLT= i.split()[2]
             print run, 'nZ2HLT', nZ2HLT
           if (line[1][0:10] == 'NofZNotIso'):
             run =  i.split()[0]
             nZNotIso= i.split()[2]
             print run, 'nZNotIso', nZNotIso
           if (line[1][0:9] == 'NofZMuSta'):
             run =  i.split()[0]
             nZMuSta= i.split()[2]
             print run, 'nZMuSta', nZMuSta  
           if (line[1][0:9] == 'NofZMuTrk'):
             run =  i.split()[0]
             nZMuTrk= i.split()[2]
             print run, 'nZMuTrk', nZMuTrk  
           # filling the table before we go to next run...
             lumiW = int(nW) * 0.24
             lumiZ = (int(nZ1HLT) + int(nZ2HLT) + int(nZNotIso)) * 2.8
             totLumiW += lumiW
             totLumiZ += lumiZ
             print 'lumiZ', lumiZ
             htmlfile.write('<TR><TH>' +run + '<TH>' + str(lumiW) + '<TH>' + str(lumiZ)) 




sqrtTotLumiW = "%.1f" % sqrt(totLumiW * 0.24)
sqrtTotLumiZ = "%.1f" %  sqrt(totLumiZ * 2.8)



htmlfile.write('<TR><TH>tot<TH>' + str(totLumiW) + '+/-' + sqrtTotLumiW  + '<TH>' + str(totLumiZ) + '+/-' + sqrtTotLumiZ  ) 


htmlfile.write('</TABLE>')
#                   <TR><TH>Males<TD>1.9<TD>0.003<TD>40%
#                   <TR><TH>Females<TD>1.7<TD>0.002<TD>43%
               
#                   )         
#['137028', 'NofW[50-200]fromEWKMuDQMv1', '7']
#['137028', 'NofZ[60-120]fromEWKMuDQMv1', '2']
#['137028', 'NofW[50-200]fromEWKMuLumiMonitorDQMv2', '8']
#['137028', 'NofZ1HLT[60-120]fromEWKMuLumiMonitorDQMv1', '0']
#['137028', 'NofZ2HLT[60-120]fromEWKMuLumiMonitorDQMv1', '0']
#['137028', 'NofZNotIso[60-120]fromEWKMuLumiMonitorDQMv1', '1']




