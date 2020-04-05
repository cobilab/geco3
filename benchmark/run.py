#!/usr/bin/env python3

import subprocess
import csv
import os

g2l16 = '-tm 1:1:0:0:0.9/0:0:0 -tm 2:1:0:0:0.78/0:0:0 -tm 3:1:0:0:0.9/0:0:0 -tm 4:1:1:0:0.78/0:0:0 -tm 5:10:1:1:0.90/0:0:0 -tm 6:1:1:0:0.85/0:0:0 -tm 7:1:1:0:0.85/0:0:0 -tm 8:1:1:0:0.91/0:0:0 -tm 9:10:0:0:0.9/0:0:0 -tm 10:10:0:0:0.9/0:0:0 -tm 11:10:0:0:0.9/0:0:0 -tm 12:20:1:1:0.94/0:0:0 -tm 13:10:1:0:0.95/0:0:0 -tm 14:50:1:1:0.95/0:0:0 -tm 16:200:1:10:0.95/1:50:0.95 -tm 17:100:1:20:0.9/3:10:0.9 -tm 20:500:1:30:0.95/2:20:0.95'

cfg1 = [[ '1', '0.06',  '8', 'DNACorpus/BuEb'],
        [ '2', '0.06', '16', 'DNACorpus/AgPh'],
        [ '3', '0.09', '24', 'DNACorpus/YeMi'],
        [ '4', '0.04', '40', 'DNACorpus/HePy'],
        [ '5', '0.04', '16', 'DNACorpus/AeCa'],
        [ '5', '0.04', '40', 'DNACorpus/HaHi'],
        [ '6', '0.03', '40', 'DNACorpus/EsCo'],
        [ '7', '0.03', '40', 'DNACorpus/PlFa'],
        [ '8', '0.03', '40', 'DNACorpus/ScPo'],
        [ '9', '0.05', '64', 'DNACorpus/EnIn'],
        ['10', '0.03', '64', 'DNACorpus/DrMe'],
        ['10', '0.03', '64', 'DNACorpus/OrSa'],
        ['10', '0.03', '64', 'DNACorpus/DaRe'],
        ['11', '0.03', '64', 'DNACorpus/GaGa'],
        ['12', '0.03', '64', 'DNACorpus/HoSa']]

cfg2 = [[ '1', '0.06',  '8', 'DNACorpus/BuEb'],
        [ '2', '0.06', '16', 'DNACorpus/AgPh'],
        [ '3', '0.09', '24', 'DNACorpus/YeMi'],
        ['16', '0.03', '64', 'DNACorpus/HePy'],
        ['16', '0.03', '64', 'DNACorpus/AeCa'],
        ['16', '0.03', '64', 'DNACorpus/HaHi'],
        ['16', '0.03', '64', 'DNACorpus/EsCo'],
        ['16', '0.03', '64', 'DNACorpus/PlFa'],
        ['16', '0.03', '64', 'DNACorpus/ScPo'],
        ['16', '0.03', '64', 'DNACorpus/EnIn'],
        ['16', '0.03', '64', 'DNACorpus/DrMe'],
        ['16', '0.03', '64', 'DNACorpus/OrSa'],
        ['16', '0.03', '64', 'DNACorpus/DaRe'],
        ['16', '0.03', '64', 'DNACorpus/GaGa'],
        ['16', '0.03', '64', 'DNACorpus/HoSa']]

cfg3 = [[ '8', '0.03',  '8',  'DNACorpus/ScPo'],
        [ '8', '0.03', '16',  'DNACorpus/ScPo'],
        [ '8', '0.03', '24',  'DNACorpus/ScPo'],
        [ '8', '0.03', '32',  'DNACorpus/ScPo'],
        [ '8', '0.03', '40',  'DNACorpus/ScPo'],
        [ '8', '0.03', '48',  'DNACorpus/ScPo'],
        [ '8', '0.03', '56',  'DNACorpus/ScPo'],
        [ '8', '0.03', '64',  'DNACorpus/ScPo'],
        [ '8', '0.03', '72',  'DNACorpus/ScPo'],
        [ '8', '0.03', '80',  'DNACorpus/ScPo'],
        [ '8', '0.03', '88',  'DNACorpus/ScPo'],
        [ '8', '0.03', '96',  'DNACorpus/ScPo'],
        [ '8', '0.03', '104', 'DNACorpus/ScPo'],

        [ '9', '0.03',  '8',  'DNACorpus/EnIn'],
        [ '9', '0.03', '16',  'DNACorpus/EnIn'],
        [ '9', '0.03', '24',  'DNACorpus/EnIn'],
        [ '9', '0.03', '32',  'DNACorpus/EnIn'],
        [ '9', '0.03', '40',  'DNACorpus/EnIn'],
        [ '9', '0.03', '48',  'DNACorpus/EnIn'],
        [ '9', '0.03', '56',  'DNACorpus/EnIn'],
        [ '9', '0.03', '64',  'DNACorpus/EnIn'],
        [ '9', '0.03', '72',  'DNACorpus/EnIn'],
        [ '9', '0.03', '80',  'DNACorpus/EnIn'],
        [ '9', '0.03', '88',  'DNACorpus/EnIn'],
        [ '9', '0.03', '96',  'DNACorpus/EnIn'],
        [ '9', '0.03', '104', 'DNACorpus/EnIn'],

        ['10', '0.03',  '8',  'DNACorpus/DrMe'],
        ['10', '0.03', '16',  'DNACorpus/DrMe'],
        ['10', '0.03', '24',  'DNACorpus/DrMe'],
        ['10', '0.03', '32',  'DNACorpus/DrMe'],
        ['10', '0.03', '40',  'DNACorpus/DrMe'],
        ['10', '0.03', '48',  'DNACorpus/DrMe'],
        ['10', '0.03', '56',  'DNACorpus/DrMe'],
        ['10', '0.03', '64',  'DNACorpus/DrMe'],
        ['10', '0.03', '72',  'DNACorpus/DrMe'],
        ['10', '0.03', '80',  'DNACorpus/DrMe'],
        ['10', '0.03', '88',  'DNACorpus/DrMe'],
        ['10', '0.03', '96',  'DNACorpus/DrMe'],
        ['10', '0.03', '104', 'DNACorpus/DrMe']]

def results(fname, cfg):
    with open(fname, 'w', newline='') as csvfile:
        lines = []
        fieldnames = ['ID', 'GeCo2 bytes', 'GeCo3 bytes', 'GeCo2 secs', 'GeCo3 secs', 'Mode', 'L.Rate', 'H.Nodes']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        totalb2 = 0
        totalb3 = 0

        totals2 = 0
        totals3 = 0

        for [l, lr, hs, f] in cfg:
            if l == '16':
                cmd = ['./GeCo2', '-lr', lr, '-hs', hs, f]
                cmd[1:1] = g2l16.split()
                out = subprocess.check_output(cmd)
            else:
                out = subprocess.check_output(['./GeCo2','-l', l, '-lr', lr, '-hs', hs, f])

            sout = str(out)

            bytes2 = sout.split('Total bytes: ')[1].split()[0]
            bytes2 = int(bytes2)
            totalb2 = totalb2 + bytes2
            bytes2 = format(bytes2, ',d')

            secs2 = sout.split('Spent ')[1].split()[0]
            secs2 = float(secs2)
            totals2 = totals2 + secs2
            secs2 = str(round(secs2, 1))


            out = subprocess.check_output(['./GeCo3','-l', l, '-lr', lr, '-hs', hs, f])
            sout = str(out)
            bytes3 = sout.split('Total bytes: ')[1].split()[0]
            bytes3 = int(bytes3)
            totalb3 = totalb3 + bytes3
            bytes3 = format(bytes3, ',d')

            secs3 = sout.split('Spent ')[1].split()[0]
            secs3 = float(secs3)
            totals3 = totals3 + secs3
            secs3 = str(round(secs3, 1))

            d = {'ID': os.path.basename(f),
                 'GeCo2 bytes': bytes2, 'GeCo3 bytes': bytes3,
                 'GeCo2 secs': secs2, 'GeCo3 secs': secs3,
                 'Mode': l, 'L.Rate': lr, 'H.Nodes': hs}
            lines.append(d)
            print(d)

        for l in reversed(lines):
            writer.writerow(l)

        totalb2 = format(totalb2, ',d')
        totalb3 = format(totalb3, ',d')

        totals2 = str(round(totals2, 1))
        totals3 = str(round(totals3, 1))

        d = {'ID': 'Total',
             'GeCo2 bytes': totalb2, 'GeCo3 bytes': totalb3,
             'GeCo2 secs': totals2, 'GeCo3 secs': totals3,
             'Mode': '', 'L.Rate': '', 'H.Nodes': ''}
        writer.writerow(d)

def resultsh(fname, cfg):
    with open(fname, 'w', newline='') as csvfile:
        lines = []
        fieldnames = ['ID', 'GeCo3 bytes', 'GeCo3 secs', 'L.Rate', 'H.Nodes']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        totalb2 = 0
        totalb3 = 0

        totals2 = 0
        totals3 = 0

        for [l, lr, hs, f] in cfg:
            out = subprocess.check_output(['./GeCo3','-l', l, '-lr', lr, '-hs', hs, f])
            sout = str(out)
            bytes3 = sout.split('Total bytes: ')[1].split()[0]
            bytes3 = int(bytes3)
            totalb3 = totalb3 + bytes3
            bytes3 = format(bytes3, ',d')

            secs3 = sout.split('Spent ')[1].split()[0]
            secs3 = float(secs3)
            totals3 = totals3 + secs3
            secs3 = str(round(secs3, 1))

            d = {'ID': os.path.basename(f),
                 'GeCo3 bytes': bytes3,
                 'GeCo3 secs': secs3,
                 'L.Rate': lr,
                 'H.Nodes': hs}
            lines.append(d)
            print(d)

        for l in reversed(lines):
            writer.writerow(l)


#results('m16.csv', cfg2)
results('original.csv', cfg1)

#resultsh('hidden.csv', cfg3)
