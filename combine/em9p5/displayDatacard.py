from __future__ import absolute_import
from __future__ import print_function
import CombineHarvester.CombineTools.ch as ch
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('-i',"--inputFile", help="Input File",default=None)
parser.add_argument('-c',"--cats", help="Categories to print",default=None)
parser.add_argument('-o',"--odir", help="Output File directory",default=None)
parser.add_argument('--all', help="Print everything",action='store_true')
parser.add_argument('--obs', help="Print all observations",action='store_true')
parser.add_argument('--procs', help="Print all Processe",action='store_true')
parser.add_argument('--pars', help="Print all Parameters",action='store_true')
parser.add_argument('--syst', help="Print all syst",action='store_true')
args = parser.parse_args()



card_to_import='v30p1/v30SMCardStd.txt'
card_to_import=args.inputFile

print("= = "*40)
print("    CARD : ",card_to_import)
print("= = "*40)
cb = ch.CombineHarvester()
cb.ParseDatacard(card_to_import, analysis='Bs2EMu', channel="emu", mass='5.36688',era='13TeV')

if args.cats:
    categoriesToPrint=args.cats.split(',')
    print("Slimming to categories : ",categoriesToPrint)
    cb.bin(categoriesToPrint)

if args.all:
    cb.PrintAll()
if args.obs:
    cb.PrintObs()
if args.procs:
    cb.PrintProcs();
if args.pars:
    cb.PrintParams();
if args.syst:
    catList = cb.bin_set()
    for cat in catList:
        procList = cb.cp().bin([cat]).process_set()
        for prc in procList:
            cb_local=cb.cp().bin([cat]).process([prc])
            print("- - "*10)
            print(f"CAT : {cat}  | PROC : {prc} ")
            cb_local.PrintSysts()
        print("- - "*10,"\n")
        print("- - "*10,"\n")
    print("= =  "*10,"\n")
    print("- - "*10,"\n")
