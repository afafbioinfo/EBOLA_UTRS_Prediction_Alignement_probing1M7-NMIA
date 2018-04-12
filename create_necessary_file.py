#rename "s/.shape/1M7Shape.txt/" *.shape
#rename "s/.shape/NMIAShape.txt/" *.shape
import os, subprocess, time
from os.path import isfile, join
def openfile(Path):
    fileIn = open(Path, "r")
    lines = fileIn.readlines()
    fileIn.close()
    return lines

def Parsefile(Path):
    fileIn = open(Path, "r")
    lines = fileIn.readlines()
    fileIn.close()
    return lines

def GetListFile(PathFile, FileExtension):
    return [os.path.splitext(f)[0] for f in os.listdir(PathFile) if isfile(join(PathFile, f)) and os.path.splitext(f)[1] == '.' + FileExtension]
def CreateFilefasta(rna,path_Shape,reagent):    
    Outfilerna=open(os.path.join(path_Shape, filz + reagent+'.fa'),'w')
    Outfilerna.write("%s \n"%(">"+filz + reagent))
    Outfilerna.write("%s"%(rna))
    Outfilerna.close()

path="fasta_files"
FileExtension="fa"
reagents=["NMIA","1M7"]
path_Shape="fasta_Shape"

for filz in GetListFile(path, FileExtension):
        print filz
        Fastafile= os.path.join(path, filz + '.' + FileExtension)
        rna = openfile(Fastafile)[1]
        for state in reagents:
        	CreateFilefasta(rna,path_Shape,state) 
