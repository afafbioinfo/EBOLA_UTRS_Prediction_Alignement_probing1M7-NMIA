import FileFunctions as FF
import os, subprocess
import conf as conf
def Adjust_structure(File,directory, IndexRna):
	rnaaligned= FF.Parsefile(File)[0]
        listdeletion=[i for i, j in enumerate(rnaaligned) if j == '_']
	lines = FF.Parsefile(File)
    	#SeqLen=len(lines[1])-1
	#print SeqLen,"seq length"
        fileout=os.path.join(directory,IndexRna+'MSA')
    	o = open(fileout, "w")
        o.write('>'+IndexRna+'MSA\n')
        for struct in lines[2:-2]:
		struc=list(struct)
		for k in listdeletion:
			#print struct
			struc[k]=""
    		o.write('%s'%("".join(struc)))
	o.close()

def StructSampling(Pathconstraints,InputAlignment, numberofsruct, Tmpr,Extension,m,b):
    print InputAlignment
    dir='OutputSamples' + str(numberofsruct)
    FF.CreateFold(dir)
    
    for Path in Pathconstraints:
        for filename in FF.GetListFile(Path, Extension):
            print "filename",filename
            Input = Path + "/" + filename
            output = dir + '/' + os.path.splitext(filename)[0]
            #Command = 'RNAsubopt  -p ' + str(numberofsruct) + ' -s -T ' + str(Tmpr)+ ' --maxBPspan '+ str(3*444/2)
            Command = 'RNAsubopt  -p ' + str(numberofsruct) + ' -s -T ' + str(Tmpr)
            if Path ==conf.PathConstrainteFile:
                Command += ' -C  --enforceConstraint'
            if Path == conf.PathConstrainteFileShape:
                ShapeFile = conf.PathConstrainteFileShape + "/" + os.path.splitext(filename)[0] + 'Shape.txt'
                #Command += ' --shape ' + ShapeFile + ' shapeMethod=Z '+ 'shapeConversion= M' # by default the shapeMethod=D (Deigan et al 2009), where the method Z is from Zarringhalm et al 2012
                #Command += ' --shape ' + ShapeFile + ' --shapeMethod="Dm3.48b-1.35"' # DMS new param
                #Command += ' --shape ' + ShapeFile #" defaults
                #Command += ' --shape ' + ShapeFile + ' --shapeMethod="Dm3.4794b-1.1598"' # for CMCT
                #print Command
                Command += ' --shape ' + ShapeFile + ' --shapeMethod="Dm'+str(m)+'b'+str(b)+'"'
            subprocess.call(Command, stdin=open(Input, 'r'), stdout=open(output, 'wb'), shell=True)
            #because the version 2.3 of rnaeval does not consider the rna, second line should be removed


            lines = []

            with open(output, 'r') as f:
                lines = f.readlines()

            with open(output, 'w') as f:
                f.writelines(lines[:1] + lines[2:])
    
    for filename in FF.GetListFile(InputAlignment, '.aln'):
	    Input= InputAlignment + "/" + filename
            print InputAlignment , Input
	    Output='Sample'+ os.path.splitext(filename)[0]
	    #For the sampling from MSA alifold
	    Command = 'RNAalifold -r -d2 --noLP  --aln --stochBT '+ str(numberofsruct)
	    subprocess.call(Command, stdin=open(Input, 'r'), stdout=open(Output, 'wb'), shell=True)
	    Adjust_structure(Output,dir,os.path.splitext(filename)[0])
    return dir


