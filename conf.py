from ConfigParser import SafeConfigParser 
import sys

#Connect to  the config file
config= SafeConfigParser()
config.read("RNAstruct.Config")
#Get parameters
rnafile=config.get("Paths","rnafile")
PathConstrainteFile=config.get("Paths","PathConstrainteFile")
PathConstrainteFileShape=config.get("Paths","PathConstrainteFileShape")
PickledData=config.get("Paths","PickledData")
numberofsruct=config.get("Conditions","numberofsruct")
MFEs=(config.get("Conditions","MFEs")).split(',')
Temperature=config.get("Conditions","Temperature")
SHAPEVis=config.get("Conditions","SHAPEVis")
percent=int(config.get("Pareto","percent"))
CutoffZcondition=float(config.get("Pareto","CutoffZcondition"))
cutoff=config.get("Conditions","cutoffBasePairs")
Psdotpath=config.get("Paths","Psdotpath")
Matrixproba=config.get("Paths","Matrixproba")
Fastaextenstion=config.get("Paths","Extension")
#constraintes=config.get("Conditions","Constraintes")
constraintes=(config.get("Conditions","constraints")).split(',')
AFIN_PROP=config.get("Clustering","AFIN_PROP")
MiniBatchKMean=config.get("Clustering","MiniBatchKMean")
Diana=config.get("Clustering","DIANA")
maxDiameterThreshold = float(config.get("Clustering","maxDiameterThreshold"))
maxAverageDiameterThreshold = float(config.get("Clustering","maxAverageDiameterThreshold"))
#function to cretae a log file
class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)