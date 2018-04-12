'''
import ThreeDplot
import Functions as mylibrary
import varna_draw as Varna
import heatmap as HMP
import plotstruct as PSt
import conf as conf 

import collections 
from collections import defaultdict
import Diana as Dn
import draw3dPareto
#rna=Functions.ParseFasta(rnafile,"RNAfold")
#numberssamples=len(constraintes)
'''

import conf , FileFunctions as FF, Sampling as SP, StructureFunctions as SF, Clustering as CL, VisualizationTools as VT, ClustersTrait as CT
import time,os,sys,pickle
from collections import defaultdict
from os.path import isfile, join
def GetListFile(PathFile, FileExtension):
    return [os.path.splitext(f)[0] for f in os.listdir(PathFile) if isfile(join(PathFile, f)) and os.path.splitext(f)[1] == '.' + FileExtension]
#Redirect all the print to Logfile.txt
sys.stdout =conf.Logger("Logfile.txt")

#rna = FF.Parsefile(conf.rnafile)[1]
#Indexe=conf.rnafile.split('.')[0]
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Loop over rna files!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

m=2.6/2
b=-0.8/2
path_Fasta = 'fasta_files'
Alignementfolder= 'Alignement'
FileExtensionFasta = 'fa'
print("Sampling Process for % s Structures" % (conf.numberofsruct))
OutputSamples = SP.StructSampling([conf.PathConstrainteFile, conf.PathConstrainteFileShape], Alignementfolder,conf.numberofsruct, conf.Temperature, conf.Fastaextenstion,m,b)

for filz in GetListFile(path_Fasta, FileExtensionFasta):
    print filz ,"Treatement "
    startimebig = time.time()
    rna = FF.Parsefile(os.path.join(path_Fasta, filz + '.' + FileExtensionFasta))[1]
    Indexe=filz
    SVMlFile="DissimilarityMatrix"+conf.numberofsruct
    listfiles = [filz + state for state in [ "NMIA","1M7","MSA"]]
    OutputSamples = 'OutputSamples' + conf.numberofsruct
    MFESnbrstruct=len(listfiles)# 1 for the case where no constraint is given
    FF.MergeFiles(OutputSamples,os.path.join(OutputSamples,'Samples.txt'),listfiles,1)
    #endtime=time.time()
    #print("Sampling done with success in  %53f\t"%(endtime-startime))

    #!!!!!!!!!!!!! Distance Matrix calculation !!!!!!!!!!!
    #
    startime = time.time()

    print("Distance Matrix generation for % d Structures started "%(int(conf.numberofsruct)*len(listfiles)+MFESnbrstruct))
    Redundant=[]
    Newnumberofsruct=defaultdict(lambda: defaultdict())
    #calculate distance
    SF.DistanceStruct(os.path.join(OutputSamples,'Samples.txt'),SVMlFile,int(conf.numberofsruct),MFESnbrstruct,listfiles)
    endtime=time.time()
    print("End of distance calculation between the structures in the sample  %53f\t"%(endtime-startime))


    #!!!!!!!!!!!!!Clustering!!!!!!!!!!!
    # load the pickled dissimilarity matrix
    DM=FF.UnpickleVariable(os.path.join(conf.PickledData,"dissmatrix.pkl"))
    Redundant=FF.UnpickleVariable(os.path.join(conf.PickledData,"Redondantestructures.pkl"))
    Redundant_Id=FF.UnpickleVariable(os.path.join(conf.PickledData,"Redondantestructures_Id.pkl"))
    Newnumberofsruct=FF.UnpickleVariable(os.path.join(conf.PickledData,"Dicnumberofsruct.pkl"))

    # clustering with DIANA Algorithm
    if conf.Diana=="true":
        startime = time.time()
        Clusters = defaultdict(list)
        structs = [i + 1 for i in range(len(DM))]
        clusters = CL.DIANA.doClustering(DM, structs,conf.maxDiameterThreshold, conf.maxAverageDiameterThreshold)
        for i in range(len(clusters)):
            Clusters[i] = clusters[i]
        endtime = time.time()
        print ("Clusters using Diana algorithm:  %s  %53f\t"%(Clusters,endtime-startime))

    #clustering with affinity propagation
    if conf.AFIN_PROP=="true":
        startime = time.time()
        startime=time.time()
        print("Start of clustering Process and Centroid (MEA) structure calculation")
        CL.AffinityPropagation(os.path.join("output",SVMlFile),Redundant)
        endtime = time.time()
        #Clusters=FF.UnpickleVariable(os.path.join(conf.PickledData,"Clusters_Aff_Prop.pkl"))
        #print ("Clusters using Affinity propagation: %s  %53f\t" % (Clusters, endtime - startime))
        print ("Clusters using Affinity propagation:   %53f\t" % ( endtime - startime))

    #!!!!!!!!!!!!!Boltzmann Probabilties calculation!!!!!!!!!!!

    #Get Clusters
    Clusters = FF.UnpickleVariable(os.path.join(conf.PickledData, "Clusters_Aff_Prop.pkl"))
    #print "here we are",Clusters


    #!!!!!!!!!!!!!Boltzmann Probabilties calculation!!!!!!!!!!!

    #Get Clusters
    #Clusters = FF.UnpickleVariable(os.path.join(conf.PickledData, "Clusters_Aff_Prop.pkl"))
    #The Normalization term is calculated relatively to the number of clusters
    Redondantestructure= FF.UnpickleVariable(os.path.join(conf.PickledData, "Redondantestructures_Id.pkl"))
    #print "redudant ", Redondantestructure
    print ("Start Boltzmann Probabilties calculation")
    startime = time.time()
    BltzProbinitial = defaultdict(lambda:defaultdict())
    BoltzmaanProbability = defaultdict(lambda: defaultdict())
    Zprobabilities = defaultdict(lambda: defaultdict())
    # Calling the Boltzmann_Calc to calculate all the probabilities needed for our study
    #print MFESnbrstruct,"by curiosity!!"
    #OutputSamples="OutputSamples12"
    #MFESnbrstruct=1
    print "Constraints!!", MFESnbrstruct
    BoltzmaanProbability=SF.Boltzmann_Calc(listfiles, OutputSamples,int(conf.numberofsruct), MFESnbrstruct,rna,Redondantestructure)
    endtime = time.time()
    print ("Boltzmann probabilities Calculated with success:   %53f\t" % ( endtime - startime))
    print "helaho",BoltzmaanProbability
    #!!!!!!!!!!!!!!!!!!!!!!!! Visualization !!!!!!!!!!!!!!!!!!!!!!!!!
    Clusters = FF.UnpickleVariable(os.path.join(conf.PickledData, "Clusters_Aff_Prop.pkl"))
    DM=FF.UnpickleVariable(os.path.join(conf.PickledData,"dissmatrix.pkl"))
    #BltzProbinitial= FF.UnpickleVariable(os.path.join(conf.PickledData, "Boltzman.pkl"))
    #BoltzmaanProbability=FF.UnpickleVariable(os.path.join(conf.PickledData,  "ConditionalBoltzman.pkl"))# This is the new probablities within the sample!
    #Todo there is a problem with the pickled Boltzann probbaility!!
    #print "helaho",BoltzmaanProbability
    #Zprobabilities=FF.UnpickleVariable(os.path.join(conf.PickledData,  "ZBolzman.pkl"))

    # visualize clusters with MFES, the function returns Clusters whithout the MFEs

    #VT.ThreeD_MatDistance_Boltzmann(DM, Clusters, BltzProbinitial, int(conf.numberofsruct), listfiles,conf.MFEs)

    Klusters =VT.FilterClusters(Clusters,int(conf.numberofsruct)*len(listfiles)+1)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!Clusters parameters calculation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Calculate centroid structure

    CentroidStructure = defaultdict(lambda:defaultdict())
    E = defaultdict()

    CentroidStructure, E, Matricdistancecentroids, ListDiameters =CT.CentroidBycluster(Klusters, os.path.join('OutputSamples' + str(conf.numberofsruct),'Samples.txt'),
                                                                                               BoltzmaanProbability,
                                                                                               int(conf.numberofsruct),
                                                                                       listfiles, float(conf.cutoff),rna)
    # We should create an intermediate  file to be sure that  RNAeval  works!!
    with open("filecentroide", "w")as Ctrdfile:
        Ctrdfile.write(">RNA")
        for elem in range(len(CentroidStructure)):
            Ctrdfile.write("\n%s" % (CentroidStructure[elem]))
            print CentroidStructure[elem]
    Centroids_Energies = SF.ENERGY_VALUES_STRUCTURES("filecentroide", rna)
    print "Centroids_Energies", Centroids_Energies


    CT.ClustersDistributions(Klusters, os.path.join("output","Clusters_details"), listfiles, int(conf.numberofsruct))
    #print("End of Clustering step, results in % s %53f\t" % ( endtime - startime))

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Election of the best structures strating!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Calculate cumulated  Boltzmaan Energy for each cluster

    CumulBE = {}
    CardinalConditions = {}
    EnergiesbyCluster = {}

    CumulBE = CT.CumulatedBoltzmanEnergiesbyCluster(Klusters, BoltzmaanProbability, int(conf.numberofsruct),listfiles)
    #epsilon= 1/n
    Epsilon = 1 / float(len(Klusters))
    CardinalConditions, EnergiesbyCluster = CT.GetCardinalConditions(Klusters, BoltzmaanProbability, listfiles,
                                                                            int(conf.numberofsruct), Epsilon)
    print 'verification pareto values'
    print 'Cardinal condition values', CardinalConditions.values()
    print 'Boltzmaan', CumulBE.values()
    print 'mean distance in a given cluster', E.values()
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Cluster analysis to elect the best cluster
    Dict = {} # a dictionary that contains the three variables characterizing clusters (Criteria of optimization)
    for ClusterNumber in Klusters:
        #Dict[ClusterNumber] = [CardinalConditions[ClusterNumber], CumulBE[ClusterNumber], -1. * E[ClusterNumber]]
        Dict[ClusterNumber] = [CardinalConditions[ClusterNumber], CumulBE[ClusterNumber]] #TODO only Two criterions are considered

    '''
    # Plot distribution of clusters in function of those criteria
    VT.plotClustercBECard(Klusters.keys(), CardinalConditions.values(), CumulBE.values(),
                                 'Number of conditions verifying Epsilon condition', '  Sum of Boltzmann energies',
                                 'Clusters_cardinal_cumulBltz.png')
    VT.plotClustercBECard(Klusters.keys(), CardinalConditions.values(), E.values(),
                                 'Number of conditions verifying Epsilon condition', ' Mean energy distance',
                                 'Clusters_Cardinal_meandistance.png')

    '''
    ListOptimalClusters = CT.Pareto(Dict)
    print "The elected clusters figuring in the Pareto front", ListOptimalClusters

    Dominatedclusters = [clusteri for clusteri in Klusters if clusteri not in ListOptimalClusters]
    #VT.plotPareto([Dict[elem] for elem in ListOptimalClusters], [Dict[elem] for elem in Dominatedclusters])

    #Fileallnonoptimalcentroides = os.path.join("output_plots",Indexe+".nonoptimals")
    #c, lista2 = VT.Drawvarna(Fileallnonoptimalcentroides, Dominatedclusters, CentroidStructure, conf.numberofsruct,
    #                                     rna, Centroids_Energies, conf.SHAPEVis)

    Filealloptimalcentroides = os.path.join("Multiprobing",Indexe+".optimals")
    c, lista = VT.Drawvarna(Filealloptimalcentroides, ListOptimalClusters, CentroidStructure, conf.numberofsruct,
                                         rna, Centroids_Energies, conf.SHAPEVis)

    #reactivities =FF.parseReactivityfile(conf.SHAPEVis)

   # VT.plotPairs(SF.OccurenceBasePairs(lista, 0), len(ListOptimalClusters))
    #VT.plotsecodnarystructures(rna, SF.OccurenceBasePairs(lista, 0), SF.OccurenceBasePairs(lista2, 0),len(ListOptimalClusters), reactivities)

    c.writeEPSfile("output_optimal_structures")

    print " RNAstruct has been run successfully"

    print ("End of RNAstruct run for a sampling of %s rna   % s structure for %s  conditions: %53f\t"%(filz,int(conf.numberofsruct) ,len(listfiles),time.time()-startimebig))
