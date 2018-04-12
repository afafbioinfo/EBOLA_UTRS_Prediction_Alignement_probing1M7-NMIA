from collections import defaultdict
import StructureFunctions as SF, FileFunctions as FF, VisualizationTools as VT
import numpy as np, scipy
import time


def CumulatedBoltzmanEnergiesbyCluster(clusters, ConditionalBoltzman, numberofsruct, constraintes):
	cBE = {}
	for ClusterNumber in clusters:
		l = 0.
		for structure in clusters[ClusterNumber]:
			ConditionNumber = int((structure - 1) / numberofsruct)
			StructureNumber = (structure - 1) - ConditionNumber * numberofsruct
			l += ConditionalBoltzman[constraintes[ConditionNumber]][StructureNumber]
		cBE[ClusterNumber] = l
	return cBE


# *************************Get cardinal for coditions verifying Epsilon test
def GetCardinalConditions(clusters, ConditionalBoltzman, constraintes, numberofsruct, Epsilon):
	EnergiesbyCluster = {}
	CardinalConditions = {}

	for ClusterNumber in clusters:
		# for each cluster, a sum overall boltzmann proba for a given condition is calculated
		# a condition counts for the cardinal of existing constrainte if the sum overall its strcutures is greater than epsilon[condition]
		l = defaultdict(lambda: 0)
		dic = {}
		for structure in clusters[ClusterNumber]:
			temp = 0
			ConditionNumber = int((structure - 1) / numberofsruct)
			StructureNumber = (structure - 1) - ConditionNumber * numberofsruct
			temp = ConditionalBoltzman[constraintes[ConditionNumber]][StructureNumber]
			l[ConditionNumber] += temp  # sum of Boltzm probabilites for a given condition
			dic[structure] = temp
		EnergiesbyCluster[ClusterNumber] = dic
		for ConditionNumber in range(len(constraintes)):
			if l[ConditionNumber] < Epsilon:
				del l[ConditionNumber]
		CardinalConditions[ClusterNumber] = len(l)
	return CardinalConditions, EnergiesbyCluster


# Cluster Diameters calculation
def ClustersDiameter(clusters, ListBPbystrcut):
	print "CLusters Diameters:"
	lista = []
	for ClusterNumber in clusters:
		if len(clusters[ClusterNumber]) > 1:
			d = max([SF.DistanceTwoBPlist(ListBPbystrcut[ClusterNumber][structure1], ListBPbystrcut[ClusterNumber][structure2])
				for structure1 in clusters[ClusterNumber] for structure2 in clusters[ClusterNumber]])
		else:
			d=0
		print "Cluster number", ClusterNumber, "diameter", d
		lista.append(d)
	if len(lista)!=0:
		print "Average distance", np.mean(lista)

	return lista


def BasePairsbyCluster(clusters, Structurefile, Boltzmann, numberofsruct, constrainte, cutoff):
	ListBPbyCluster = defaultdict()
	ListBPbyStruct = defaultdict(lambda: defaultdict())
	Zcluster = defaultdict(lambda: defaultdict())
	BoltzmannOverPairsbyCluster = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
	Fiftypercent = defaultdict()
	for ClusterNumber in clusters:
		liste = []
		for structure in clusters[ClusterNumber]:
			#print structure, 'clusters[ClusterNumber]',clusters[ClusterNumber]
			ConditionNumber = int((structure - 1) / numberofsruct)
			StructureNumber = (structure - 1) - ConditionNumber * numberofsruct
			#print 'ConditionNumber', ConditionNumber,'StructureNumber',StructureNumber, Boltzmann[constrainte[ConditionNumber]][StructureNumber]

			if(Boltzmann[constrainte[ConditionNumber]][StructureNumber]!=0):
				liste.append(Boltzmann[constrainte[ConditionNumber]][StructureNumber])
		print "LILSALISALISA",ClusterNumber,liste
		Zcluster[ClusterNumber] = sum(liste)

		list1 = []
		for structure in clusters[ClusterNumber]:

			ConditionNumber = int((structure - 1) / numberofsruct)
			StructureNumber = (structure - 1) - ConditionNumber * numberofsruct
			for (i, j) in SF.ListBasePairsFromStruct(FF.GetlinefromFile(Structurefile, structure - 1)):
				BoltzmannOverPairsbyCluster[ClusterNumber][i][j] += Boltzmann[constrainte[ConditionNumber]][
																		StructureNumber] / Zcluster[ClusterNumber]

			ListBPbyStruct[ClusterNumber][structure] = SF.ListBasePairsFromStruct(
				FF.GetlinefromFile(Structurefile, structure - 1))

		ListBPbyCluster[ClusterNumber] = list1
		Fiftypercent[ClusterNumber] = len(clusters) / 2
	return ListBPbyStruct, ListBPbyCluster, Fiftypercent, BoltzmannOverPairsbyCluster


def ClustersDistances(clusters, Boltzmanprobabilities, ListBPbystrcut, numberofsruct, condition):
	Emd = defaultdict()
	for ClusterNumber in clusters:
		liste = []
		for structure1 in clusters[ClusterNumber]:
			for structure2 in clusters[ClusterNumber]:
				# TODO Change the conversion to  be considered as fucntion!!!
				ConditionNumber1 = int((structure1 - 1) / numberofsruct)
				StructureNumber1 = (structure1 - 1) - ConditionNumber1 * numberofsruct
				ConditionNumber2 = int((structure2 - 1) / numberofsruct)
				StructureNumber2 = (structure2 - 1) - ConditionNumber2 * numberofsruct

				liste.append(Boltzmanprobabilities[condition[ConditionNumber1]][StructureNumber1] *
							 Boltzmanprobabilities[condition[ConditionNumber2]][
								 StructureNumber2] * SF.DistanceTwoBPlist(
					ListBPbystrcut[ClusterNumber][structure1], ListBPbystrcut[ClusterNumber][structure2]))
		#print "verify", ClusterNumber,2 * sum(liste) / float(len(clusters[ClusterNumber]) * (len(clusters[ClusterNumber]) - 1)),liste
		if len(clusters[ClusterNumber])!=0:
			Emd[ClusterNumber] = 2 * sum(liste) / float(len(clusters[ClusterNumber]) * (len(clusters[ClusterNumber]) - 1))
		else:
			Emd[ClusterNumber] =0
	#print "Mean distance in the clusters:", Emd
	return Emd


def BackTracing(W, BPp, rna, P_ss, i, j, pair):
    if i < j:
        # print i, j ,"what is the problem???",P_ss[j],W[i][j-1]
        if W[i][j] == P_ss[i] + W[i + 1][j]:
            # print "SOS1",i,j,W[i][j],P_ss[i],W[i+1][j]
            BackTracing(W, BPp, rna, P_ss, i + 1, j, pair)

        elif W[i][j] == P_ss[j] + W[i][j - 1]:
            # print "SOS2",i,j,W[i][j],P_ss[j],W[i][j-1]
            BackTracing(W, BPp, rna, P_ss, i, j - 1, pair)

        elif W[i][j] == (2 * BPp[i][j] + W[i + 1][j - 1]):
            pair.append((i, j))
            BackTracing(W, BPp, rna, P_ss, i + 1, j - 1, pair)
        else:
            for k in range(i, j):
                if W[i][j] == W[i][k] + W[k + 1][j]:
                    BackTracing(W, BPp, rna, P_ss, i, k, pair)
                    BackTracing(W, BPp, rna, P_ss, k + 1, j, pair)
                    break

    return pair

def MEA(BPp, rna, pair):
	startime = time.time()
	P_ss = defaultdict(lambda: 0)
	W = defaultdict(lambda: defaultdict(lambda: 0))
	for i in range(len(rna)):
		P_ss[i] = 1 - sum([BPp[min(i, j)][max(i, j)] for j in range(len(rna))])
	W = MEA_algo(P_ss, BPp, rna)
	endtime = time.time()
	print("End of MEA %53f\t" % (endtime - startime))
	# print BackTracing(W,BPp,rna,P_ss,0,len(rna)-1)
	# print fromPairsToStruct(rna, BackTracing(W,BPp,rna,P_ss,0,len(rna)-1,pair))
	return SF.fromPairsToStruct(rna, BackTracing(W, BPp, rna, P_ss, 0, len(rna) - 1, pair)), BackTracing(W, BPp, rna, P_ss,
																									  0, len(rna) - 1,
																									  pair)


def MEA_algo(P_ss, BPp, rna):
	startime = time.time()
	nb = len(rna)
	W = defaultdict(lambda: defaultdict(lambda: 0))
	for i in range(0, nb):
		W[i][i] = P_ss[i]
	for d in range(1, nb):
		for j in range(d, nb):
			i = j - d
			Res1 = P_ss[i] + W[i + 1][j]
			Res2 = P_ss[j] + W[i][j - 1]
			Res3 = 2 * BPp[i][j] + W[i + 1][j - 1]
			lista = []
			for k in range(i, j):
				lista.append(W[i][k] + W[k + 1][j])
			Res4 = np.max(lista)

			W[i][j] = np.max([Res1, Res2, Res3, Res4])

			# W[i][j]=np.max([P_ss[i]+W[i+1][j],P_ss[j]+W[i][j-1],2*BPp[i][j]+W[i+1][j-1],np.max([W[i][k]+W[k+1][j] for k in range(i,j)])])

	endtime = time.time()
	#print("End of W matrix fill in %53f\t" % (endtime - startime))
	return W


def CentroidBycluster(clusters, StructFile, Boltzmann, numberofsruct, constrainte, cutoff, rna):
	E = defaultdict()
	mycentroid = defaultdict()
	listpairscentroid = defaultdict(lambda: defaultdict())
	Myproba = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
	ListBPbystrcut, ListBP, Fiftypercent, Myproba = BasePairsbyCluster(clusters, StructFile, Boltzmann, numberofsruct,
																	   constrainte, cutoff)
	# print Diameters
	ListDiameters = ClustersDiameter(clusters, ListBPbystrcut)
	E = ClustersDistances(clusters, Boltzmann, ListBPbystrcut, numberofsruct, constrainte)
	for ClusterNumber in clusters:
		pair = []

		print "MEA algorithm is running for ClusterNumber", ClusterNumber
		mycentroid[ClusterNumber], listpairscentroid[ClusterNumber] = MEA(Myproba[ClusterNumber], rna, pair)

	print "Centroids calculation is done"
	# herein we add teh distance between clusters:
	print "Cluster distances"
	MatriceDistanceCentroids = scipy.zeros([len(clusters), len(clusters)])
	MatriceDistanceClustersEucld = scipy.zeros([len(clusters), len(clusters)])
	for ClusterNumber in clusters:
		for ClusterNumber2 in clusters:
			if ClusterNumber2 > ClusterNumber :
				l = SF.DistanceTwoBPlist(listpairscentroid[ClusterNumber], listpairscentroid[ClusterNumber2])
				# print "distance between clusters comparing the centroide's distances",l
				MatriceDistanceCentroids[ClusterNumber][ClusterNumber2] = l
				MatriceDistanceCentroids[ClusterNumber2][ClusterNumber] = l
				# print "distance between clusters comparing the means distances", ClusterNumber, ClusterNumber2, np.abs(E[ClusterNumber]-E[ClusterNumber2]),np.sqrt(abs(pow(E[ClusterNumber],2)-pow(E[ClusterNumber2],2)))
				print E
				l = np.sqrt(abs(pow(E[ClusterNumber], 2) - pow(E[ClusterNumber2], 2)))
				MatriceDistanceClustersEucld[ClusterNumber][ClusterNumber2] = l
				MatriceDistanceClustersEucld[ClusterNumber2][ClusterNumber] = l
				# print "distance between clusters compring the centroide's distances", ClusterNumber, ClusterNumber2, DistanceTwoBPlist(ListBPbystrcut[ClusterNumber][listCentroidStructure[ClusterNumber][0]],ListBPbystrcut[ClusterNumber2][listCentroidStructure[ClusterNumber2][0]])
	#VT.plotDistanceClusters(MatriceDistanceCentroids, clusters, "blue", " Base pair distance between centroids")
	#VT.plotDistanceClusters(MatriceDistanceClustersEucld, clusters, "red", "Eucledian distance between structures")
	return mycentroid, E, MatriceDistanceCentroids, ListDiameters


# count the occurence of present conditions in a given cluster
def SamplePresentInClusters(origine, occ, clusters, numberofsruct):
	for ClusterNumber in clusters:
		for StructureNumber in clusters[ClusterNumber]:
			origine[ClusterNumber][StructureNumber] = GetOrigineofStructure(StructureNumber, numberofsruct)
			# calculate occurence of conditions present whithin a cluster
	for ClusterNumber in clusters:
		for ConditionInCluster in origine[ClusterNumber].values():
			occ[ClusterNumber][ConditionInCluster] = origine[ClusterNumber].values().count(ConditionInCluster)

		return occ


# for a given structure characterized by a 'StructureNumber' return the condition represented by this structure
def GetOrigineofStructure(StructureNumber, numberofsruct):
	if (StructureNumber % numberofsruct != 0):
		return int(StructureNumber / numberofsruct) + 1
	else:
		return int(StructureNumber / numberofsruct)


def ClustersDistributions(clusters, Filename, constraintes, numberofsruct):
	origine = defaultdict(lambda: defaultdict(lambda: 0))
	occ = defaultdict(lambda: defaultdict(lambda: 0))
	it = defaultdict(lambda: 0)
	numberssamples = len(constraintes)
	o = open(Filename, "w")
	o.write("Cluster \t structures \t")
	for j in range(0, numberssamples):
		o.write("constraint %i = %s \t" % (j + 1, constraintes[j]))
	o.write("Number of structures \t Number of groups\n")
	occ = SamplePresentInClusters(origine, occ, clusters, numberofsruct)

	for elem in clusters:
		o.write("%i \t %s \t" % (elem + 1, clusters[elem]))
		for j in range(1, numberssamples + 1):
			if (occ[elem][j] != 0):
				it[elem] += 1
			o.write("%i\t" % (occ[elem][j]))
		o.write("%i\t%i\t" % (len(clusters[elem]), it[elem]))
		o.write("\n")

	o.write("Cluster(s) with  high number of  present conditions is(are) : %s" % (
		[v + 1 for v in it.keys() if it[v] == max(it.values())]))


# *************************!!!!!!!!!!!!!Pareto front*************************!!!!!!!!!!!!!!!!!!
def dominates(row, rowCandidate):
	return all(r >= rc for r, rc in zip(row, rowCandidate))


def Pareto(Dico):
    cleareCluster = []
    remaining = Dico.keys()

    while remaining:
        candidateKey = remaining.pop()
        candidate = Dico[candidateKey]

        new_remaining = []
        for other in remaining:
            if not dominates(candidate, Dico[other]):
                new_remaining.append(other)

        if not any(dominates(Dico[other], candidate) for other in new_remaining):
            cleareCluster.append(candidateKey)
        remaining = new_remaining
        #print len(remaining)
    # return cleareCluster,cleared
    return cleareCluster