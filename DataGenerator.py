import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plotClusters(data):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(data[:, 0], data[:, 1], data[:, 2], c='r', marker='o')
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
	plt.show()


def generateRandomCluster(dimensions, shift, sample_size):
	#variance = np.random.random_sample(dimensions)
	variance = shift * np.random.random_sample(3)
	covarianceMatrix = np.eye(dimensions)*variance
	mean = np.random.rand(dimensions)
	cluster = np.random.multivariate_normal(mean, covarianceMatrix, sample_size) + np.random.randint(low=0,high=shift, size=dimensions)
	return cluster

def generateClusters(cluster_count, dimensions, total_datapoints):
	data = []
	data_per_cluster = int(total_datapoints/cluster_count)
	for i in range(cluster_count):
		np.random.seed(i)
		random_shift = np.random.randint(20)
		cluster = generateRandomCluster(dimensions, random_shift, data_per_cluster)
		if data != []:
			data = np.concatenate((data, cluster), axis=0)
		else:
			data = cluster
	return data

def createDataFile_csv(data):
	#data_size = ','.join(map(str,data.shape))
	#np.savetxt("data.csv", data,fmt='%10.8f', delimiter=",",header=data_size)
	f = open(filename_csv, "w")
	f.write(','.join(map(str,data.shape))+"\n")
	for row in data:
		f.write(','.join(map(str,row))+"\n")


def createDataFile_binary(data):
	import struct
	sizeList = list(data.shape)
	out_file = open(filename_bin,"wb")
	s = struct.pack('i'*len(sizeList),*sizeList)
	out_file.write(s)
	for row in data:
		s = struct.pack('d'*len(row),*row)
		out_file.write(s)
	out_file.close()

def readBinaryFile():
	import numpy as npb
	memberships_parallel = np.fromfile(open("memberships_parallel.bin", "r"), dtype=np.uint32)
	memberships_serial = np.fromfile(open("memberships_serial.bin", "r"), dtype=np.uint32)
	memberships_serial = np.array(memberships_serial)
	print(len(memberships_serial))
	print(len(memberships_parallel))
	print(sum(memberships_parallel == memberships_serial)/len(memberships_serial))

def createFiles(numCluster, numDims, numDatapoints ,plot=False):
	np.random.seed(60)
	data = np.random.uniform(0,30,size=(numDatapoints,numDims))
	#data = generateClusters(numCluster, numDims, numDatapoints)
	createDataFile_binary(data)
	if plot:
		plotClusters(data)

##Create a function that creates many clusters and returns a vector with all the clusters
#
filename_bin = "uniform_data_3_10000.bin"
filename_csv = "uniform_data_3_10000.csv"
createFiles(5,3,10000000)
#readBinaryFile()
