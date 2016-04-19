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
	variance = np.random.random_sample(dimensions)
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

def createDataFile(data):
	data_size = ','.join(map(str,data.shape))
	np.savetxt("data.csv", data,fmt='%10.8f', delimiter=",",header=data_size)

def createDataFile2(data):
	print(list(data.shape))
	np.insert(data, 0, list(data.shape))
	np.array(list(data.shape)).tofile("data.bin",sep="")
	np.array(list(data.shape)).tofile(open('data2.bin', 'wb'))
	'''
	newfile = open("data.dat",'wb')
	data_size = ','.join(map(str,data.shape))
	newfile.write(bytes(data_size+'\n', 'UTF-8'))
	for row in data:
		newfile.write(bytes(','.join(str(v) for v in value_list)+'\n', 'UTF-8'))
	newfile.close()
	'''

def createDataFile3(data):
	import struct
	sizeList = list(data.shape)
	out_file = open("data.bin","wb")
	s = struct.pack('i'*len(sizeList),*sizeList)
	out_file.write(s)
	for row in data:
		s = struct.pack('f'*len(row),*row)
		out_file.write(s)
	out_file.close()


##Create a function that creates many clusters and returns a vector with all the clusters
#
#data = np.random.uniform(0,1,size=(200,3))
data = generateClusters(2,3,10)
#print(data[0])
createDataFile3(data)
createDataFile(data)
print(data)
#plotClusters(data)