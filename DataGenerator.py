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

##Create a function that creates many clusters and returns a vector with all the clusters
#
#data = np.random.uniform(0,1,size=(200,3))
data = generateClusters(2,3,10)
#print(data[0])
createDataFile(data)
#plotClusters(data)