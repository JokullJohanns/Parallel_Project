#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>

double** datapoints;  // The datapoints to cluster
unsigned char* membership;      // Which cluster each datapoint belongs to
double** centroids;    // Means of each cluster
double** newCentroids; // Sum of the dimensions of datapoints belonging to each cluster
int* newCentroidsSize;  // Number of datapoints belonging to each cluster

int numDatapoints = 0;
int numDims = 0;
int K = 2;


void printCentroids() {
    for(int i = 0; i < K; i++) {
        for(int j = 0; j < numDims; j++) {
            printf("%.8f, ", centroids[i][j]);
        }
        printf("\n");
    }
}


double** allocMatrix(int rows, int cols, int useCalloc)
{
    int i;
    int header= rows *sizeof(double *);
    int data=rows*cols*sizeof(double);
    double ** rowptr;
    if(useCalloc == 1) {
        rowptr=(double **)malloc(header+data);
    }
    else {
        rowptr=(double **)malloc(header+data);    
    }
    if(rowptr==NULL){
        return NULL;
    }
    double * buf=(double*)(rowptr+rows);
    for(i=0;i<rows;i++){
        rowptr[i]=buf+i*cols;
    } 
    return rowptr;
 }


char* getfield(char* line, int num){
    char* tok;
    for (tok = strtok(line, ","); tok && *tok; tok = strtok(NULL, ",\n"))
    {
        if (!--num)
            return tok;
    }
    return NULL;
}

void readFile(){
    FILE* stream = fopen("data.csv", "r");

    // First extract the number of datapoints and their dimension
    char firstLine[1024];
    fgets(firstLine, 1024, stream);
    numDatapoints = atoi(getfield(strdup(firstLine), 1)); // num datapoints
    numDims = atoi(getfield(strdup(firstLine), 2)); // 

    // Allocate array to store matrix
    datapoints = allocMatrix(numDatapoints, numDims, 0);

    // Read the datapoints, one line at a time
    char line[1024];
    int row = 0;
    while (fgets(line, 1024, stream))
    {
        for(int column = 0; column < numDims; column++){
            char* tmp = strdup(line);
            double value;
            sscanf(getfield(tmp, column+1), "%lf", &value);
            datapoints[row][column] = value;
            free(tmp);
        }
        row++;
    }
}

void initVars() {
    membership  = (unsigned char *) malloc(numDatapoints * sizeof(unsigned char));
    centroids    = allocMatrix(K, numDims, 0);
    newCentroids = allocMatrix(K, numDims, 1);
    newCentroidsSize = (int *) calloc(numDatapoints, sizeof(int));

    for(int k = 0; k < K; k++) {
        for(int i = 0; i < numDims; i++) {
            centroids[k][i] = datapoints[k][i];
        }
    }
}

double calcDistance(int pointIndex, int centroidIndex) {
    double sum = 0.0;
    for(int i = 0; i < numDims; i++) {
        sum = sum + pow(datapoints[pointIndex][i] - centroids[centroidIndex][i], 2);
    }
    return sqrt(sum);
}

void sumDatapointAndNewCentroid(pointIndex, centroidIndex) {
    for(int i = 0; i < numDims; i++) {
        newCentroids[centroidIndex][i] += datapoints[pointIndex][i];
    }
}

void calcNewCentroids() {
    for(int k = 0; k < K; k++) {
        for(int d = 0; d < numDims; d++) {
            centroids[k][d] = newCentroids[k][d] / newCentroidsSize[k];
            //printf("newCentroid %d %d = %lf / %d = %lf\n", k, d, newCentroids[k][d], newCentroidsSize[k], centroids[k][d]);
            newCentroids[k][d] = 0.0;
        }
        newCentroidsSize[k] = 0;
    }
}

void kMeans() {
    int numChanged = 0; // Tracks how many datapoints changed clusters each iteration
    double minDistance;
    double distance = 0;
    int curCluster;

    for(int i = 0; i < 10; i++) { // nubmer of iterations
        numChanged = 0; // Tracks how many datapoints changed clusters each iteration
        for(int pointIndex = 0; pointIndex < numDatapoints; pointIndex++) { // iterate through datapoints
            minDistance = LONG_MAX;
            for(int k = 0; k < K; k++) { // iterate throught cluster centroids
                distance = calcDistance(pointIndex, k);
                if(distance < minDistance) {
                    minDistance = distance;
                    curCluster = k;
                }
            }
            if(membership[pointIndex] != curCluster) {
                numChanged += 1;
                membership[pointIndex] = curCluster;
            }
            sumDatapointAndNewCentroid(pointIndex, curCluster);
            newCentroidsSize[curCluster] = newCentroidsSize[curCluster] + 1;
        }

        // calculate new cluster centers
        calcNewCentroids();
        //printCentroids();
        printf("Datapoints that changed clusters: %d\n", numChanged);
    }
}



void printDatapoints() {
    for(int i = 0; i < numDatapoints; i++) {
        for(int j = 0; j < numDims; j++) {
            printf("%.8f, ", datapoints[i][j]);
        }
        printf("\n");
    }

    printf("******************\n");

    printCentroids();
}

void printDatapointClusters() {
    for(int i = 0; i < numDatapoints; i++) {
        printf("%d\n", membership[i]);
    }
}


int main(int argc, char **argv)
{
    readFile();
    initVars();
    printf("Centroids at start:\n");
    printCentroids();
    printf("Starting k-Means\n");
    //printDatapoints();
    kMeans();
    printDatapointClusters();
    return 0;
}