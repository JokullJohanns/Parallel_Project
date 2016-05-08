#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <sys/time.h>

double** datapoints;  // The datapoints to cluster
int* membership;      // Which cluster each datapoint belongs to
double** centroids;    // Means of each cluster
double** newCentroids; // Sum of the dimensions of datapoints belonging to each cluster
int* newCentroidsSize;  // Number of datapoints belonging to each cluster

int numDatapoints = 0;
int numDims = 0;
int K = 2;
float threshold;

void printCentroids() {
	int i,j;
    for(i = 0; i < K; i++) {
        for(j = 0; j < numDims; j++) {
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

void readFile(char* filename){
    FILE* stream = fopen(filename, "r");

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
		int column;
        for(column = 0; column < numDims; column++){
            char* tmp = strdup(line);
            double value;
            sscanf(getfield(tmp, column+1), "%lf", &value);
            datapoints[row][column] = value;
            free(tmp);
        }
        row++;
    }
}

void readBinary(char* filename){
    FILE *ptr;

    ptr = fopen(filename,"rb");  // r for read, b for binary

    fread(&numDatapoints, sizeof(int),1,ptr); // read 10 bytes to our buffer
    fread(&numDims, sizeof(int),1,ptr);
    datapoints = allocMatrix(numDatapoints, numDims, 0);
	printf("Datapoints(%d), dims(%d)\n", numDatapoints, numDims);
    fread(datapoints[0], sizeof(float)*numDatapoints*numDims,1,ptr);
}

void initVars() {
    //membership  = (unsigned char *) malloc(numDatapoints * sizeof(unsigned char));
    membership = (int *) calloc(numDatapoints, sizeof(int));
    centroids    = allocMatrix(K, numDims, 0);
    newCentroids = allocMatrix(K, numDims, 1);
    newCentroidsSize = (int *) calloc(numDatapoints, sizeof(int));
	int k, i;
    for(k = 0; k < K; k++) {
        for(i = 0; i < numDims; i++) {
            centroids[k][i] = datapoints[k][i];
        }
    }
}

double calcDistance(int pointIndex, int centroidIndex) {
    double sum = 0.0;
	int i;
    for(i = 0; i < numDims; i++) {
        sum = sum + pow(datapoints[pointIndex][i] - centroids[centroidIndex][i], 2);
    }
    return sqrt(sum);
}

void sumDatapointAndNewCentroid(pointIndex, centroidIndex) {
	int i;    
	for(i = 0; i < numDims; i++) {
        newCentroids[centroidIndex][i] += datapoints[pointIndex][i];
    }
}

void calcNewCentroids() {
	int k,d;
    for(k = 0; k < K; k++) {
        for(d = 0; d < numDims; d++) {
            centroids[k][d] = newCentroids[k][d] / newCentroidsSize[k];
            //printf("newCentroid %d %d = %lf / %d = %lf\n", k, d, newCentroids[k][d], newCentroidsSize[k], centroids[k][d]);
            newCentroids[k][d] = 0.0;
        }
        newCentroidsSize[k] = 0;
    }
    /*
    for(int i = 0; i < K; i++){
        printf("Rank %d has global_centroids(%d) [%.8lf, %.8lf, %.8lf]\n",0,i, centroids[i][0], centroids[i][1], centroids[i][2]);
    }
    */
}

void kMeans() {
    float numChanged = 0; // Tracks how many datapoints changed clusters each iteration
    double minDistance;
    double distance = 0;
    int curCluster;
    int itercount = 0;
    do {
        numChanged = 0.0;
		int pointIndex;
        for(pointIndex = 0; pointIndex < numDatapoints; pointIndex++) { // iterate through datapoints
            minDistance = LONG_MAX;
			int k;
            for(k = 0; k < K; k++) { // iterate throught cluster centroids
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
        //printf("Datapoints that changed clusters: %f\n", numChanged/numDatapoints);
    } while((numChanged/numDatapoints > threshold) && itercount++ < 500);
}



void printDatapoints() {
	int i,j;
    for(i = 0; i < numDatapoints; i++) {
        for(j = 0; j < numDims; j++) {
            printf("%.8f, ", datapoints[i][j]);
        }
        printf("\n");
    }

    printf("******************\n");

    printCentroids();
}

void printDatapointClusters() {
    FILE *write_ptr;
    remove("memberships_serial.bin");
    write_ptr = fopen("memberships_serial.bin","wb");  // w for write, b for binary

    fwrite(membership,sizeof(int)*numDatapoints,1,write_ptr); // write 10 bytes to our buffer
    /*
    for(int i = 0; i < numDatapoints; i++) {
        printf("%d,", membership[i]);
    }
    */
}

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


int main(int argc, char **argv)
{
	if(argc < 2){
		printf("No input file given");
		exit(0);
	}
	K = 5;
    readBinary(argv[1]);
    //readFile("uniform_data_16_1000.csv");
    
    initVars();
    threshold = 0.001;
    //printf("Centroids at start:\n");
    //printCentroids();
    //printf("Starting k-Means\n");
    //printDatapoints();
    //time_t start = time(NULL);
    double startTime = get_wall_time();
    kMeans();
    double endTime = get_wall_time();
    printf("%.2f\n", endTime-startTime);
    
    printDatapointClusters();
    
    return 0;
    
}
