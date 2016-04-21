#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

MPI_Status status;
double** local_datapoints;  // The datapoints to cluster
int* local_membership;      // Which cluster each datapoint belongs to
double** local_newCentroids; // Sum of the dimensions of datapoints belonging to each cluster
int* local_newCentroidsSize;  // Number of datapoints belonging to each cluster
double** global_centroids;    // Means of each cluster
int* global_centroidsSize;    

int global_datapoint_count, processCount,Ip,rank;
int local_datapoint_count;
int numDims;
int K, tag;
float threshold;

double** allocMatrix(int rows, int cols, int useCalloc)
{
    int i;
    int header= rows *sizeof(double *);
    int data=rows*cols*sizeof(double);
    double ** rowptr;
    if(useCalloc == 1) {
        rowptr = (double **) calloc(header+data, sizeof(double));
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


void mpi_read(char* filename)
{
    
    int  i, len, divd, rem;
    MPI_Status status;
    int            err;
    MPI_Offset     disp;
    MPI_Datatype   filetype;
    MPI_File       fh;

    err = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (err != MPI_SUCCESS) {
        char errstr[MPI_MAX_ERROR_STRING];
        int  errlen;
        MPI_Error_string(err, errstr, &errlen);
        printf("Error at opening file %s (%s)\n",filename,errstr);
        MPI_Finalize();
        exit(1);
    }

    MPI_File_read(fh, &global_datapoint_count,   1, MPI_INT, &status);
    MPI_File_read(fh, &numDims, 1, MPI_INT, &status);

    if (global_datapoint_count <= 0 || numDims <= 0) {
        printf("Error: file format (%s)\n",filename);
        MPI_Finalize();
        exit(1);
    }
    if (global_datapoint_count < processCount) {
        printf("Error: number of data points must be larger than the number of MPI processes.\n");
        MPI_Finalize();
        exit(1);
    }

    divd = global_datapoint_count / processCount;
    rem  = global_datapoint_count % processCount;
    len  = (rank < rem) ? rank*(divd+1) : rank*divd + rem;
    disp = 2 * sizeof(int) + len * numDims * sizeof(double);

    //Calculate the local size
    local_datapoint_count = (rank < rem) ? divd+1 : divd;

    //allocate a matrix to store the datapoints
    local_datapoints = allocMatrix(local_datapoint_count, numDims, 1);


    MPI_Type_contiguous(local_datapoint_count, MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File_set_view(fh, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    MPI_File_read_all(fh, local_datapoints[0], local_datapoint_count*numDims, MPI_DOUBLE, &status);
    MPI_Type_free(&filetype);
    MPI_File_close(&fh);
    /*
    printf("numObjs %d\n", global_datapoint_count);
    printf("numCoords %d\n", numDims);
    for(int i = 0; i < local_datapoint_count; i++){
        printf("Rank %d has [%.8lf, %.8lf, %.8lf]\n",rank, local_datapoints[i][0], local_datapoints[i][1], local_datapoints[i][2]);
    }
    */
}

double calcDistance(int pointIndex, int centroidIndex) {
    double sum = 0.0;
    for(int i = 0; i < numDims; i++) {
        sum = sum + pow(local_datapoints[pointIndex][i] - global_centroids[centroidIndex][i], 2);
    }
    return sqrt(sum);
}

void sumDatapointAndNewCentroid(pointIndex, centroidIndex) {
    for(int i = 0; i < numDims; i++) {
        local_newCentroids[centroidIndex][i] += local_datapoints[pointIndex][i];
    }
}

void calcNewCentroids() {

    MPI_Allreduce(local_newCentroids[0], global_centroids[0], K*numDims, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_newCentroidsSize, global_centroidsSize, K, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //printf("Rank %d has centroidsSize [%d, %d]\n",rank, global_centroidsSize[0], global_centroidsSize[1]);

    for(int k = 0; k < K; k++) {
        for(int d = 0; d < numDims; d++) {
            global_centroids[k][d] = global_centroids[k][d] / global_centroidsSize[k];
            local_newCentroids[k][d] = 0.0;
        }
        local_newCentroidsSize[k] = 0;
    }
}

void kMeans() {
    float delta = 0; // Tracks how many datapoints changed clusters each iteration
    float totalNumChanged = 0;
    double minDistance;
    double distance = 0;
    unsigned char curCluster;
    int itercount = 0;
    //for(int i = 0; i < 1; i++) { // nubmer of iterations
    do {
        delta = 0.0; // Tracks how many datapoints changed clusters each iteration
        totalNumChanged = 0;
        for(int pointIndex = 0; pointIndex < local_datapoint_count; pointIndex++) { // iterate through datapoints
            minDistance = LONG_MAX;

            for(int k = 0; k < K; k++) { // iterate throught cluster centroids
                distance = calcDistance(pointIndex, k);

                if(distance < minDistance) {

                    minDistance = distance;
                    curCluster = k;
                }
            }
            
            if(local_membership[pointIndex] != curCluster) {
                delta += 1.0;
                local_membership[pointIndex] = curCluster;
            }
            sumDatapointAndNewCentroid(pointIndex, curCluster);
            local_newCentroidsSize[curCluster] = local_newCentroidsSize[curCluster] + 1;
        }

        // calculate new cluster centers
        calcNewCentroids();
        //printCentroids();
        MPI_Allreduce(&delta, &totalNumChanged, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        delta = totalNumChanged / global_datapoint_count;

        //printf("Datapoints that changed clusters: %f\n", delta);
    } while((delta > threshold) && itercount++ < 500);
}
int mpi_write(char *filename)
{
    int        divd, rem, len, err;
    int        i, j, k;
    char       outFileName[1024], fs_type[32], str[32], *delim;
    MPI_File   fh;
    MPI_Status status;

    remove(filename);
    err = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    if (err != MPI_SUCCESS) {
        char errstr[MPI_MAX_ERROR_STRING];
        int  errlen;
        MPI_Error_string(err, errstr, &errlen);
        printf("Error at opening file %s (%s)\n", filename,errstr);
        MPI_Finalize();
        exit(1);
    }

    
    divd = global_datapoint_count / processCount;
    rem  = global_datapoint_count % processCount;
    len  = (rank < rem) ? rank*(divd+1) : rank*divd + rem;

    MPI_Offset disp = (len ) * sizeof(int);
    MPI_File_seek(fh, disp, MPI_SEEK_SET);
    MPI_File_write(fh, local_membership, local_datapoint_count, MPI_INT, &status);
    
    MPI_File_close(&fh);
    
    return 1;
}


void printMemberships(){
    int imDone = 0;
    int prevProcessDone = 0;
    if(rank == 0){
        for(int i = 0; i < local_datapoint_count; i++){
            //printf("Rank %d -> %d\n",rank, local_membership[i]);
            printf("%d",rank);
        }
        
        //printf("(%d)\n",rank);
        MPI_Send(&imDone, 1, MPI_INT, rank+1, tag, MPI_COMM_WORLD);
    }else{
        
        MPI_Recv(&prevProcessDone, 1, MPI_INT, rank-1, tag, MPI_COMM_WORLD, &status);
        //printf("(%d)\n",rank);
        for(int i = 0; i < local_datapoint_count; i++){
            //printf("Rank %d -> %d\n",rank, local_membership[i]);
            printf("%d",rank);
        }
        if(rank != processCount-1){
            MPI_Send(&prevProcessDone, 1, MPI_INT, rank+1, tag, MPI_COMM_WORLD);
        }
    }
}

void initVars() {
    local_membership = (int *) calloc(local_datapoint_count, sizeof(int));
    global_centroids    = allocMatrix(K, numDims, 1);
    local_newCentroids = allocMatrix(K, numDims, 1);
    local_newCentroidsSize = (int *) calloc(K, sizeof(int));
    global_centroidsSize = (int *) calloc(K, sizeof(int));
    if(rank == 0){
        for(int k = 0; k < K; k++) {
            for(int i = 0; i < numDims; i++) {
                global_centroids[k][i] = local_datapoints[k][i];
            }
        }
    }
    MPI_Bcast(global_centroids[0], K*numDims, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*
    for(int i = 0; i < K; i++){
        printf("Rank %d has centroids(%d) [%.8lf, %.8lf, %.8lf]\n",rank,i, global_centroids[i][0], global_centroids[i][1], global_centroids[i][2]);
    }
    */
}




int main(int argc, char *argv[])
{
    int rc;
    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &processCount);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    K = 2;
    threshold = 0.001;

    mpi_read("data.bin");
    initVars();
    MPI_Barrier(MPI_COMM_WORLD);
    double startTime = MPI_Wtime();
    kMeans();
    double endTime = MPI_Wtime();
    double global_startTime;
    double global_endTime;

    MPI_Reduce(&startTime, &global_startTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&endTime, &global_endTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if(rank == 0){
        printf("Parallel processes took %f\n", global_endTime-global_startTime);
    }
    //printMemberships();
    mpi_write("memberships_parallel.bin");
    /*
    double** local_datapoints;  // The datapoints to cluster
    unsigned char* local_membership;      // Which cluster each datapoint belongs to
    double** local_newCentroids; // Sum of the dimensions of datapoints belonging to each cluster
    int* local_newCentroidsSize;  // Number of datapoints belonging to each cluster
    double** global_centroids;    // Means of each cluster
    int* global_centroidsSize;
    */

    free(local_datapoints);
    free(local_membership);
    free(local_newCentroidsSize);
    free(local_newCentroids);
    free(global_centroids);
    free(global_centroidsSize);

    MPI_Finalize();
    //printf("Total sum is %d \n", totalSum);
    //printf("process %d just finished \n", rank);
    exit(0);

}