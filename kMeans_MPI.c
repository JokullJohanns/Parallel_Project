#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


double** local_datapoints;  // The datapoints to cluster
unsigned char* local_membership;      // Which cluster each datapoint belongs to
double** local_newCentroids; // Sum of the dimensions of datapoints belonging to each cluster
int* local_newCentroidsSize;  // Number of datapoints belonging to each cluster
double** global_centroids;    // Means of each cluster


int global_datapoint_count, processCount,Ip,rank;
int local_datapoint_count;
int numDims = 0;
int K = 2;

float** allocMatrix(int rows, int cols, int useCalloc)
{
    int i;
    int header= rows *sizeof(float *);
    int data=rows*cols*sizeof(float);
    float ** rowptr;
    if(useCalloc == 1) {
        rowptr=(float **)malloc(header+data);
    }
    else {
        rowptr=(float **)malloc(header+data);    
    }
    if(rowptr==NULL){
        return NULL;
    }
    float * buf=(float*)(rowptr+rows);
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
    disp = 2 * sizeof(int) + len * numDims * sizeof(float);

    //Calculate the local size
    local_datapoint_count = (rank < rem) ? divd+1 : divd;

    //allocate a matrix to store the datapoints
    float **local_datapoints = allocMatrix(local_datapoint_count, numDims, 1);


    MPI_Type_contiguous(local_datapoint_count, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File_set_view(fh, disp, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    MPI_File_read_all(fh, local_datapoints[0], local_datapoint_count*numDims, MPI_FLOAT, &status);
    MPI_Type_free(&filetype);
    MPI_File_close(&fh);

    printf("numObjs %d\n", global_datapoint_count);
    printf("numCoords %d\n", numDims);
    for(int i = 0; i < local_datapoint_count; i++){
        printf("Rank %d has [%f, %f, %f]\n",rank, local_datapoints[i][0], local_datapoints[i][1], local_datapoints[i][2]);
    }

}


int main(int argc, char *argv[])
{
    MPI_Status status;
    int rc;
    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &processCount);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    mpi_read("data.bin");
    MPI_Finalize();
    //printf("Total sum is %d \n", totalSum);
    //printf("process %d just finished \n", rank);
    exit(0);

}