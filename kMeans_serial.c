#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

double** allocMatrix(int rows,int cols)
{
    int i;
    int header= rows *sizeof(double *);
    int data=rows*cols*sizeof(double);
    double ** rowptr=(double **)malloc(header+data);
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

    char line[1024];
    char firstLine[1024];
    fgets(firstLine, 1024, stream);
    int size = atoi(getfield(strdup(firstLine), 1));
    int dims = atoi(getfield(strdup(firstLine), 2));

    
    double** arr = allocMatrix(1000,3);
    arr[0][0] = 1.01;
    arr[0][1] = 1;
    arr[0][2] = 1;
    int row = 0;
    while (fgets(line, 1024, stream))
    {
        //printf("Datapoint %d \n",counter);
        for(int column = 0; column < dims; column++){
            char* tmp = strdup(line);
            double value;
            sscanf(getfield(tmp, column+1), "%lf", &value);
            arr[row][column] = value;
            free(tmp);
        }
        row++;
    }
    for(int i = 0; i < size; i++){
        printf("%.8f,%.8f,%.8f\n", arr[i][0], arr[i][1], arr[i][2]);
    }

}


int main(int argc, char **argv)
{
    readFile();
    return 0;
}