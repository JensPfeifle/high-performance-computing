#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    float sum = 0;
    #pragma omp parallel for reduction( +:sum)
    for (int n=0; n<N;n++)
    {
        sum = sum + n;
    }
    printf("sum = %f \n",sum);
}