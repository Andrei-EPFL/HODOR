#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

// int MPI_Get_library_version(char *version, int *resultlen)



int main() {
   /* my first program in C */
   printf("Hello, World! \n");

   int l = 0;
   int n = 4000;
   char version[400];
//    version = (char*)malloc(n*sizeof(char));
   
   printf("before %d\n", l);

   MPI_Get_library_version(version, &l);
   printf("after %d\n", l);

    for (int i = 0; i < l; ++i) {
        printf("%c", version[i]);
    }

   return 0;
}