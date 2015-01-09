/* mpicc Program.c -o Program.x
 * This program is to measure the performance of 
 * built-in MPI_Bcast() function, and compare to
 * MPI_Bcast_DIY.c
*/

#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAX_MSG_SIZE 1048576
#define REPEAT_NUM 100
#define SIZE_NUM 5
int main(int argc, char** argv) {
  int rank, i, proc_size, msg_size[] ={1, 1024, 4096, 65536, 1048576}, size_seq = 0, repeat = 0;
  double t[SIZE_NUM], t1 = 0.0, t2 = 0.0, *all_times, *time_all_sizes, total_t = 0.0;
  char buf[MAX_MSG_SIZE]; 
  const int root=0;
  char filename[50] = "lab0z_n", file_proc[8];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
  
  sprintf(file_proc, "%d", proc_size);
  strcat(filename, file_proc);
  strcat(filename, ".txt");
  
  if(rank == root) {
    for(i = 0; i < MAX_MSG_SIZE; i++){
      buf[i] = 'o';
    }
    all_times = (double *)malloc(proc_size * sizeof(double));
    time_all_sizes = (double *)malloc(SIZE_NUM * sizeof(double));
  }
  else{
    for(i = 0; i < MAX_MSG_SIZE; i++){
        buf[i] = 'x';
    }
  }
  for(i = 0; i < SIZE_NUM; i++){
    if(root == rank) time_all_sizes[i] = 0.0;
    t[i] = 0.0;
  }
  
  for(size_seq = 0; size_seq < SIZE_NUM; size_seq++){
    for(repeat = 0; repeat < REPEAT_NUM; repeat++){
      t1 = MPI_Wtime();
      MPI_Bcast(buf, msg_size[size_seq], MPI_CHAR, root, MPI_COMM_WORLD);
      t2 = MPI_Wtime();
      t[size_seq] += t2 - t1;
    }
  }
  // Gather the time information to the root process to find the total time
  for(size_seq = 0; size_seq < SIZE_NUM; size_seq++){
    MPI_Gather(&t[size_seq], 1, MPI_DOUBLE, all_times, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    if( root == rank){
      total_t = 0.0;
      for(i = 0; i < proc_size; i++){
        total_t += all_times[i];
      }
      time_all_sizes[size_seq] = total_t / ((double)REPEAT_NUM) / (double)proc_size;
    }
  }
  
  if(rank == root){
    printf(" Number of Processes: %d\n", proc_size);
    printf(" Message size (next line) && time (next next line)\n");
    for(i = 0; i < SIZE_NUM; i++){
        printf("  %-10d", msg_size[i]);
    }
    printf("\n");
    for(i = 0; i< SIZE_NUM; i++){
        printf("  %-10.4e", time_all_sizes[i]);
    }
    printf("\n");
  }
  MPI_Finalize();
  
  if(rank == root){
      FILE *fp;
      fp = fopen(filename, "w");
      fprintf(fp, " Number of Processes: %d\n", proc_size);
      fprintf(fp, " Message size (next line) && time (next next line)\n");
      for(i = 0; i < SIZE_NUM; i++){
        fprintf(fp, "  %-10d", msg_size[i]);
      }
      fprintf(fp, "\n");
      for(i = 0; i< SIZE_NUM; i++){
        fprintf(fp, "  %-10.4e", time_all_sizes[i]);
      }
      fprintf(fp, "\n");
      fclose(fp);
      printf("\n *** Info also saved in Filename \"%s\" ***\n\n", filename);
  }
  return 0;
}
