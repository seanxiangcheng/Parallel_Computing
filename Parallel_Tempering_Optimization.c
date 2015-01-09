/*   mpicc Program.c -lgsl -lgslcblas -lm -o program.x
 * 
 * Keywords: Open MPI, Parallel Tempering, Constraint Optimization, 
 *           Complex Networks, Jamming
 * Open MPI Implementation of Parallel Tempering in HN3 with l = 1 
 * The probelm is a NP-hard Constraint Optimization Problem.
 * The goal is to find the group state / state with the most number of particles in the complex network
 * It is designed to run any number of copies of the systems at different temperature.
 * The suggested largest system size to run is 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <term.h>
#include <ncurses.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#define MAX_SWAPS 1000000
#define MAX_SOL 256
#define CHECK_EVERY_MC 5
#define MAX_MU 7.0
#define MIN_MU 1.2
#define ROOT 0
#define TAG 0

int L = 256, K = 8, S_Count = 0, Nsolutions = 10;
int Constraint = 5; // usually it is a L-dim vector; in this case, it is the same for all, so we set it as a constant; degree + l + 1
int Site_Const = 3; // the constraint is site_const + the neighbor particles of a certain site < constraint; site_const = degree
int Check_Point_Steps, Con_Exe = 0, Max_State_World=0;
int Proc_Num = 1, Rank = 0, *Current_States;
double *Mus, Waste_Steps = 0.0;

/* 	integer power of an integer (with 0^0 = 1) 	*/
int Intpow(int base, int power){
  int x=1, i;
  for(i=0; i<power; i++) {x *= base;}
  return(x);
}


/* 	integer log2 of an integer	*/
int Intlog2(int num){
  int power=0;
  while(num>1){
    num=num/2;
    power++;
  }
  return(power);
}

int *Findij(int n){
  static int ij[2];
  ij[0]=1;
  ij[1]=0;
  while(n%2==0) {
    n=n/2;
    ij[0]++;
  }
  ij[1]=(n-1)/2;
  return(ij);
}

/*	Find the neighbors of site n in a network of length L with Periodic Boundary Conditions	  */
/*	The input is 1, 2, ...,L; output is a pointer to an array of 3 */
/*	e.g. Findneighbor(1, 8), returns *neighbors = [2, 3, 8]	       */
int * FindNeighbors(n){ 
  	static int neighbors[3];
  	int *ij, halfL=L/2;
  	if(n<1 || n>L){
    		fprintf(stderr,"\nError:  Cannot find the neighbor of site %d in HN3 with length %d! \n", n, L);
  	}

  	if(n==1){
    		neighbors[0]=2;
    		neighbors[1]=3;
    		neighbors[2]=L;
  	}
  	else if(n==halfL){
    		neighbors[0]=n-1;
    		neighbors[1]=n+1;
    		neighbors[2]=L;
  	}
  	else if(n==L){
    		neighbors[0]=1;
    		neighbors[1]=halfL;
    		neighbors[2]=L-1;
  	}
  	else{
    		ij=Findij(n);
    		neighbors[0]=n-1;
    		neighbors[1]=n+1;
    		if(ij[1]%2==0){
      			neighbors[2]=Intpow(2,ij[0]-1)*(2*(ij[1]+1)+1);
    		}
   		else{
      			neighbors[2]=Intpow(2,ij[0]-1)*(2*(ij[1]-1)+1);
    		}
  	}
  	return(neighbors);
}

//   Initialize(xs, neighbors, solutions, n_proc, Nsolutions, mus, current_state, max_state_found, total_steps;
/*	Initialize the StateInfo, neighbors and other parameters*/
void Initialize(int xs[], int neighbors[], int solutions[],  int Nsolutions, double *mu, int *current_state, int *max_state_found){
  int i;
  int *neigh;

  for(i = 0; i < L; i++){
    xs[i] = 0;
	}
  for(i = 0; i < L; i++){
    neigh=FindNeighbors(i+1);
    neighbors[3*i]=neigh[0]-1;
    neighbors[3*i+1]=neigh[1]-1;
    neighbors[3*i+2]=neigh[2]-1;
  }
  for(i = 0; i < L*Nsolutions; i++){
      solutions[i] = 0;
  }
  *mu = (double)MAX_MU- (double)(MAX_MU- MIN_MU)* Rank / (double)(Proc_Num-1);
  *current_state = 0;
  *max_state_found =  9*L/16; 
  Max_State_World = *max_state_found;
  Check_Point_Steps = CHECK_EVERY_MC * L;
  return;
}


/* Print sites*/
void Print_Sites(int xs[]){
    int i, rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    printf("  Process # %d\n", rank);
    for(i = 0; i < L; i++){
        printf(" %-5d  :  %-5d \n", i+1, xs[i]);
    }
    return;
}


/*	Check whether to add a particle to an empty site i	*/
int Add_or_Not(int xs[], int i, int neighbors[]){
	int j, neighloc_begin, neighloc_end;
	int neig_num=0, site_1_particle;
	if(xs[i] != 0){/*Could Delete it*/
		fprintf(stderr, "\n \n Error: The site tried to add particle has 0 particle: site[%d]=%d\n", i, xs[i]); 
    Print_Sites(xs);
    exit(0);	
	}
	neighloc_begin=3*i;
	neighloc_end=neighloc_begin+3;
	/*	First Check its own neighbors	*/
	for(j = neighloc_begin; j < neighloc_end; j++){
		neig_num += xs[neighbors[j]];
    if(xs[neighbors[j]] == 1){
			site_1_particle = neighbors[j];
    }
	}
	if(neig_num > 1){
		return(0);
	}
	else if(neig_num==0){
		return(1);
	}
	else{
		neighloc_begin = site_1_particle * 3;
		neighloc_end = neighloc_begin + 3;
		neig_num=0;
		for(j = neighloc_begin; j < neighloc_end; j++){
			neig_num += xs[neighbors[j]];
		}
		if(neig_num == 0){return(1);}
 		else{return(0);}
	}
	return(0);
}

/*	Check whether to switch a particle with an empty site i	*/
int Exchange_or_not(int xs[], int i, int j, int neighbors[]){
	int neighloc_begin, neighloc_end, l;
	int neig_num = 0, site_1_particle = 0;
	neighloc_begin = 3*j;
	neighloc_end = neighloc_begin + 3;
	/*	First Check its own neighbors	*/
	for(l = neighloc_begin; l < neighloc_end; l++){
		neig_num += xs[neighbors[l]];
	}
	if(neig_num>2){
		return(0);
	}
	else if(neig_num==1){
		return(1);
	}
	else if(neig_num==2){
		for(l = neighloc_begin; l < neighloc_end; l++){
			if(xs[neighbors[l]]==1 && neighbors[l] != i ){
				site_1_particle = neighbors[l];
				break;
			}
		}
		neighloc_begin = site_1_particle*3;
		neighloc_end = neighloc_begin+3;
		neig_num = 0;
		for(l = neighloc_begin; l < neighloc_end; l++){
			if(neighbors[l] != i){
				neig_num += xs[neighbors[l]];
			}
		}
		if(neig_num==0){return(1);}
 		else{return(0);}
	}
	else{
		fprintf(stderr, "\n  Error: neig_num is not  1, 2, 3 but = %d!; i: site[%d]=%d; j: site[%d]=%d, \n", neig_num, i, xs[i], j, xs[j]); 
        	exit(0);		
	}
	return(0);
}
//  Update(xs, neighbors, i, exp_mu_inv, current_state, max_state_found)
int Update(int xs[], int neighbors[], double exp_mu_inv, int *current_state, int *max_state_found){
	int i, j, l, return_value = 0;
	i = rand()%L;
	if(xs[i] == 1)
  {
			if((double)rand()/(double)RAND_MAX < exp_mu_inv)
      {
				xs[i] = 0;
				*current_state = *current_state - 1;
      }
	}
	else
  {
		if(Add_or_Not(xs, i, neighbors))
    {
			xs[i] = 1;
      *current_state = *current_state + 1;
			if(*current_state >= *max_state_found)	return_value = 1;
		}
		else  Waste_Steps += 1.0;
  }
	
	for(j = 0; j < 6; j++)
  {
    i = rand()%L;
    if(xs[i] == 1) break;
  }
	if(xs[i])
  {
    for(l = 0; l < 6; l++)
    {
      j = neighbors[3*i + rand()%3];
      if(xs[j] == 0) break;
    }
    if(xs[j]==0){
      if(Exchange_or_not(xs, i, j, neighbors))
      {
        xs[i] = 0;
        xs[j] = 1;
      }
    }
	}
  return(return_value);
}

int Reset0(int solutions[], double sol_time[]){
  int i, S_Count_L = S_Count * L;
  for(i = 0; i < S_Count_L; i++){
    solutions[i] = 0;
  }
  for(i = 0; i< S_Count; i++){
    sol_time[i] = -1.0;
  }
  return(0);
}

int Config2Solution(int xs[], int solutions[], double sol_time[], double start, double end){
  int i, j, diff_count = 0, same_solu = -1, S_Count_L;
  if(S_Count >= Nsolutions) return 0;
  end = MPI_Wtime();
  if(S_Count == 0){
      S_Count = 1;
      for(i = 0; i < L; i++){
        solutions[i] = xs[i];
      }
      sol_time[0] = end - start;
      printf("     First solution in Rank #%d! \n", Rank);
  }
  else{
      for(i = 0; i < S_Count; i++){
          for(j = 0; j < L; j++){
            if(xs[j] != solutions[i*L + j]) {
              diff_count++;
              break;
            }
            else if(j == L-1){
              same_solu = i;
              break;
            }
          }
          if(same_solu != -1) break;
      }
  }
  if(diff_count == S_Count){
    S_Count_L = S_Count * L;
    for(i = 0; i < L; i++){
      solutions[S_Count_L + i] = xs[i];
    }
    sol_time[S_Count] = end - start; 
    S_Count++;
  }
  return(0);
}

int Swaps_or_Not(int cs1, int cs2, double mu1, double mu2)
{
    if((double)rand()/(double)RAND_MAX < exp((mu2 - mu1)*(cs2-cs1)))
        return(1);
    else return(0);
}



int Check_Config(int xs[], int neighbors[], int cs){
  int i, neig_num, cs_actual=0;
  for(i = 0; i < L; i++){
      cs_actual +=xs[i];
      neig_num = xs[i]*3 + xs[neighbors[3*i]] + xs[neighbors[3*i + 1]] + xs[neighbors[3*i + 2]];
      if(neig_num > 4){
          printf(" Constraint violated! x[%d] = %d; neighbors: x[%d] = %d, x[%d] = %d, x[%d] = %d\n", i, xs[i],neighbors[3*i], xs[neighbors[3*i]], neighbors[3*i + 1], xs[neighbors[3*i + 1]], neighbors[3*i + 2], xs[neighbors[3*i + 2]]);
          exit(0);
      }
  }
  if(cs_actual != cs)
  {
      printf(" Actual current state = %d;\n in program, current state=%d\n", cs_actual, cs);
      exit(0);  
  }
  return(0);
}

/*	Function to read input parameters	*/
void Commandlineparse(int argc, char **argv, int *L,int *K, int *name, int *nsolution, int *seed){
  	int i;
  	for (i = 1; i < argc; i++){  //Start at i = 1 to skip the command name.
    		if (argv[i][0] == '-'){
      			switch (argv[i][1]){
      				case 'L':       *L = atoi(argv[++i]);
        			break;
      				case 'k':       *K = atoi(argv[++i]);
        			break;
      				case 'n':       *name=atoi(argv[++i]);
        			break;
              case 'm':       *nsolution = atoi(argv[++i]);
              break;
              case 's':       *seed = atoi(argv[++i]);
              break;
      				case 'h':
        			fprintf(stderr,"\n Available options: \n\
              -L = system size, default 256; it has to be 2^k, k=2,3,...;\n\
              -k = system size is 2^k; default 8; \n\
              -m = number of solutions to output (It has to be less than system size)\n\
              -n = name: the last digit of the output file name, default 0;\n\
              e.g. input of 'n' is 0, filename is solution_k8_r0 \n\
              which is the solution for system size 2^8\n\
              \n");
        			exit(0);
      			}//switch
    		}//if
  	}//for
}



int main(int argc, char *argv[]){
  int i, j, l, temp, seed = 1;
  int  name = 0;
  char s_filename[30] = "HN3Sol_mpi_";
  char buf[6];
  MPI_Status mpi_stat;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
  MPI_Comm_size(MPI_COMM_WORLD, &Proc_Num);
  if(Rank == 0)
  {
    printf("\n    ****** Parallel Tempering (MPI version) Running ****** \n");
    printf("      argc: %d\n      ",argc);
    for(j = 0; j<argc; j++)  printf("%s ",argv[j]);
    printf("\n\n");
  }
  Commandlineparse(argc, argv, &L, &K, &name, &Nsolutions, &seed);
  
  if(K != 8 ){
   		L=Intpow(2,K);
  }
  else if(L != 256){
      K=Intlog2(L);
  }
  if(Rank == 0)
  {
    if(argc < 2) printf("    --- Default settings due to no arguments --- \n     (\".\\ProgramName -h\" for all available options) \n\n");
    printf("     System size: L = %d; \n", L);
    printf("     Max MC Steps: %.4e \n", (double)(MAX_SWAPS*L)*(double)CHECK_EVERY_MC);
    if(Nsolutions > MAX_SOL)
    {
      printf("     Number of solutions to output %d is too big; \n     it should be less than %d", Nsolutions, MAX_SOL);
      exit(0);
    }
  }
  srand((unsigned int)(Rank*3 + Rank*Rank + 87 + Rank%2 + Rank%3 + name*73));
  int  total_swaps = 0, check_swaps = MAX_SWAPS / 10000, extra_swaps = 0; // check_swaps is the extra random steps to see whether it is the real solution
  double  sol_time[Nsolutions], total_time = 0.0, initial_time = 0.0;
  
  //int max_swaps_abs = MAX_STEPS * L;
  int xs[L], xs_send[L], *xs_copy, swap_rank[2], *swap_l0r; // number of sites / variables
  int neighbors[3*L]; // locations of constraint sites
  int solutions[Nsolutions * L], *solutions_proc_copy, *sol_nums;
  int current_state, max_state_found;
  double mu, exp_mu_inv, t0, t1;
  int in_steps = 0;
  if(check_swaps < 1) check_swaps = 1;
  if(Rank == ROOT)
  {
      Current_States = (int *)malloc(Proc_Num * sizeof(int));
      Mus = (double *)malloc(Proc_Num*sizeof(double));
      sol_nums = (int *)malloc(Proc_Num*sizeof(int));
      solutions_proc_copy = (int *)malloc(Nsolutions *L *sizeof(int));
      xs_copy = (int *)malloc(L*sizeof(int));
      swap_l0r = (int *)malloc(2*Proc_Num*sizeof(int));
      for(i = 0; i < Proc_Num; i++)
      {
        Current_States[i] = 0;
        Mus[i] = (double)MAX_MU - (double)(MAX_MU-MIN_MU)*i/(double)(Proc_Num-1);      
      }
      sprintf(buf, "%d", K);
      strcat(s_filename, buf);
      strcat(s_filename, "_r");
      sprintf(buf,"%d", name);
      strcat(s_filename, buf);
  }
  t0 = MPI_Wtime();
  Initialize(xs, neighbors, solutions, Nsolutions, &mu, &current_state, &max_state_found);
  //printf("     Initialization DONE: Proc #%d. \n", Rank);
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  initial_time = t1 - t0;
  exp_mu_inv = exp(-mu);
  while (total_swaps < MAX_SWAPS)
  {
      total_swaps++;
      in_steps = 0;
      while(in_steps < Check_Point_Steps)
      {
          if(Update(xs, neighbors, exp_mu_inv, &current_state, &max_state_found)) //if current_state >= max_state_found
          {
              if(current_state > max_state_found)
              {
                  max_state_found = current_state;
                  Max_State_World = max_state_found;
                  Reset0(solutions, sol_time);
                  S_Count = 0;
                  extra_swaps = 0;
                  Config2Solution(xs, solutions, sol_time, t0, t1);
              }
              else
              {
                  if(S_Count < Nsolutions)
                  Config2Solution(xs, solutions, sol_time, t0, t1);
              }
          }
          in_steps++;
          //Check_Config(xs, neighbors, current_state);
      }
      //printf("Proc #%d: CS=%d : Nsol = %d\n", Rank, current_state, S_Count); //del
      
      if(extra_swaps > check_swaps)
      {
          if(Rank == ROOT)
            printf("     Number of solutions reached %d \n", Nsolutions);
          break;
      }
      MPI_Allreduce(&max_state_found, &Max_State_World, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if(Max_State_World > max_state_found)
      {
          max_state_found = Max_State_World;
          Reset0(solutions, sol_time);
          S_Count = 0;
          extra_swaps = 0;
      }
      MPI_Gather(&S_Count, 1, MPI_INT, sol_nums, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
      MPI_Gather(&current_state, 1, MPI_INT, Current_States, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
      if(Rank == ROOT)
      {
        for(i = 0; i < Proc_Num - 1; i++)
        {
            // Receive the solutions from  all processors except root
            MPI_Recv(solutions_proc_copy, sol_nums[i+1]*L, MPI_INT, i+1, TAG, MPI_COMM_WORLD, &mpi_stat);
            // Merge the solutions into Nsolutions 
            for(j = 0; j < sol_nums[i+1]; j++)
            {
                if(S_Count < Nsolutions)
                {
                    for(l = 0; l < L; l++)  xs_copy[l] = solutions_proc_copy[j*L + l];
                    Config2Solution(xs_copy, solutions, sol_time, t0, t1);
                }
            }
            //printf("  sol merging for Proc #%d in Root #%d: N sol = %d\n", i+1, Rank, S_Count); //del
        }
        MPI_Bcast(&S_Count, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(solutions, S_Count*L, MPI_INT, ROOT, MPI_COMM_WORLD);
        // Check to see whether to exchange with some adjacent processor
        for(i=0; i < Proc_Num-1; i++)
        {
          temp = Swaps_or_Not(Current_States[i], Current_States[i+1], Mus[i], Mus[i+1]);
          swap_l0r[2*i + 1] = temp;
          swap_l0r[2*(i+1)] = temp;
        }
        swap_l0r[0] = 0;
        swap_l0r[2*Proc_Num - 1] = 0;
        MPI_Scatter(swap_l0r, 2, MPI_INT, swap_rank, 2, MPI_INT, ROOT, MPI_COMM_WORLD);
      }
      else
      {
            MPI_Send(solutions, S_Count*L, MPI_INT, ROOT, TAG, MPI_COMM_WORLD);
            MPI_Bcast(&S_Count, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(solutions, S_Count*L, MPI_INT, ROOT, MPI_COMM_WORLD);
            MPI_Scatter(swap_l0r, 2, MPI_INT, swap_rank, 2, MPI_INT, ROOT, MPI_COMM_WORLD);
      }
      
      // Exchange configurations based on the swap_rank 
      for(i=0; i<L; i++) xs_send[i] = xs[i]; // a copy of its configuration
      if(Rank%2 == 0 && Rank != (Proc_Num-1) && swap_rank[1]==1)
      {
          MPI_Send(xs_send, L, MPI_INT, Rank+1, TAG, MPI_COMM_WORLD);
          MPI_Recv(xs, L, MPI_INT, Rank+1, TAG, MPI_COMM_WORLD, &mpi_stat);
      }
      else if(Rank%2 == 1 && Rank != 0 && swap_rank[0]==1)
      {
          MPI_Recv(xs, L, MPI_INT, Rank-1, TAG, MPI_COMM_WORLD, &mpi_stat);
          MPI_Send(xs_send, L, MPI_INT, Rank-1, TAG, MPI_COMM_WORLD);
      }
          
      if(Rank%2 == 0 && Rank != 0 && swap_rank[0]==1)
      {
          MPI_Send(xs_send, L, MPI_INT, Rank-1, TAG, MPI_COMM_WORLD);
          MPI_Recv(xs, L, MPI_INT, Rank-1, TAG, MPI_COMM_WORLD, &mpi_stat);
      }
      else if(Rank%2 == 1 && Rank != 0 && swap_rank[1]==1)
      {
          MPI_Recv(xs, L, MPI_INT, Rank+1, TAG, MPI_COMM_WORLD, &mpi_stat);
          MPI_Send(xs_send, L, MPI_INT, Rank+1, TAG, MPI_COMM_WORLD);
      }
      current_state = 0;
      for(i = 0; i<L; i++) current_state += xs[i];
      //printf("  Swap # %d Done: Proc # %d\n", total_swaps, Rank); //del
      if(S_Count == Nsolutions) extra_swaps++;    
    /* Check current states of all and try to exchange with its neighbor */
    // Fill the code for processes communciation and exchange with each other 
    /*                 */
      //Proc_Swaps(xs, current_state, n_proc, mus);
  }
  
  
  t1 = MPI_Wtime();
  total_time = t1 - t0; 
  if(total_swaps == MAX_SWAPS && Rank == ROOT){
     printf("     Max Swaps reached: %.2e \n", (double)MAX_SWAPS*(double)Check_Point_Steps*Proc_Num);
     printf("     Total time spent: %.6e \n", total_time);
   }
  if(Rank == ROOT)
  {
      FILE *fp;
      fp=fopen(s_filename,"w");
      fprintf(fp, "L          : %d\n", L);
      fprintf(fp, "Total time : %.6e\n", total_time);
      fprintf(fp, "Total Stpes: %.6e\n", (double)total_swaps*(double)Proc_Num*(double)Check_Point_Steps);
      fprintf(fp, "Nsolutions : %d\n", S_Count);
      fprintf(fp, "Num Proc   : %d\n", Proc_Num);
      fprintf(fp, "Total Exch : %d\n", total_swaps);
      fprintf(fp, "\nTime spent to find the %d solutions:\n", S_Count);
      for(i = 0; i < S_Count; i++)
        fprintf(fp, "%-.6e \n", sol_time[i]);
      fprintf(fp, "\nInit time  : %.6e\n", initial_time);
      fprintf(fp, "Max State  : %d\n", Max_State_World);
      fprintf(fp, "\nMus      CS\n");
      for(i = 0; i < Proc_Num; i++)
        fprintf(fp, "%-8.4f %-8d\n", Mus[i], Current_States[i]);
      fprintf(fp, "\nSoltions (one column one solution):\n");
      for(i = 0; i < L; i++)
      {  
          for(j=0; j < S_Count; j++)
              fprintf(fp, "%-4d", solutions[i + j*L]);
          fprintf(fp, "\n");
      }
      printf("     Nsolutions = %d; maxstatefound = %d; \n     current_state = %d, %d, %d, %d;\n     Total Config Exchanges: %d\n\n", S_Count, max_state_found, Current_States[0], Current_States[1], Current_States[2],Current_States[3], total_swaps);
  }
  /* Free memeory and finialize MPI program */
  if(Rank == ROOT)
  {
    free(Current_States); 
    free(Mus);
    free(sol_nums);
    free(xs_copy);
    free(swap_l0r);
    free(solutions_proc_copy);
  }
  MPI_Finalize();
  return(0);  
}
