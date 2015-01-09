/* This is a code to sample the system and find proper parameters 
*  for Parallel_Tempering_Optimization_MPI.c
*  It is designed to run on single CPU
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <term.h>
#include <ncurses.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include <mpi.h>
#define MAX_STEPS 10000
#define CHECK_EVERY_MC 10

int L = 256, K = 8, S_Count = 0;
int Constraint = 5; // usually it is a L-dim vector; in this case, it is the same for all, so we set it as a constant; degree + l + 1
int Site_Const = 3; // the constraint is site_const + the neighbor particles of a certain site < constraint; site_const = degree
int Check_Point_Steps, Con_Exe = 0;
double Waste_Steps = 0.0;

gsl_rng *rnd;	// RNG: random number


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

//   Initialize(xs, neighbors, solutions, n_proc, nsolutions, mus, current_state, max_state_found, total_steps;
/*	Initialize the StateInfo, neighbors and other parameters*/
void Initialize(int xs[], int neighbors[], int solutions[], int n_proc,  int nsolutions, double mus[], int current_state[], int *max_state_found){
  int i;
  int *neigh;

  for(i = 0; i < L * n_proc; i++){
    xs[i] = 0;
	}
  for(i = 0; i < L; i++){
    neigh=FindNeighbors(i+1);
    neighbors[3*i]=neigh[0]-1;
    neighbors[3*i+1]=neigh[1]-1;
    neighbors[3*i+2]=neigh[2]-1;
  }
  for(i = 0; i < L*nsolutions; i++){
      solutions[i] = 0;
  }
  for(i = 0; i < n_proc; i++){
      mus[i] = 6.8 - 6.8 * i / n_proc;
      current_state[i] = 0;
  }
  *max_state_found =  9*L/16; 
  Check_Point_Steps = CHECK_EVERY_MC * L;
  return;
}


/* Print sites*/
void Print_Sites(int xs[], int proc){
    int i;
    printf("  Process # %d\n",proc);
    for(i = 0; i < L; i++){
        printf(" %-5d  :  %-5d \n", i+1, xs[proc*L + i]);
    }
    return;
}


/*	Check whether to add a particle to an empty site i	*/
int Add_or_Not(int xs[], int i, int proc, int neighbors[]){
	int j, neighloc_begin, neighloc_end;
	int neig_num=0, site_1_particle;
  int procL = proc * L;
	if(xs[procL + i] != 0){/*Could Delete it*/
		fprintf(stderr, "\n \n Error: The site tried to add particle has %d particle!\n Proc = %d; ProcL = %d \n", xs[procL + i], proc, procL); 
    Print_Sites(xs, proc);
    exit(0);	
	}
	neighloc_begin=3*i;
	neighloc_end=neighloc_begin+3;
	/*	First Check its own neighbors	*/
	for(j = neighloc_begin; j < neighloc_end; j++){
		neig_num += xs[procL + neighbors[j]];
    if(xs[procL + neighbors[j]] == 1){
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
			neig_num += xs[procL + neighbors[j]];
		}
		if(neig_num == 0){return(1);}
 		else{return(0);}
	}
	return(0);
}

/*	Check whether to switch a particle with an empty site i	*/
int Exchange_or_not(int xs[], int i, int j, int proc, int neighbors[]){
	int neighloc_begin, neighloc_end, l;
	int neig_num = 0, site_1_particle = 0;
  int procL = proc * L;
	neighloc_begin = 3*j;
	neighloc_end = neighloc_begin + 3;
	/*	First Check its own neighbors	*/
	for(l = neighloc_begin; l < neighloc_end; l++){
		neig_num += xs[procL + neighbors[l]];
	}
	if(neig_num>2){
		return(0);
	}
	else if(neig_num==1){
		return(1);
	}
	else if(neig_num==2){
		for(l = neighloc_begin; l < neighloc_end; l++){
			if(xs[procL + neighbors[l]]==1 && neighbors[l] != i ){
				site_1_particle = neighbors[l];
				break;
			}
		}
		neighloc_begin = site_1_particle*3;
		neighloc_end = neighloc_begin+3;
		neig_num = 0;
		for(l = neighloc_begin; l < neighloc_end; l++){
			if(neighbors[l] != i){
				neig_num += xs[procL + neighbors[l]];
			}
		}
		if(neig_num==0){return(1);}
 		else{return(0);}
	}
	else{
		fprintf(stderr, "\n  Error: neig_num is not  1, 2, 3 but = %d!; i: site[%d]=%d; j: site[%d]=%d, \n", neig_num, i, xs[procL + i], j, xs[procL + j]); 
        	exit(0);		
	}
	return(0);
}
//  Update(xs, neighbors, i, exp_mu_inv, current_state, max_state_found)

int Update(int xs[], int neighbors[], int proc, double exp_mu_inv, int current_state[], int max_state_found){
	int i, j, l, return_value = 0;
  int procL = proc*L;
	i=(int)(L*gsl_rng_uniform(rnd));
	if(xs[i + procL] == 1){
			if(gsl_rng_uniform(rnd) < exp_mu_inv){
				xs[i + procL] = 0;
				current_state[proc]--;
      }
	}
	else{
		if(Add_or_Not(xs, i, proc, neighbors)){
			xs[i + procL] = 1;
			current_state[proc]++;
			if(current_state[proc] >= max_state_found){
				return_value = 1;
			}
		}
		else{
      Waste_Steps += 1.0;
    }
  }
	
	for(j = 0; j < 6; j++){
    i = (int)(L*gsl_rng_uniform(rnd));
    if(xs[i + procL] == 1) break;
  }
	if(xs[i + proc * L]){
    for(l = 0; l < 6; l++){
      j = neighbors[3*i + (int)(3*gsl_rng_uniform(rnd))];
      if(xs[procL + j] == 0) break;
    }
    if(xs[ procL + j]==0){
      if(Exchange_or_not(xs, i, j, proc, neighbors)){
        xs[i + proc * L] = 0;
        xs[j + proc * L] = 1;
      }
    }
	}
  return(return_value);
}

int Reset0(int solutions[]){
  int i, S_Count_L = S_Count * L;
  for(i = 0; i < S_Count_L; i++){
    solutions[i] = 0;
  }
  return(0);
}

int Config2Solution(int xs[], int proc, int solutions[]){
  int i, j, procL = proc * L, diff_count = 0, same_solu = -1, S_Count_L;
  if(S_Count == 0){
      S_Count = 1;
      for(i = 0; i < L; i++){
        solutions[i] = xs[procL + i];
      }
  }
  else{
      for(i = 0; i < S_Count; i++){
          for(j = 0; j < L; j++){
            if(xs[procL + j] != solutions[i*L + j]) {
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
      solutions[S_Count_L + i] = xs[procL + i];
    }
    S_Count++;
  }
  
  return(0);
}

int Proc_Swaps(int xs[], int current_state[], int n_proc, double mus[]){
  int i, j, temp;
  for(i = 1; i < n_proc; i++){
    if(gsl_rng_uniform(rnd) < exp((mus[i] - mus[i-1])*(current_state[i] - current_state[i-1]))){
      for(j = 0; j < L; j++){
        temp = xs[i * L + j];
        xs[i*L + j] = xs[(i-1)*L + j];
        xs[(i-1)*L + j] = temp;
      }
      temp = current_state[i];
      current_state[i] = current_state[i-1];
      current_state[i-1] = temp;
      Con_Exe++;
    }  
  }
  return(0);
}

int Check_Config(int xs[], int neighbors[], int n_proc){
  int i, j, neig_num;
  for(i = 0; i < n_proc; i++){
    for(j = 0; j < L; j++){
      neig_num = xs[i*L + j]*3 + xs[i*L + neighbors[3*j]] + xs[i*L + neighbors[3*j + 1]] + xs[i*L + neighbors[3*j + 2]];
      if(neig_num > 4){
          printf(" Constraint violated! x[%d] = %d; neighbors: x[%d] = %d, x[%d] = %d, x[%d] = %d", j+1, xs[i*L + j],neighbors[3*j], xs[i*L + neighbors[3*j]], neighbors[3*j + 1], xs[i*L + neighbors[3*j + 1]], neighbors[3*j + 2],xs[i*L + neighbors[3*j + 2]]);
          exit(0);
      }
    }
  }
    
  return(0);
}

/*	Function to read input parameters	*/
void Commandlineparse(int argc, char **argv, int *L,int *K, int *n_proc, int *seed, int *name, int *nsolution){
  	int i;
  	*seed = getpid() * 7;
  	for (i = 1; i < argc; i++){  //Start at i = 1 to skip the command name.
    		if (argv[i][0] == '-'){
      			switch (argv[i][1]){
      				case 'L':       *L = atoi(argv[++i]);
        			break;
      				case 's':       *seed = atoi(argv[++i]);
        			break;
      				case 'p':       { 
                                *n_proc = atoi(argv[++i]);
                                if (*n_proc < 2){ 
                                  printf("The number of processes has to be >= 2. Otherwise, this is not Parallel Tempering.\n");
                                  exit(0);
                                }
                                break;
                              }
      				case 'k':       *K = atoi(argv[++i]);
        			break;
      				case 'n':       *name=atoi(argv[++i]);
        			break;
              case 'm':       *nsolution = atoi(argv[++i]);
              break;
      				case 'h':
        			fprintf(stderr,"\n Available options: \n\
              -L = system size, default 256; it has to be 2^k, k=2,3,...;\n\
              -k = system size is 2^k; default 8; \n\
              -p = number of processes (default 4, max = 32 for memory of 8G);\n\
              -s = seed: random seed  (default=pid); \n\
              -m = number of solutions to output (It has to be less than system size)\n\
              -n = name: the last 3 digits of the output file name, default 108;\n\
              e.g. input of 'n' is 108, filename is solution108 \n\
              which is the solution for system size 2^8\n\
              \n");
        			exit(0);
      			}//switch
    		}//if
  	}//for
}



int main(int argc, char *argv[]){
  int i, j;
  int n_proc = 4, name = 108, seed, nsolutions = 8;
  char s_filename[20] = "HN3Sol";
  char buf[6];
  printf("\n    ************ Parallel Tempering Running ************ \n");
  printf("\n      argc: %d\n      ",argc);
  for(j = 0; j<argc; j++){
    printf("%s ",argv[j]);
  }
  printf("\n\n");
  Commandlineparse(argc, argv, &L, &K, &n_proc, &seed, &name, &nsolutions);
  
  rnd =  gsl_rng_alloc(gsl_rng_mt19937);
 	gsl_rng_set(rnd, (unsigned long int)seed);
  for(i = 0; i < 20; i++){
    gsl_rng_uniform(rnd);
  }
  if(K != 8 ){
   		L=Intpow(2,K);
  }
  else if(L != 256){
      K=Intlog2(L);
  }
  if(argc <2) printf("    --- Default settings due to no arguments --- \n     (\".\\ProgramName -h\" for all available options) \n\n");
  printf("     System size: L = %d; \n", L);
  printf("     Number of processes: n_proc = %d; \n     Max MC Steps: %.4e \n", n_proc, (double)(MAX_STEPS*L)*(double)CHECK_EVERY_MC);
  if(nsolutions > L){
      printf("     Number of solutions to output %d is too big; \n     it should be less than system size %d", nsolutions, L);
      exit(0);
  }
  int  check_swaps = MAX_STEPS / 10, extra_swaps = 0;
  //int max_swaps_abs = MAX_STEPS * L;
  int xs[L * n_proc]; // number of sites / variables
  int neighbors[3*L]; // locations of constraint sites
  
  int solutions[nsolutions * L];
  int current_state[n_proc], max_state_found;
  double mus[n_proc], exp_mu_inv;
  int  total_swaps = 0;  
  int in_steps = 0;
  
  sprintf(buf,"%d", name);
 	strcat(s_filename, buf);
  
  Initialize(xs, neighbors, solutions, n_proc, nsolutions, mus, current_state, &max_state_found);
  printf("     Initialization DONE. \n");
  while (total_swaps < MAX_STEPS){
    total_swaps++;
    //printf(" swaps: %d\n",total_swaps); //del
    for(i = 0; i < n_proc; i++){
        //printf("         proc: %d \n", i); //del
        exp_mu_inv = exp(-mus[i]);
        in_steps = 0;
        while(in_steps < Check_Point_Steps){
          if(Update(xs, neighbors, i, exp_mu_inv, current_state, max_state_found)){
              if(current_state[i] > max_state_found){
                  max_state_found = current_state[i];
                  Reset0(solutions);
                  S_Count = 0;
                  extra_swaps = 0;
                  Config2Solution(xs, i, solutions);
              }
              else{
                    if(S_Count == nsolutions){
                      extra_swaps++;
                      break;
                    }
                    Config2Solution(xs, i, solutions);
              }
          }
          Check_Config(xs, neighbors, n_proc);
          in_steps++;
        }
        if(S_Count == nsolutions){
            extra_swaps++;
            break;
        }
    }
    if(extra_swaps > check_swaps){
      printf("     Number of solutions reached %d \n", nsolutions);
      break;
    }
    Check_Config(xs, neighbors, n_proc);
    Proc_Swaps(xs, current_state, n_proc, mus);
  }
  
  if(total_swaps == MAX_STEPS) printf("     Max Swaps reached: %.2e \n", (double)MAX_STEPS*(double)Check_Point_Steps*n_proc);
  
  FILE *fp;
  fp=fopen(s_filename,"w");
  if(total_swaps == MAX_STEPS) fprintf(fp, " Max Swaps reached: %.2e", (double)MAX_STEPS*(double)Check_Point_Steps*n_proc);
  fprintf(fp, "Total Steps: %.4e; Number of solutions: %d; MaxState: %d \n; Total Config Exchanges: %d",(double)total_swaps*(double)n_proc*(double)Check_Point_Steps, S_Count, max_state_found, Con_Exe);
  for(i = 0; i < L; i++){  
    for(j=0; j < S_Count; j++){
      fprintf(fp, "%-4d", solutions[i + j*L]);
    }
    fprintf(fp, "\n");
  }
  Check_Config(xs, neighbors, n_proc);
  printf("     nsolutions = %d; maxstatefound = %d; \n     current_state = %d, %d, %d, %d;\n     Total Config Exchanges: %d", S_Count, max_state_found, current_state[0], current_state[1], current_state[2],current_state[3], Con_Exe);
  return(0);  
}




