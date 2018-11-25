
# include <cstdlib>
# include <iostream>
# include <fstream> 
# include <cilk/cilk.h>
# include <vector>
# include <string>
# include <math.h>  
# include <algorithm>
# include "mpi.h"

using namespace std;

//Matrix Array
char **globalmatrix;
char **mySubMatrix;
char ****submatrices;
char *continousSubMatrices;
char *continousSubMatrix;

//Debug Variables
bool debug = true;

void declareGlobalMatrix(int n) {
	
	globalmatrix = (char **)malloc(sizeof(char *) * n); 
	for(int i = 0; i < n; i++)
    {
        globalmatrix[i] = (char (*))malloc(sizeof(char) * n);
    }
}

void declareMySubMatrix(int n, int p1, int p2) {
	int my_n = (n*n)/(p1*p2);
	
	mySubMatrix = (char **)malloc(sizeof(char *) * my_n); 
	for(int i = 0; i < my_n; i++)
    {
        mySubMatrix[i] = (char (*))malloc(sizeof(char) * my_n);
    }
    
    continousSubMatrix = (char *)malloc(sizeof(char *) * my_n);
}

void declareSubMatrices(int n, int p1, int p2){
	
	submatrices = (char ****)malloc(sizeof(char ***) * p1);
	for(int i = 0; i < p1; i++)
	{
		submatrices[i] = (char ***)malloc(sizeof(char **) * p2);
		for(int j = 0; j < p2; j++)
		{
			submatrices[i][j] = (char **)malloc(sizeof(char *) * n);
			for(int k = 0; k < n; k++)
			{
				submatrices[i][j][k] = (char (*))malloc(sizeof(char) * n);
			}
		}
	}
	
	continousSubMatrices = (char *)malloc(sizeof(char *) * n);
}

void fileToMatrix(string file, int n) {
		
		ifstream inputFile;
		inputFile.open(file);

		if (!inputFile) {
		cout << "FILE ERROR: fileToMatrix()" << endl;
		return;
		}
		
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				inputFile >> globalmatrix[i][j];
				if(debug==true){cout << globalmatrix[i][j];}
			}
			if(debug==true){cout << endl;}
		}

			
		inputFile.close();
}

void printGMatrix(int n){
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cout << globalmatrix[i][j];
			}
			cout << endl;
		}

}




void matrixDivider(int p1, int p2, int n){
	int subi = 0;
	int subj = 0;
	
	cout << p1 << endl;
	
	for(int a = 1; a<= p1; a++){
	for(int i = (((a-1)*n)/a); i < n/(p1 - (a-1)); i++){
		
		for(int b = 1; b<= p2; b++){
		cout << a << "|" << b << "  :";
		for(int j = (((b-1)*n)/b); j < n/(p2 -(b-1)); j++){
			
			cout << "[" << i << "," << j << "]=[" << subi << "," << subj << "]  ";
			submatrices[a-1][b-1][subi][subj] = globalmatrix[i][j];
			subj++;
		}
		subj=0;
		cout << endl;}
		subi++;

	}
	subi=0;
	}	
	
	//Make continous NOT SURE IF THIS WORKS WITH 2*1 :(
	int index = 0;
	for(int a = 1; a<= p1; a++){
	for(int b = 1; b<= p2; b++){
		
	for(int i = 0; i < ((n*n)/(p1*p2))/(n/p1); i++){
		for(int j = 0; j < ((n*n)/(p1*p2))/(n/p2); j++){
			
			cout << submatrices[a-1][b-1][i][j];
			continousSubMatrices[index] = submatrices[a-1][b-1][i][j] ;
			index++;
		}
		cout << endl;
	}

	}}	
}



int main(int argc, char *argv[])
{
  int rank;
  int p;
  double wtime;
  
  //MPI set-up
  MPI::Init(argc, argv);
  p = MPI::COMM_WORLD.Get_size(); 
  rank = MPI::COMM_WORLD.Get_rank(); 
  int *buf;
  buf = (int *) malloc(p*sizeof(int));
  int const ROOT = 0;
  
  //Assignment Variables
  int p1, p2; // p1, p2: boardgame partitioned into p1 x p2 subrects
  int k;      // k: number of evolutoins
  int N;      // N: boardgame is size NxN
  int m;      // m: at each m-th step proc-0 collects subrects from 
              //    other procs and prints the current configuration
              //    into the output file.
  
  
  char rbuf[10][10];
  
  //Step 0: Processor-0 reads in Assigment Variables
	if (rank == 0){
		//HARDCODED
		N = 10;
		p1 = 2;
		p2 = 2;
		k  = 100;
		m  = 25;
		
		declareGlobalMatrix(N);
		declareSubMatrices(N, p1, p2);
	}
	
	//NEED TO BROADCAST N AND P1, P2
	declareMySubMatrix(10, 2, 2);
  
  //Step 1: Processor-0 reads in NxN binary matrix from input.txt
	if (rank == 0) {
		string file = "input1.txt";
		fileToMatrix(file, N);
	}
   
   //Step 2: Subdivide Matrix
	if (rank == 0){
		matrixDivider(p1, p2, N);
		
		for(int i=0; i<=5; i++){
			for(int j=0; j<5; j++){
				cout << submatrices[0][0][i][j];
			}
			cout << endl;
		}
	}
    
   
   
   
	if (MPI_Scatter(continousSubMatrices, (10*10)/(2*2), MPI_CHAR,
					continousSubMatrix, (10*10)/(2*2), MPI_CHAR,
					ROOT, MPI_COMM_WORLD) != MPI_SUCCESS){
					cout << "----------SCATTER ERROR UH-OH---------" << endl;
				}
   
    
    if (rank == 3){
		for(int i = 0; i<25; i++){
			cout << continousSubMatrix[i];
			if ((((i+1) % 5) == 0) && (i != 0)){
				cout << endl;
			}
		}
	}
	cout << endl;
   
   
   
   
   
   
  
  
  /*
  //Step 2: Each Processor determines his subrectangle
  //HARDCODED
  int startx, endx;
  int starty, endy;
  
  
  if (rank == 0){
	 buf[0] = N;
  }
  MPI_Bcast(&buf, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
 
  
  int nn = buf[0];
  
  if (rank==0){
	  startx = 0;
	  endx   = N/2;
	  starty = 0;
	  endy   = N/2;
  }
  else if (rank==1){
	  startx = nn/2;
	  endx   = nn;
	  starty = 0;
	  endy   = nn/2;
  }
  else if (rank==2){
	  startx = nn/2;
	  endx   = 0;
	  starty = 0;
	  endy   = nn/2;
  }
  else if (rank==3){
	  startx = nn/2;
	  endx   = 0;
	  starty = nn/2;
	  endy   = 0;
  }
  
  //Testing
  /*
  char c = 2;
  if (rank==3){
	  for(int i=starty; i<endy; i++){
		  for(int j=startx; i<endx; i++){
			  
			  globalmatrix[i][j] = c;
		  }
	  }
	  cout << endl;
	  printMatrix(nn);
  }
  */
  
  
  
  
  
  
  //Assi-Step 3: ALL procs execute k evolutionary steps of the game in 
  //             SYNCHRONOUS fashion
  
		//Assi-Step 4.1: At each m-th step, Proc-0 collects subrects from other 
		//               procs and prints whole game board to output.txt WITH
		//               runntime

  
  
  
  //Wrap-UP
  MPI::Finalize();
  free(globalmatrix);
  
  
  
  return 0;
}

