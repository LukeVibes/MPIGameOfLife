
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
int *continousSubMatrices;
int *continousSubMatrix;

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
    
    continousSubMatrix = (int *)malloc(sizeof(int *) * my_n);
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
	
	continousSubMatrices = (int *)malloc(sizeof(int *) * n);
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
			continousSubMatrices[index] = (int) submatrices[a-1][b-1][i][j] ;
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
    
   
   
	//Step 3: Send out SubMatrices
	if (MPI_Scatter(continousSubMatrices, (10*10)/(2*2), MPI_INT,
					continousSubMatrix, (10*10)/(2*2), MPI_INT,
					ROOT, MPI_COMM_WORLD) != MPI_SUCCESS){
					cout << "----------SCATTER ERROR UH-OH---------" << endl;
				}
   
   
  //Step 4: Conduct k evoluntionary steps in SYNC
  for(int i=0; i<k; i++){
	  gameOfLifeRUles();
  }
  
  
  /*
   * - prefrom rules on own area
   * 	- use your data
   * 	- but also take into factor the border data
   * - check neighbors borders for changes
   * - send neightbors your borders
   * 
   * 
   * Proc-1 [____, right wall, bottom wall, _______]
   * Proc-2 [left wall, _____, _______, bottom wall]
   * Proc-3 [upper wall, _____, ______, right wall ]
   * Proc-4 [_______, upperwall, left wall,________]
   * 
   */
   
	
  

  
  
  
  //Wrap-UP
  MPI::Finalize();
  free(globalmatrix);
  
  
  
  return 0;
}

void gameOfLifeRules(int **arr, int a, int b){
	
	for(int i=0; i<a; i++){
		for(int j=0; j<b; j++){
			
			if( arr[i][j] == 1){
			
				//Lives
				if(sumOfNeighboors(arr, i, j) == 2 or sumOfNeighboors(arr, i, j) == 3){
						//update temp array
				}
				
				//Dies
				if(sumOfNeighboors(arr, i, j) < 1 or sumOfNeighboors(arr, i, j) > 4){
						//update temp array
				}
			}
			else if (arr[i][j] == 0){
				if(sumOfNeighboors == 3){
						//update temp array
				}
			}
			
		}
	}
	
	//Make actual array temp array.
	//clear temp array?
}

int sumOfNeighboors(int **arr, int a, int b){
	
	
	//Ideas for dealing with edges:
		//- need to include neighboor wall
			//HOW TO ASSOCIATE NEIGHBOR WALL
				//1x1 = none
				//2x1 = two
				//2x2 = eight
				//therefore each submatrix has one bound each row/column
				//2x1 means 2 rows, 1 column, so each sub has 1
				//2x2 means 2 rows, 2 columns so each sub has 2
				//we create an array walls[p1xp2][4], 4 is the max...
				//check wall[submatrix][i] where it has a value
					//-each index represents different location
					//-sooo if you are needing to check values from top wall check that index...
					
		//- if a hard edge, ignore
			//hard edge = no neighbor edge
	
	int sum = 0;
	for(int i= a-1; i<=a+1; a++){
		for(int j=b-1; j<=b+1; b++){
			
			if (!(i=a and j=b)){
				sum += arr[i][j];
			}
		} 
	}
	
	return sum;
}
	

