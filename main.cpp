
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

//1. REMEBER TO TAKE A LOOK AT MYSUBMATRIX--> N/P1  N/P2 RAHTHER THAN MY_N
//2. make 2d arrays int 1d
//3. to declare nieghbor arrays somewhere (neighborEdges, recvNeighborEdges) down there you use p......why? 



///Matrix Array
int **globalmatrix;
int **mySubMatrix;
int ****submatrices;
int *continousSubMatrices;
int *continousSubMatrix;  
int *temp_contiousSubMatrix;  
int *neighborEdges; 
int *recvNeighborEdges;      

///Debug Variables
bool debug = true;

///Global Consts
int OUTOFBOUNDSVALUE = 0;

void declareGlobalMatrix(int n) {
	
	globalmatrix = (int **)malloc(sizeof(int *) * n); 
	for(int i = 0; i < n; i++)
    {
        globalmatrix[i] = (int (*))malloc(sizeof(int) * n);
    }
}

void declareMySubMatrix(int n, int p1, int p2) {
	int my_n = (n*n)/(p1*p2);
	
	mySubMatrix = (int **)malloc(sizeof(int *) * my_n/(n/p1)); 
	for(int i = 0; i < my_n/(n/p1); i++)
    {
        mySubMatrix[i] = (int (*))malloc(sizeof(int) * my_n/(n/p1));
    }
    
    continousSubMatrix = (int *)malloc(sizeof(int *) * my_n);
}

void declareSubMatrices(int n, int p1, int p2){
	
	submatrices = (int ****)malloc(sizeof(int ***) * p1);
	for(int i = 0; i < p1; i++)
	{
		submatrices[i] = (int ***)malloc(sizeof(int **) * p2);
		for(int j = 0; j < p2; j++)
		{
			submatrices[i][j] = (int **)malloc(sizeof(int *) * n);
			for(int k = 0; k < n; k++)
			{
				submatrices[i][j][k] = (int (*))malloc(sizeof(int) * n);
			}
		}
	}
	
	continousSubMatrices = (int *)malloc(sizeof(int *) * (n*n));
}

void fileToMatrix(string file, int n) {
		
		ifstream inputFile;
		inputFile.open(file);

		if (!inputFile) {
		cout << "FILE ERROR: fileToMatrix()" << endl;
		return;
		}
		
		char c;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				inputFile >> c;
				globalmatrix[i][j] = (int) c - '0';
			
				if(debug==true){cout << globalmatrix[i][j] << " ";}
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
	
	/// This For Loop loops poplates eacgh submatrix, using the "submatrices"
	/// 4D array we creates where: 
	///      p1- being the row and 
	///      p2- being the colm
	/// of the grid, of the divided global matrix.
	for(int a = 1; a<= p1; a++){
	for(int i = (((a-1)*n)/a); i < n/(p1 - (a-1)); i++){
		
		for(int b = 1; b<= p2; b++){
		if (debug ==true) {cout << a << "|" << b << "  :";}
		for(int j = (((b-1)*n)/b); j < n/(p2 -(b-1)); j++){
			
			if (debug ==true) {cout << "[" << i << "," << j << "]=[" << subi << "," << subj << "]  ";}
			submatrices[a-1][b-1][subi][subj] = globalmatrix[i][j];
			subj++;
		}
		subj=0;
		cout << endl;}
		subi++;

	}
	subi=0;
	}	
	
	
	/// Now we turn this 4D array into a continous 1D array as MPI does 
	/// not like mulit-dimensional arrays. 
	int index = 0;
	for(int a = 1; a<= p1; a++){
	for(int b = 1; b<= p2; b++){
		
	for(int i = 0; i < ((n*n)/(p1*p2))/(n/p1); i++){
		for(int j = 0; j < ((n*n)/(p1*p2))/(n/p2); j++){

			continousSubMatrices[index] = submatrices[a-1][b-1][i][j];
			index++;
		}
	}

	}}	
	
}

int sumOfNghbrs(int rank, int p, int q, int n1, int n2){

	int i = p-1;
	int j = q-1;
	
	int sum = 0;
	int rankdplce = rank * (n1 + n2);
	
	//WHEN YOU COME BACK CHECK WHAT N1 AND N2 ARE, AND UNDERSTAND HOW P1, P2 MAP TO i, j
	
	
	int index;
	int jndex;
	while(i <= p+1){
		while(j <= q+1){
			///(see Left Edge for comments that apply to all Edges)
			
			///Left Edge
			if(j < 0){
				if(recvNeighborEdges[rankdplce+i] != -1){  ///does this processor have a nighbr at this out of bound (ie has it been sent wall data?)
					sum += recvNeighborEdges[rankdplce+i];
				}
				else{
					sum += OUTOFBOUNDSVALUE; ///if not, it must be at the edge of the actual global graph
				}
			}
			
			
			///Right Edge
			else if(j >= n2){
				if(recvNeighborEdges[rankdplce+n2+i] != -1){ 
					sum += recvNeighborEdges[rankdplce+n2+i];
				}
				else{
					sum += OUTOFBOUNDSVALUE; 
				}
			}
			
			
			///Top Edge
			else if(i < 0){
				if(recvNeighborEdges[rankdplce+n2+n2+j] != -1){ 
					sum += recvNeighborEdges[rankdplce+n2+n2+j];
				}
				else{
					sum += OUTOFBOUNDSVALUE; 
				}
			}
			
			
			///Bottom Edge
			else if(i >= n1){
				if(recvNeighborEdges[rankdplce+n2+n2+n1+j] != -1){ 
					sum += recvNeighborEdges[rankdplce+n2+n2+n1+j];
				}
				else{
					sum += OUTOFBOUNDSVALUE; 
				}
			}
			
			
			///Normal location
			else {
				if( !(i==p and j==q)){
					
					sum += mySubMatrix[i][j];
				}
				
			}
			
			j++;
		}
		i++;
		j = q-1;
	}


}

void gameOfLifeRules(int rank, int a, int b, int p1, int p2, int n){
	
	int sum=0;
	for(int i=0; i<a; i++){
		for(int j=0; j<b; j++){
			
			
			sum = sumOfNghbrs(rank, i, j, a, b);
			if( mySubMatrix[i][j] == 1){
			
				
				///Lives
				if(sum == 2 or sum == 3){
						
						//temp_contiousSubMatrix[(i*(n/p2))+j] = 1; ///should already be 1
				}
				
				///Dies
				if(sum < 1 or sum > 4){
						
						//temp_contiousSubMatrix[(i*(n/p2))+j] = 0;
				}
			}
			else if (mySubMatrix[i][j] == 0){
				if(sum == 3){
						
						//temp_contiousSubMatrix[(i*(n/p2))+j] = 1;
				}
			}
			
		}
	}
	
	//Make actual array temp array.
	//clear temp array?
}

int twoDSubMatrix(int n, int p1, int p2){
	int count = 0;
	
	for(int i=0; i<((n*n)/(p1*p2))/(n/p1); i++){
		for(int j=0; j<((n*n)/(p1*p2))/(n/p2); j++){
			mySubMatrix[i][j] = continousSubMatrix[count];
			count++;
		}
	}
}

void freeAll(int rank, int n, int p1, int p2){
	
	if(rank == 0){
		///globalmatrix
		for(int i = 0; i<n; i++){
			free(globalmatrix[i]);
		}
		free(globalmatrix);
		
		for(int i = 0; i < p1; i++){
			for(int j = 0; j < p2; j++){
				for(int k = 0; k < n; k++){
					free(submatrices[i][j][k]);
				}
				free(submatrices[i][j]);
			}
			free(submatrices[i]);
		}
		free(submatrices);
	
		///continous submatrices
		free(continousSubMatrices);
		
		
	}
	
	
	///mysubmatrix
	for(int i = 0; i < ((n*n)/(p1*p2))/(n/p1); i++){
		free(mySubMatrix[i]);
	}
	free(mySubMatrix);

	
	///continoussubmatrix
	free(continousSubMatrix);
	
	
	///neighboorsWalls
	free(neighborEdges);
}

void fill_neighborEdges(int rank, int p, int n1, int n2){
	int rankdplce = rank * (n1 + n2);
	
	///Set blanc values
	for(int i=0; i<(p*(n1 + n2)); i++){
		neighborEdges[i]= -1;
	}

	///Set actual values
	/// Thanks to the forloops ordering all values will automatically
	/// be added in the correct order!
	for(int i = 0; i< n1; i++){
		for(int j = 0; j< n2; j++){
			///Left Edge
			if(j==0){
				neighborEdges[rankdplce+j] = mySubMatrix[i][j];
			}
			
			///Right Edge
			if(j==(n2-1)){
				neighborEdges[rankdplce+n2+j] = mySubMatrix[i][j];
			}
			
			///Top Edge
			if(i==0){
				neighborEdges[rankdplce+n2+n2+j] = mySubMatrix[i][j];
			}
			
			///Bottom Edge
			if(i==(n1-1)){
				neighborEdges[rankdplce+n2+n2+n1+j] = mySubMatrix[i][j];
			}		
		}
	}
}

int main(int argc, char *argv[])
{
	int rank;
	int p;
	double wtime;

	///MPI set-up
	MPI::Init(argc, argv);
	p = MPI::COMM_WORLD.Get_size(); 
	rank = MPI::COMM_WORLD.Get_rank(); 
	int const ROOT = 0;


	///Assignment Variables
	int p1, p2; /// p1, p2: boardgame partitioned into p1 x p2 subrects
	int k;      /// k: number of evolutoins
	int N;      /// N: boardgame is size NxN
	int m;      /// m: at each m-th step proc-0 collects subrects from 
			  ///    other procs and prints the current configuration
			  ///    into the output file.
 
 
	///Inital Bcast Variables
	int *sendbuf;
	sendbuf = (int (*))malloc(sizeof(int) * 5);
  
//Step 0: Processor-0 reads in Assigment Variables
	if (rank == 0){
		N = 10;
		p1 = 2;
		p2 = 2;
		k  = 100;
		m  = 25;
		
		sendbuf[0] = N;
		sendbuf[1] = p1;
		sendbuf[2] = p2;
		sendbuf[3] = k;
		sendbuf[4] = m;
		
		declareGlobalMatrix(N);
		declareSubMatrices(N, p1, p2);
	}
	
	//---Step 0.1: send Assignment Variables to all processors
	MPI_Bcast(sendbuf, 5, MPI_INT, 0, MPI_COMM_WORLD);
	
	//---Step 0.2: have all processors locally store Assigment Variables (and use them)
	if(rank != 0){
		N  = sendbuf[0];
		p1 = sendbuf[1];
		p2 = sendbuf[2];
		k  = sendbuf[3];
		m  = sendbuf[4];
	}
	free(sendbuf);

	declareMySubMatrix(N, p1, p2);
 
  
//Step 1: Processor-0 reads in NxN binary matrix from input.txt
	if (rank == 0) {
		string file = "input1.txt";
		fileToMatrix(file, N);
	}
   
   
//Step 2: Subdivide Matrix
	if (rank == 0){
		matrixDivider(p1, p2, N);
		
		if(debug==true){
			for(int i=0; i<5; i++){
				for(int j=0; j<5; j++){
					cout << submatrices[0][0][i][j];
				}
				cout << endl;
			}
		}
	}
    
   
 //Step 3: Send out SubMatrices
	if (MPI_Scatter(continousSubMatrices, ((N*N)/(p1*p2)), MPI_INT,
					continousSubMatrix, ((N*N)/(p1*p2)), MPI_INT,
					ROOT, MPI_COMM_WORLD) != MPI_SUCCESS){
					cout << "----------SCATTER ERROR UH-OH---------" << endl;
				}
   
    
 //Step 4: Conduct k evoluntionary steps in SYNC
	twoDSubMatrix(N, p1, p2);
  
	temp_contiousSubMatrix = (int *)malloc(sizeof(int *) * ((N/p1) * (N/p2)));  ///traverse each i-> (+N/p2)
  
	neighborEdges = (int *)malloc(sizeof(int *) * (p1*p2)*((2*(N/p1)) + (2*(N/p2))));
	recvNeighborEdges = (int *)malloc(sizeof(int *) * (p1*p2)*((2*(N/p1)) + (2*(N/p2))));
	
	///recvNeighborBoard = (int *)malloc(sizeof(int *) 
  
	for(int i=0; i<k; i++){
		fill_neighborEdges(rank, p, (((N*N)/(p1*p2))/(N/p1)), (((N*N)/(p1*p2))/(N/p2)));
		
		MPI_Alltoall(neighborEdges, ((2*(N/p1)) + (2*(N/p2))), MPI_INT,
					 recvNeighborEdges, ((2*(N/p1)) + (2*(N/p2))), MPI_INT, 
					 MPI_COMM_WORLD);
		
		//gameOfLifeRules(rank, ((N*N)/(p1*p2))/(N/p1), ((N*N)/(p1*p2))/(N/p2), p1, p2, N);
		gameOfLifeRules(rank, (N/p1), (N/p2), p1, p2, N);
		
//Step 5: Send each processors board to Processor-0 to display, every m cycles

		///MPI_Gather();
	}
	

	
	/// Proc-1 [____, right wall, bottom wall, _______]
	/// Proc-2 [left wall, _____, _______, bottom wall]
	/// Proc-3 [upper wall, _____, ______, right wall ]
	/// Proc-4 [_______, upperwall, left wall,________]
	




	///Wrap-UP
	MPI::Finalize();
	freeAll(rank, N, p1, p2);



	return 0;
}

