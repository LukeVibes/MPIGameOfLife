
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
//4. improve print functions
//5. ad p1=1 p2=1 functionality to fill_edges




///Matrix Array
int **globalmatrix;
int **mySubMatrix;
int ****submatrices;
int *continousSubMatrices;
int *continousSubMatrix;  
int *temp_contiousSubMatrix;  
int *neighborEdges; 
int *recvNeighborEdges;      
int *finalMatrixSend;
int *recvFinalMatrix;
int *printArr;

///Debug Variables
bool debug = true;

///Global Consts
int OUTOFBOUNDSVALUE = 0;

int bigger(int a, int b){
	int r = -1;
	if (a >= b){
		r = a;
	}
	else{
		r = b;
	}
	return r;
}

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

int sumOfNghbrs(int rank, int p, int x, int y, int n1, int n2, int p1, int p2){

	int i = x-1;
	int j = y-1;
	
	int sum = 0;
	
	int edgelength = bigger(n1, n2);
	int rankdplce = rank * edgelength;
	int crank = 1;
	
	/*
	if(rank ==0){

			cout << "recvNghbrEdges:::  [" ;
			
			
			for(int i=0; i<p*edgelength; i++){
				if(i%5==0){cout << " | ";}
				cout << recvNeighborEdges[i] << ", ";
			}
			cout << "] " << endl<<endl;
	}
	*/


	bool printcond = false;
	
	if(rank == 1 and x==4 and y==0){printcond = true;}

	int index;
	int jndex;
	while(i <= x+1){
		while(j <= y+1){
			///(see Left Edge for comments that apply to all Edges)
			
			///Left Edge
			if(j < 0){
				
				
				
				//left top
				if(i<0){
						if(p1 > 1){
						if(rank-2>=0){
							sum += recvNeighborEdges[(abs(3-rank)*edgelength)+j];
							if (printcond==true){cout << "flag ac1 " << sum << endl;}
						}
						else{
							sum += OUTOFBOUNDSVALUE; 
							if (printcond==true){cout << "flag ac " << sum << endl;}
						}
						}
						else{
							sum += OUTOFBOUNDSVALUE; 
							if (printcond==true){cout << "flag ac " << sum << endl;}
						}
			  }
			  
			  //left under
			  else if(i>=n1){
					  if(p1 > 1){
						if(rank-2>=0){
							sum += OUTOFBOUNDSVALUE; 
							if (printcond==true){cout << "flag ad " << sum << endl;}
						
						}
						else{
							sum += recvNeighborEdges[(abs(1+rank)*edgelength)+j];
							if (printcond==true){cout << "flag ad1 " << sum << endl;}
						}
					}
					else{
							sum += OUTOFBOUNDSVALUE; 
							if (printcond==true){cout << "flag ad " << sum << endl;}
						} 
			  }
				

				
				else{
				if(p2 > 1){
					if(rank%2==0){
						sum+=OUTOFBOUNDSVALUE;
						if (printcond==true){cout << "flag a " << sum << endl;}
						}
					else{
						sum+=recvNeighborEdges[(abs(rank-1)* edgelength)+i];
						if (printcond==true){cout << "flag a1 " << sum << endl;}
						}
				}
				else{
					sum+=OUTOFBOUNDSVALUE;
					if (printcond==true){cout << "flag a " << sum << endl;}
					}	
				}
			}
			
			
			///Right Edge
			else if(j >= n2){
				
				//up-right
				if(i<0){
						if(p1 > 1){
						if(rank-2>=0){
							sum += recvNeighborEdges[(abs(p2-rank)*edgelength)+j];
							if (printcond==true){cout << "flag bc1 " << sum << endl;}
						}
						else{
							sum += OUTOFBOUNDSVALUE; 
							if (printcond==true){cout << "flag bc " << sum << endl;}
						}
						}
						else{
							sum += OUTOFBOUNDSVALUE; 
							if (printcond==true){cout << "flag bc " << sum << endl;}
						}
			  }
			    //under-right
				else if(i>= n1){
						if(p1 > 1){
						if(rank-2>=0){
							sum += OUTOFBOUNDSVALUE; 
							if (printcond==true){cout << "flag bd " << sum << endl;}

						}
						else{
							sum += recvNeighborEdges[(abs(p2+rank)*edgelength)+j];
							if (printcond==true){cout << "flag bd1 " << sum << endl;}
						}
						}
						else{
							sum += OUTOFBOUNDSVALUE; 
							if (printcond==true){cout << "flag bd " << sum << endl;}
						}
					
				}
				
				
				
				else{
				if(p2 > 1){
					if(rank%2==0){
						sum+=recvNeighborEdges[(abs(rank-1)* edgelength)+i];
						if (printcond==true){cout << "flag b1 " << sum << endl;}
						}
					else{
						sum+=OUTOFBOUNDSVALUE;
						if (printcond==true){cout << "flag b " << sum << endl;}
						
						}
				}
				else{
					sum+=OUTOFBOUNDSVALUE;
					if (printcond==true){cout << "flag b " << sum << endl;}
				}
				}
			}
			
			
			///Top Edge
			else if(i < 0){
				if(p1 > 1){
					if(rank-2>=0){
						sum += recvNeighborEdges[(abs(p2-rank)*edgelength)+j];
						if (printcond==true){cout << "flag c1 " << sum << endl;}
					}
					else{
						sum += OUTOFBOUNDSVALUE; 
						if (printcond==true){cout << "flag c " << sum << endl;}
					}
				}
				else{
						sum += OUTOFBOUNDSVALUE; 
						if (printcond==true){cout << "flag c " << sum << endl;}
					}
			}
			
			
			///Bottom Edge
			else if(i >= n1){
				if(p1 > 1){
					if(rank-2>=0){
						sum += OUTOFBOUNDSVALUE; 
						if (printcond==true){cout << "flag d " << sum << endl;}
					
					}
					else{
						sum += recvNeighborEdges[(abs(p2+rank)*edgelength)+j];
						if (printcond==true){cout << "flag d1 " << sum << endl;}
					}
				}
				else{
						sum += OUTOFBOUNDSVALUE; 
						if (printcond==true){cout << "flag d " << sum << endl;}
					}
			}
			
			
			///Normal location
			else {
				if( !(i==x and j==y)){
					
					sum += mySubMatrix[i][j];
						if (printcond==true){cout << "flag e " << sum << endl;}
					
				}
				
			}
			
			j++;
		}
		i++;
		j = y-1;
	}

	return sum;
}

void gameOfLifeRules(int rank, int a, int b, int p1, int p2, int n){
	
	int sum=0;
	for(int i=0; i<a; i++){
		for(int j=0; j<b; j++){
			
			
			sum = sumOfNghbrs(rank, p1*p2, i, j, a, b, p1, p2);
			if(rank == 1){
						cout << "[" << i << "," << j << "]   sum->" << sum <<  endl;
						

					}
			if( mySubMatrix[i][j] == 1){
			
				
				///Lives
				if(sum == 2 or sum == 3){
						
						temp_contiousSubMatrix[(i*(n/p2))+j] = 1; ///should already be 1
				}
				
				///Dies
				if(sum <= 1 or sum >= 4){
						
						temp_contiousSubMatrix[(i*(n/p2))+j] = 0;
				}
			}
			else if (mySubMatrix[i][j] == 0){
					
				if(sum == 3){
						
						temp_contiousSubMatrix[(i*(n/p2))+j] = 1;
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
		free(continousSubMatrices);
		
		
	}
	
	
	///mysubmatrix
	for(int i = 0; i < ((n*n)/(p1*p2))/(n/p1); i++){
		free(mySubMatrix[i]);
	}
	
	free(mySubMatrix);
	free(continousSubMatrix);
	
	free(neighborEdges);
	free(temp_contiousSubMatrix);
	free(recvNeighborEdges);
	
	free(finalMatrixSend);
	free(recvFinalMatrix);
	free(printArr);
}

void fill_neighborEdges(int rank, int p, int n1, int n2, int p1, int p2){
	///n1= (((N*N)/(p1*p2))/(N/p1))----> N/p1
	///n2= (((N*N)/(p1*p2))/(N/p2))----> N/p2
	///p= p
	
	int edgelength = bigger(n1, n2);
	int rankdplce = (rank * edgelength);
	
	///Set blanc values
	
	for(int i=0; i<(p*edgelength); i++){
		neighborEdges[i]= -1;
	}

	///Set actual values
	/// Thanks to the forloops ordering all values will automatically
	/// be added in the correct order!
	for(int i = 0; i< n1; i++){
		for(int j = 0; j< n2; j++){
			if(p2 > 1){
				///Left Edge
				if(j==0){
						if((abs(p2-rank)%2)==0){
							///Has no left neighboor
						}
						else{
							neighborEdges[(abs(rank-1)* edgelength)+i] = mySubMatrix[i][j];
							
							if(i==(n1-1)){
								for(int i=0; i<p; i++){
									neighborEdges[
								}
							}
						}
					
				}
				
				///Right Edge
				if(j==(n2-1)){
					if(rank==0){
						
					}
					if((abs(p2-rank)%2)==0){
		
							neighborEdges[((rank+1)* edgelength)+i] = mySubMatrix[i][j];
						}
						else{
							///Has no right neighboor
						}
				}
			}
			
			if(p1 > 1){
				///Top Edge
				if(i==0){
					if(rank >= ((p1*p2)/2)){
						neighborEdges[(abs(p2-rank)*edgelength)+j] = mySubMatrix[i][j];
					}
					else{
						///Has no top neighboor
					}
					
				}
				
				///Bottom Edge
				if(i==(n1-1)){
					if(rank >= ((p1*p2)/2)){
						///Has no top neighboor
					}
					else{
						neighborEdges[(abs(p2+rank)*edgelength)+j] = mySubMatrix[i][j];
					}
				}	
			}
			
			
			if(p1 == 1 and p2 == 1){
				
			}	
		}
	}
	
	
	
}

void updateSubMatrix(int rank, int n, int p1, int p2){
	int count = 0;
	
	for(int i=0; i<n/p1; i++){
		for(int j=0; j<n/p2; j++){
			
			if(temp_contiousSubMatrix[count] == 0 or temp_contiousSubMatrix[count] == 1){
				mySubMatrix[i][j] = temp_contiousSubMatrix[count];
				
			}
			count++;
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
					cout << submatrices[0][1][i][j];
				}
				cout << endl<<endl;;
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
	int edgelength = bigger(N/p1, N/p2);
  
	temp_contiousSubMatrix = (int *)malloc(sizeof(int *) * ((N/p1) * (N/p2)));  ///traverse each i-> (+N/p2)
  
	///neighborEdges = (int *)malloc(sizeof(int *) * (p1*p2)*((2*(N/p1)) + (2*(N/p2))));
	///recvNeighborEdges = (int *)malloc(sizeof(int *) * (p1*p2)*(((N/p1)) + ((N/p2)))); ///made changes, use to look like one above
	recvNeighborEdges = (int *)malloc(sizeof(int *) * (p1*p2) * edgelength);  
	neighborEdges = (int *)malloc(sizeof(int *) * (p1*p2) * edgelength);  
	finalMatrixSend = (int *)malloc(sizeof(int *) * (N*N));
	recvFinalMatrix = (int *)malloc(sizeof(int *) * (N*N));
	printArr		= (int *)malloc(sizeof(int *) * (N*N));
	
	///recvNeighborBoard = (int *)malloc(sizeof(int *) 
	for (int i=0; i<(p * edgelength); i++){
		recvNeighborEdges[i] = -1;
	}
	
	
	if(rank == 1){
			cout << endl;
			for(int i = 0; i<(N/p1); i++){
				for(int j = 0; j<(N/p2); j++){
					cout << mySubMatrix[i][j] << " ";
				}
				cout << endl;
			}
		
		}
  
 
	for(int w=0; w<3; w++){
		fill_neighborEdges(rank, p, N/p1, N/p2, p1, p2);
		
		/*
		if(rank ==0){

			cout << "sendNghbrEdges:  [" ;
			
			
			for(int i=0; i<p*edgelength; i++){
				if(i%5==0){cout << " | ";}
				cout << neighborEdges[i] << ", ";
			}
			cout << "] " << endl;
		}
		* */
		
		MPI_Alltoall(neighborEdges, edgelength, MPI_INT,
					 recvNeighborEdges, edgelength, MPI_INT, 
					 MPI_COMM_WORLD);
		
		
		if(rank ==0){

			/*
			cout << "recvNghbrEdges:  [" ;
			
			
			for(int i=0; i<p*edgelength; i++){
				if(i%5==0){cout << " | ";}
				cout << recvNeighborEdges[i] << ", ";
			}
			cout << "] " << endl<<endl;
			*/ 
		}
		
		///gameOfLifeRules(rank, ((N*N)/(p1*p2))/(N/p1), ((N*N)/(p1*p2))/(N/p2), p1, p2, N);
		
		gameOfLifeRules(rank, (N/p1), (N/p2), p1, p2, N);
		updateSubMatrix(rank, N, p1, p2);
		
		if(rank == 1){
			cout << endl;
			for(int i = 0; i<(N/p1); i++){
				for(int j = 0; j<(N/p2); j++){
					cout << mySubMatrix[i][j] << " ";
				}
				cout << endl;
			}
			
		}
		
		
		
//Step 5: Send each processors board to Processor-0 to display, every m cycles
/*
		if(((w+1) % m)==0){
			///Prep finalMatrixSend
			for(int i = 0; i<N/p1; i++){
				for(int j=0; j<N/p2; j++){
				
					finalMatrixSend[(i*(N/p2)) + j] = mySubMatrix[i][j]; ///Each Proc puts it in location 0, so only Proc-0 gets data sent to i B^)
				}
			}
			
			///MPI_Gather?
			MPI_Alltoall(finalMatrixSend, ((N/p1) * (N/p2)), MPI_INT,
						 recvFinalMatrix, ((N/p1) * (N/p2)), MPI_INT,
						 MPI_COMM_WORLD);
					 
			if(rank==0){
				///Print RecvFinalMatrix
				///HARDCODED FOR P1*P2=4!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				for(int subV=0; subV<(p1*p2); subV++){
					for(int i = 0; i < N/p1; i++){
						cout << (subV * (N/p2)) + i << " ";
						printArr[(subV * (N/p2)) + i] = recvFinalMatrix[(subV*((N/p1) * (N/p2)))];
					}
					cout << endl;
				}
				cout << endl;
				
				for(int i=0; i<N*N; i++){
					
					cout << printArr[i] << " ";
					if(((i+1) % N) == 0){ cout << endl; }
				}
				
				cout << endl << endl;
			}
			
		}
		* */
		
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



