
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


//. Have all values send to Proc-0                 14:30 DONE
//. clean up code										 DONE
//. Have Proc-0 print to file final array          16:15 DONE
//. Add runttime timer							   16:45 DONE
//. GET PROPER EDGE SENDING!                       22:30
		//. Design								19:10
		//. Test								22:10
		//. Final Code							23:10

//. add p1=1 p2=1 functionality to fill_edges
//. test with p1 != p2 values
//. test with N=odd values
//. get rid of nD arrays



///Matrix Array
int **globalmatrix;
int **mySubMatrix;
int ****submatrices;
int *continousSubMatrices;
int *continousSubMatrix;  
int *temp_contiousSubMatrix;  
int *neighborEdges; 
int *Edges;
int *recvEdges;
int *recvNeighborEdges;      
int *sendFinal;
int *recvFinal;
int *printArr;

///Debug Variables
bool debug = true;

///Global Consts
int OUTOFBOUNDSVALUE = 0;

///Timmer Values
double wtime;

///bigger: used for simple max calculations
int bigger(int a, int b){
	int r = a;
	if (a >= b){
		r = a;
	}
	else{
		r = b;
	}
	return r;
}

///sAdd: when doing addition on rank values, if sum is greater than greatest
///      rank value, returns intial rank value. (Needed for scalable algebra in
///       edge detection and edge surrounding cell addition.
int sAdd(int rank, int v, int maxrank){
	if( (rank + v) > maxrank ){ return rank;}
	else{					   return rank+v;}
}

///sSub: similar to sAdd but if the rank - v is less than zero, returns the
///      intial rank.
int sSub(int rank, int v){
	if( (rank-v) < 0 ){ return rank;}
	else{				return rank-v;}
}

int special_i(int i, int N, int rank, int p1){
	int sp_i;
	if(rank < p1){
		sp_i = ((N/p1) -1);
	}
	else if(rank >= p1){
		sp_i = 0;
	}
	
	return sp_i;
}

int special_j(int j, int N, int rank, int p2){
	int sp_j;

	if(rank%2==0){
		sp_j = ((N/p2) -1);
	}
	else{
		sp_j = 0;
	}

	return sp_j;
}

///declareGlobalMatrix: simple malloc on gobalmatrix dataset
void declareGlobalMatrix(int n) {
	
	globalmatrix = (int **)malloc(sizeof(int *) * n); 
	for(int i = 0; i < n; i++)
    {
        globalmatrix[i] = (int (*))malloc(sizeof(int) * n);
    }
}

///declareMySubMatrix: simple malloc on mySubMatrix and continousSubMatrix
void declareMySubMatrix(int n, int p1, int p2) {
	int my_n = (n*n)/(p1*p2);
	
	mySubMatrix = (int **)malloc(sizeof(int *) * my_n/(n/p1)); 
	for(int i = 0; i < my_n/(n/p1); i++)
    {
        mySubMatrix[i] = (int (*))malloc(sizeof(int) * my_n/(n/p1));
    }
    
    continousSubMatrix = (int *)malloc(sizeof(int *) * my_n);
}


///declareSubMatrices; simple malloc on submatrices and continoussubmatrices
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


///fileToMatrix: turns input.txt into gobalmatrix
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


///printGMatrix: simply prints gobalmatrix
void printGMatrix(int n){
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cout << globalmatrix[i][j];
			}
			cout << endl;
		}

}


///matrixDivider: turns the 2D globalmatrix into the 4D submatrices and then,
///               turns submtrices into the 2D contionussubmatrices 
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


///printFinaltoFile: simply prints recvFinal to output.txt, along with runntime
void printFinaltoFile(int N, int p1, int p2, double time, int p){
				int rowlength = N/p2;
				int submtrxSize = (N/p1) * (N/p2);  
			
				///File Streaming Set-up
				fstream ofs;
				ofs.open("output.txt", ios::out | ios::trunc);
				
				
				for(int i=0; i<N; i++){
					for(int j=0; j<N; j++){
						if(i<N/p1){
							if(j<N/p2){
								ofs << recvFinal[ (rowlength*i) + j ] << " ";
							}
							
							else if(j>=N/p2){
								ofs << recvFinal[ (rowlength*i) + submtrxSize + (j-(N/p2)) ] << " ";
							}
						}
						
						
						else if(i>=N/p1){
							if(j<N/p2){
								ofs << recvFinal[ (rowlength*(i-N/p1)) + (p1*submtrxSize) + j ] << " ";
							}
							
							else if(j>=N/p2){
								ofs << recvFinal[ (rowlength*(i-N/p1)) + (p1*submtrxSize) + submtrxSize + (j-(N/p2)) ] << " ";
							}
						}
					
					}
					ofs << endl;
				}
				ofs << endl;
				
				ofs << "Time            : " << time << endl;
				ofs << "Time/#processors: " << time/p << endl;
				
				///File Stream Wrap-up
				ofs.close();
				
				
			

}


///sumOfNghbrs: given and (x,y) from mySubMatrix, gives the sum of the surrounding
///             8 values, utiliing the neighbors edge values aswell
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
	
	if(rank == 5 and x==44 and y==44){printcond = true;}

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



int sumOfNeighbors(int x, int y, int rank, int N, int p1, int p2){
	//NOT OPTIMIZES FOR P1xP2==2
	
	int sum = 0;
	int max = bigger(N/p1, N/p2);
	
	bool testable = false;
	
	//if(rank==2 and x==1 and y==4){testable=true;}
		

					
	
	if(p1*p2==4){
				if(x==special_i(x, N, rank, p1) and y==special_j(y, N, rank, p2)){
					sum += recvEdges[ ((3-rank)*max)  +  max-1 ];
					if(testable){cout << "corner hit  sum: "<< sum << endl;}
					//continue; //TEST
					
				}
			}
	
	for(int i=x-1; i<=x+1; i++){
		for(int j=y-1; j<=y+1; j++){
			
			if(!(i==x and j==y))
			{
			
			///Left Case
			if(j<0){
				if(rank%2==0){sum += 0;}
				else{sum += recvEdges[ (sSub(rank, 1) * max) + i ];}
				if(testable){cout << "left hit  location: "<< (sSub(rank, 1) * max) + i << "  sum: " << sum << endl;}
			}
			
			///Right Case
			else if(j>((N/p2)-1)){
				if(!(rank%2==0)){sum += 0;}
				else{sum += recvEdges[ (sAdd(rank, 1,  (p1*p2)-1) * max) + i ];}
				if(testable){cout << "right hit  sum: "<< sum <<  "  location: " << (sAdd(rank, 2,  (p1*p2)-1) * max) + i << " sAdd:" <<sAdd(rank, 2,  (p1*p2)-1)  << endl;}
			}
			
			///Top Case
			else if(i<0){
				if(rank>=((p1*p2)/2)){sum += recvEdges[ (sSub(rank, 2) * max) + j ];}
				else{sum += 0;}
				if(testable){cout << "top hit  sum: "<< sum << "  location: " <<  (sSub(rank, 2) * max) + j << endl;}
			}
			
			///Bottom Case
			else if(i>((N/p1)-1)){
				if(rank<((p1*p2)/2)){sum += recvEdges[ (sAdd(rank, 2, (p1*p2)-1) * max) + j ];}
				else{ sum += 0;}
				if(testable){cout << j << " bottom hit  sum: "<< sum << "  location: " <<  (sAdd(rank, 2, max-1) * max) + j << "  value " <<  recvEdges[ (sAdd(rank, 2, max-1) * max) + j ] << endl;}
			}
			
			///Regular Case
			else{
				sum+= mySubMatrix[i][j];
				if(testable){cout << "regular hit  sum: "<< sum << endl;}
			}
		
			}	
		}
	} 
	
	
	if(testable){cout << "FINAL SUM: " << sum << endl;}
	return sum;
}


///gameOfLifeRules: using the game of life rules decided whether to kill or birth or
///                 keep alive values in mySubMatrix
void gameOfLifeRules(int p1, int p2, int n, int rank){
	
	int sum=0;
	int max = bigger(n/p1, n/p2);
	
	for(int i=0; i<n/p1; i++){
		for(int j=0; j<n/p2; j++){
			
			
			sum = sumOfNeighbors(i, j, rank, n, p1, p2);
			if(rank == 2){
						//cout << "[" << i << "," << j << "]   sum->" << sum <<  endl;
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


///twoDSubMatrix: simply turns the 1D contiousSubMatrix to the 2D mySubMatrix
int twoDSubMatrix(int n, int p1, int p2){
	int count = 0;
	
	for(int i=0; i<((n*n)/(p1*p2))/(n/p1); i++){
		for(int j=0; j<((n*n)/(p1*p2))/(n/p2); j++){
			mySubMatrix[i][j] = continousSubMatrix[count];
			count++;
		}
	}
}

///freeAll; simmply frees all malloc-ed datasets
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
	
	free(Edges);
	free(recvEdges);
	free(temp_contiousSubMatrix);
	//free(recvNeighborEdges);
	
	free(sendFinal);
	free(recvFinal);
	free(printArr);
}


///fill_neighborEdges: populates neighborEdges with the appropraite edge values
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

///updateSubMatrix: simply updates the 2D mySubMatrix with the value from the
///                 1D temp_continousSubMatrix
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



void populateMyEdges(int p1, int p2, int N, int rank){
	
	bool testable=false;
	
	//if(rank==3){testable=true;}
	
	if(testable){cout << "function hit" << endl;}
	
	
	
	int max = bigger(N/p1, N/p2);
	for(int i=0; i<N/p1; i++){
		for(int j=0; j<N/p2; j++){
			
				if(p1*p2==2){
					///Only left/right edges needed
					if(p1>1){
						if(i==N/p1 or i==0){
							Edges[(abs(rank-1)*N) + j] = mySubMatrix[i][j];
						}
					}

					///Only bottom/top edges needed
					else if(p2>1){
						if(i==N/p2 or i==0){
							Edges[(abs(rank-1)*N) + i] = mySubMatrix[i][j];
						}
					}
				}


				if(p1*p2>2){
					///Corner Cases
					if(i==special_i(i, N, rank, p1) and j==special_j(j, N, rank, p2)){
						Edges[ ((3-rank) * max) + (max - 1)] = mySubMatrix[i][j];
						if(testable){cout<<"special hit [" << i << "," << j << "]   value: " <<  mySubMatrix[i][j] << "  location : " <<  ((3-rank) * max) + (max - 1) << endl;}
					}
					
					///Left Edge
					if(j == 0){
						if(!(rank%2==0)){
							Edges[ ((rank-1)*max) + i ] = mySubMatrix[i][j];
							if(testable){cout<<"left hit [" << i << "," << j << "] " << endl;}
						}
					}
					
					///Right Edge
					if(j == ((N/p2) -1)){
						if(rank%2==0){
							Edges[ ((rank+1)*max) + i ] = mySubMatrix[i][j];
							if(testable){cout<<"right hit [" << i << "," << j << "]  location: " <<(rank+1) + i  << endl;}
						}
					}

					///Top Edge
					if(i == 0){
						if(rank>=(p1*p2/2)){
							Edges[ ((rank-2) * max) + j ] = mySubMatrix[i][j];
							if(testable){cout<<"top hit [" << i << "," << j << "]  value: "<<  mySubMatrix[i][j] << endl;}
						}
					}
					
					///Bottom Edge
					if(i == ((N/p2) -1)){
						if(rank<(p1*p2/2)){
							Edges[ ((rank+2) * max) + j ] = mySubMatrix[i][j];
							if(testable){cout<<"bottom hit [" << i << "," << j << "]  location: " << ((rank+2) * max) + j << endl;}
						}
					}

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
	
	//~~~Start Timmer~~~
	if (rank == 0) 
	{
		wtime = MPI::Wtime();
	}
  
//Step 0: Processor-0 reads in Assigment Variables
	if (rank == 0){
		N = 20;
		p1 = 2;
		p2 = 2;
		k  = 100;   ////Should be: 100
		m  = 25;   ////Should be: 25
		
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
		string file = "input3.txt";
		fileToMatrix(file, N);
	}
   
   
//Step 2: Subdivide Matrix
	if (rank == 0){
		matrixDivider(p1, p2, N);
	}
	
 //Step 3: Send out SubMatrices
	if (MPI_Scatter(continousSubMatrices, ((N*N)/(p1*p2)), MPI_INT,
					continousSubMatrix, ((N*N)/(p1*p2)), MPI_INT,
					ROOT, MPI_COMM_WORLD) != MPI_SUCCESS){
					cout << "----------SCATTER ERROR UH-OH---------" << endl;
	}
   
    
 //Step 4: Conduct k evoluntionary steps in SYNC
	///Bring continous matrix back to 2D Form (for convience)
	twoDSubMatrix(N, p1, p2);
	int max = bigger(N/p1, N/p2);
	
	///Datasets needed to communicate Neighbors Edges to each other
	temp_contiousSubMatrix = (int *)malloc(sizeof(int *) * ((N/p1) * (N/p2)));  ///traverse each i-> (+N/p2)
	recvEdges = (int *)malloc(sizeof(int *) * (p1*p2) * max);  
	Edges     = (int *)malloc(sizeof(int *) * (p1*p2) * max);  
	
	///Datasets needed for printing final matrix to output.txt
	sendFinal		= (int *)malloc(sizeof(int *) * (N*N));
	recvFinal 		= (int *)malloc(sizeof(int *) * (N*N));

	///Initalize Edge Dataset to avoid errors (and for testing convience)
	for (int i=0; i<(p * max); i++){
		Edges[i] = 0;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
			cout << "rank " << rank << " SubMatrix" << endl;
			for(int i=0; i<N/p1; i++){
				for(int j=0; j<N/p2; j++){
					cout << mySubMatrix[i][j] << " ";
				}
				cout << endl;
			}
		}
	
	///Evolve Game Of Life board!
	for(int w=0; w<k; w++){
		
		///Collect processors submatrix edge values
		populateMyEdges(p1, p2, N, rank);
	
		//Edge print
		if(1==0){
					MPI_Barrier(MPI_COMM_WORLD);
					if(rank == 0){
						cout<< "rank " << rank <<  " Edges[]: "<< endl;
						for(int i=0; i<(p1*p2) * max; i++){
							if(i%5==0){cout<<" | ";}
							cout << Edges[i] << " ";
							
						}
						cout << endl<< endl;
					}
					
					MPI_Barrier(MPI_COMM_WORLD);
					if(rank == 1){
						cout<< "rank " << rank <<  " Edges[]: "<< endl;
						for(int i=0; i<(p1*p2) * max; i++){
							if(i%5==0){cout<<" | ";}
							cout << Edges[i] << " ";
						}
						cout << endl<< endl;
					}
					
					MPI_Barrier(MPI_COMM_WORLD);
					if(rank == 2){
						cout<< "rank " << rank <<  " Edges[]: "<< endl;
						for(int i=0; i<(p1*p2) * max; i++){
							if(i%5==0){cout<<" | ";}
							cout << Edges[i] << " ";
						}
						cout << endl<< endl;
					}
					
					MPI_Barrier(MPI_COMM_WORLD);
					if(rank == 3){
						cout<< "rank " << rank <<  " Edges[]: "<< endl;
						for(int i=0; i<(p1*p2) * max; i++){
							if(i%5==0){cout<<" | ";}
							cout << Edges[i] << " ";
						}
						cout << endl<< endl;
					}
		}
		
		
		///Begin with communicating processors with neighbors edge values
		MPI_Alltoall(Edges, max, MPI_INT,
					 recvEdges, max, MPI_INT, 
					 MPI_COMM_WORLD);
		
		
		///Conduct Game Of Life Rules onto processors submatrix
		gameOfLifeRules(p1, p2, N, rank);
		updateSubMatrix(rank, N, p1, p2);
		
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank==0){
			cout << "rank " << rank << " SubMatrix" << endl;
			for(int i=0; i<N/p1; i++){
				for(int j=0; j<N/p2; j++){
					cout << mySubMatrix[i][j] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}
		
		
//Step 5: Send each processors board to Processor-0 to display, every m cyclesS
		
		if(((w+1) % m)==0){
			///Populate sendFinal with processors subMatrix values
			int rowlength = N/p2;
			for(int i = 0; i<N/p1; i++){
				for(int j=0; j<N/p2; j++){
				
					/// Each processor will save data in the "Proc-0 section", that way
					/// Proc-0 gets is the only one recvieving values (transpose logic).
					sendFinal[(i*rowlength) + j] = mySubMatrix[i][j];
				}
			}
			
			///Send values to Proc-0
			MPI_Alltoall(sendFinal, ((N/p1) * (N/p2)), MPI_INT,
						 recvFinal, ((N/p1) * (N/p2)), MPI_INT,
						 MPI_COMM_WORLD);
			

			//~~~End Timmer~~~
			if (rank == 0)
			{
				wtime = MPI::Wtime() - wtime;
			}

			///Print final matrix to output.txt
			if(rank==0){
				printFinaltoFile(N, p1, p2, wtime, p);	
			}
			
			////SubMatrix by SubMatrix print for testing
			if(1==0){
				MPI_Barrier(MPI_COMM_WORLD);
				if(rank == 0){
					cout << "Rank 0 mySubMatrix" << endl;
					for(int i = 0; i<(N/p1); i++){
						for(int j = 0; j<(N/p2); j++){
							cout << sendFinal[(i*rowlength) + j] << " ";
						}
						cout << endl;
					}
					
				}
				
				MPI_Barrier(MPI_COMM_WORLD);
				if(rank == 1){
					cout << "Rank 1 mySubMatrix" << endl;
					for(int i = 0; i<(N/p1); i++){
						for(int j = 0; j<(N/p2); j++){
							cout << sendFinal[(i*rowlength) + j] << " ";
						}
						cout << endl;
					}
					
				}
				
				MPI_Barrier(MPI_COMM_WORLD);
				if(rank == 2){
					cout << "Rank 2 mySubMatrix" << endl;
					for(int i = 0; i<(N/p1); i++){
						for(int j = 0; j<(N/p2); j++){
							cout << sendFinal[(i*rowlength) + j] << " ";
						}
						cout << endl;
					}
				}
				
				MPI_Barrier(MPI_COMM_WORLD);
				if(rank == 3){
					cout << "Rank 3 mySubMatrix" << endl;
					for(int i = 0; i<(N/p1); i++){
						for(int j = 0; j<(N/p2); j++){
							cout << sendFinal[(i*rowlength) + j] << " ";
						}
						cout << endl;
					}
				}
				
			}
			
		}
		
	}
	
	///Wrap-UP
	MPI::Finalize();
	freeAll(rank, N, p1, p2);



	return 0;
}



