
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

//Debug Variables
bool debug = true;

void delcareGlobalMatrix(int n) {
	globalmatrix = (char **)malloc(sizeof(char *) * n); 
	for(int i = 0; i < n; i++)
    {
        globalmatrix[i] = (char (*))malloc(sizeof(char) * n);
    }
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

int main(int argc, char *argv[])

{
  int rank;
  int p;
  double wtime;
  
  //MPI set-up
  MPI::Init(argc, argv); 
  p = MPI::COMM_WORLD.Get_size(); 
  rank = MPI::COMM_WORLD.Get_rank(); 
  
  //Assignment Variables
  int p1, p2; // p1, p2: boardgame partitioned into p1 x p2 subrects
  int k;      // k: number of evolutoins
  int N;      // N: boardgame is size NxN
  int m;      // m: at each m-th step proc-0 collects subrects from 
              //    other procs and prints the current configuration
              //    into the output file.
             

  //Step 0: Processor-0 reads in Assigment Variables
	if (rank == 0){
		//HARDCODED
		N = 10;
		p1 = 2;
		p2 = 2;
		k  = 100;
		m  = 25;
		
		delcareGlobalMatrix(N);
	}
  
  //Step 1: Processor-0 reads in NxN binary matrix from input.txt
    if (rank == 0) {
	  string file = "input1.txt";
	  fileToMatrix(file, N);
   }
  
  //Step 2: Each Processor determines his subrectangle
  //HARDCODED
  int startx, endx;
  int starty, endy;
  
  int buf;
  buf = N;
  const root=0;
  MPI_Brcast(&buf, 1, MPI_INT, root, buf);
  
  int nn = buff;
  if (rank==0){
	  startx = 0;
	  endx   = N/2;
	  starty = 0;
	  endy   = N/2;
  }
  if (rank==1) else{
	  startx = nn/2;
	  endx   = nn;
	  starty = 0;
	  endy   = nn/2;
  }
  if (rank==2) else{
	  startx = nn/2;
	  endx   = 0;
	  starty = 0;
	  endy   = nn/2;
  }
  if (rank==3) else{
	  startx = nn/2;
	  endx   = 0;
	  starty = nn/2;
	  endy   = 0;
  }
  
  
  //Assi-Step 2: Proc-0 sends each proc its inital (NxN)/p subrect
  
  
  
  
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

