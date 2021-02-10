#include <iostream>
#include <cmath>
#include <mpi.h>

using namespace std;

static const double PI = 3.1415926536; 

int main(int argc, char* argv[]){

int M = 101;  // M length intervals
int N = 10000; // N time intervals
double T = 0.01; //atof(argv[1]);  // get final time from input argument
int J =27;
double partU[J];
double partU_new[J];// stores the numerical values of function U; two rows to also store values of previous time step 
double U_final[M+1]; 
double Usol[M+1];  // stores true solution 
double dt = T/N;
double dx = 1./M;
double dtdx = dt/(dx*dx);

int rank;
int size = 4;

//initialise MPI

MPI_Init(NULL,NULL);
MPI_Comm_rank( MPI_COMM_WORLD, &rank);
MPI_Comm_size( MPI_COMM_WORLD, &size);

if (rank==0){
cout<< "\ndx="<<dx<<", dt="<<dt<<", dt/dx²="<< dtdx<<endl;
}

// initialize numerical array with given conditions
//start timing
MPI_Barrier(MPI_COMM_WORLD);
double t_start = MPI_Wtime();

for (int j =1; j<=N; j++){
if (rank == 0){
partU[0]=0;
partU_new[0]=0;
if (j==1){
for(int m=1; m<J+1; ++m){
	partU[m] = sin(2*PI*m*dx) + 2*sin(5*PI*m*dx) + 3*sin(20*PI*m*dx);
}
}
for(int m=1; m<J; ++m){
	partU_new[m] = partU[m] + dtdx* (partU[m-1] - 2*partU[m] + partU[m+1]);
}
//set new values to old values
for (int m=1; m<J;++m){
partU[m]=partU_new[m];
}
MPI_Send(&partU[J-1],1,MPI_DOUBLE,rank+1,23,MPI_COMM_WORLD);
MPI_Recv(&partU[J],1,MPI_DOUBLE,rank+1,23,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}
else if (rank == size-1){
partU[-1]=0;
partU_new[-1]=0;
if (j==1){
for(int m=0; m<J; ++m){
	int k=rank*(J-2)-1+m;
	partU[m] = sin(2*PI*k*dx) + 2*sin(5*PI*k*dx) + 3*sin(20*PI*k*dx);
}
}
for(int m=1; m<J; ++m){
        partU_new[m] = partU[m] + dtdx* (partU[m-1] - 2*partU[m] + partU[m+1]);
}
//set new values to old values
for (int m=1; m<J;++m){
partU[m]=partU_new[m];
}
MPI_Send(&partU[1],1,MPI_DOUBLE,rank-1,23,MPI_COMM_WORLD);
MPI_Recv(&partU[0],1,MPI_DOUBLE,rank-1,23,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

else{
if (j==1){
for(int m=0; m<J+1; ++m){
	int k=rank*(J-2)-1+m;
        partU[m] = sin(2*PI*k*dx) + 2*sin(5*PI*k*dx) + 3*sin(20*PI*k*dx);
}
}
for(int m=1; m<J; ++m){
        partU_new[m] = partU[m] + dtdx* (partU[m-1] - 2*partU[m] + partU[m+1]);
}
//set new values to old values
for (int m=1; m<J;++m){
partU[m]=partU_new[m];
}
MPI_Send(&partU[J-1],1,MPI_DOUBLE,rank+1,23,MPI_COMM_WORLD);
MPI_Send(&partU[1],1,MPI_DOUBLE,rank-1,23,MPI_COMM_WORLD);
MPI_Recv(&partU[J],1,MPI_DOUBLE,rank+1,23,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
MPI_Recv(&partU[0],1,MPI_DOUBLE,rank-1,23,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

MPI_Barrier(MPI_COMM_WORLD);
}
for(int m=0; m<J; ++m){
        MPI_Send(&partU[m],1,MPI_DOUBLE,0,m,MPI_COMM_WORLD);
}
//end timing
MPI_Barrier(MPI_COMM_WORLD);
double t_end = MPI_Wtime();


if (rank==0){
for (int i=1; i<J;++i){
U_final[i]=partU[i];
}
U_final[0]=0;
U_final[M]=0;
for (int k =1; k<size; k++){
	for (int j =0; j<J; j++){
		MPI_Recv(&U_final[k*(J-2)+j-1],1,MPI_DOUBLE, k, j ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

}
}
//print out time of computation
cout << "time: " << t_end-t_start << endl;
// print out array entries of numerical solution next to true solution
cout << "\nTrue and numerical values at M="<<M<<" space points at time T="<<T <<":"<<endl;
cout << "\nTrue values           Numerical solutions\n"<<endl;
for(int m=0; m<=M; ++m){
	Usol[m] = exp(-4*PI*PI*T)*sin(2*PI*m*dx) + 2*exp(-25*PI*PI*T)*sin(5*PI*m*dx) + 3*exp(-400*PI*PI*T)*sin(20*PI*M*dx);
	//Usol[m] = exp(-4*PI*PI*0)*sin(2*PI*m*dx) + 2*exp(-25*PI*PI*0)*sin(5*PI*m*dx) + 3*exp(-400*PI*PI*0)*sin(20*PI*M*dx);		
	cout << Usol[m] << "          ;  " << U_final[m] << endl;
	// note that we did not really need to store the true solution in the array just to print out the values.
}
}
MPI_Finalize();
return 0;
}
