//
// Based in https://github.com/dmitru/mpi-1d-heat-equation/blob/master/solve.c
// and also in https://people.sc.fsu.edu/~jburkardt/cpp_src/heat_mpi/heat_mpi.html
//

#include <iostream>
#include <mpi.h>
#include <cmath>        // M_PI is Pi
#include <cstdlib>

using namespace std;


double initial_temperature(double x){
    return sin(2 * M_PI * x) + 2 * sin(5 * M_PI * x) + 3 * sin(20 * M_PI * x);
}

double compute_update(double Dh, double u_1, double u_2, double u_3) {
    return u_2 + Dh * (u_1 - 2 * u_2 + u_3);
}

int main(int argc, char* argv[]){
    int rank,           // Task identifier
        size,           // Number of tasks
        tag,            // Interaction tag
        n;              //
    
    int M,              // M length intervals
        N,              // N time intervals
        J;              // J grid points for each process
    double T = atof(argv[1]);  // Final time is an input argument
    
    /* Settings for MPI */
    MPI_Comm comm;
    MPI_Status status;
    comm = MPI_COMM_WORLD;
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    if (size < 2){
        cout << "There is only one task.\n";
        cout << "We need at least one worker task!";
        MPI_Abort(comm, 1);
        exit(1);
    }
    
    // Define quantities
    J = 1000;
    M = size * (J - 2) + 2;
    
    // We want dt/dx² < 0.5, thus we can take
    // dt/dx² = T/N * M^2 < 0.5
    //
    N = (int)ceil(3 * T * M * M);
    
    // Discretisation steps
    double dt = T/N;
    double dx = 1.0/M;
    double Dh = dt/(dx*dx);
    
    // 1. Master process fills the buffer with initial data
    double *U = NULL;
    
    if (rank == 0){
        U = (double *) calloc(M, sizeof(double));
        cout << "\ndx = " << dx << ", dt = " << dt << ", dt/dx² = "<< Dh << endl;
        
        //double U[M+1];
        
        for (unsigned int i = 1; i < M; i++){
            U[i] = initial_temperature( i * dx );
        }
        U[0] = 0.0, U[M] = 0.0;
    }
    
    // 2. Processes determine which the subintervals they are responsible for
    int *si_begin = (int *) calloc(size, sizeof(int));
    int *si_end   = (int *) calloc(size, sizeof(int));
    
    
    for (unsigned int ranked = 0; ranked < size; ranked++){
        si_begin[ranked] = ranked * (J-1);
        si_end[ranked] = (ranked+1) * (J-1);
    }
    
    int mine_begin = si_begin[rank];
    int mine_end   = si_end[rank];
    
    
    // 3. Master process distributes the initial values among the processes
    double *buf_1 = (double *) calloc(J, sizeof(double));
    double *buf_2 = (double *) calloc(J, sizeof(double));
    double *my_u = buf_1;
    double *my_u_temp = buf_2;
    
    if (rank == 0){
        // 0 does not have to send to itself
        for (unsigned int ranked = 1; ranked < size; ranked++){
            MPI_Send(U + si_begin[ranked], J, MPI_FLOAT,
                     ranked, 1, comm);
        }
        for (unsigned int i = mine_begin; i < mine_end; i++){
            my_u[i] = U[i];
        }
    } else{
        MPI_Recv(my_u, J, MPI_FLOAT, 0, 1, comm, MPI_STATUS_IGNORE);
    }
    
    
    // 4. Processes start computing the solution
    for (unsigned int i = 0; i < N; i++){
    //for (unsigned int i = 0; i < 10000; i++){
        
        int first_index_to_compute = 1;
        int last_index_to_compute  = J - 2;
        
        for (int j = first_index_to_compute; j <= last_index_to_compute; j++) {
            my_u_temp[j] = compute_update(Dh, my_u[j - 1], my_u[j], my_u[j + 1]);
        }
        
        // Pass the border points to neighbor processes
        // To improve this, we can process one way with even indexes
        // and then process odd indeces
        if (rank % 2 == 0){
            if (rank < size - 1){
                MPI_Send(&my_u_temp[last_index_to_compute], 1, MPI_FLOAT, rank + 1,
                         2, comm);
            }
            if (rank > 0){
                MPI_Send(&my_u_temp[first_index_to_compute], 1, MPI_FLOAT, rank - 1,
                         2, comm);
            }
            if (rank < size - 1) {
                MPI_Recv(&my_u_temp[J - 1], 1, MPI_FLOAT, rank + 1, 2, comm, MPI_STATUS_IGNORE);
            }
            if (rank > 0) {
                MPI_Recv(&my_u_temp[0], 1, MPI_FLOAT, rank - 1, 2, comm, MPI_STATUS_IGNORE);
            }
        } else {
            if (rank > 0){
                MPI_Recv(&my_u_temp[0], 1, MPI_FLOAT, rank - 1,
                         2, comm, MPI_STATUS_IGNORE);
            }
            if (rank < size - 1) {
                MPI_Recv(&my_u_temp[J - 1], 1, MPI_FLOAT, rank + 1,
                         2, comm, MPI_STATUS_IGNORE);
            }
            if (rank > 0){
                MPI_Send(&my_u_temp[first_index_to_compute], 1, MPI_FLOAT, rank - 1,
                         2, comm);
            }
            if (rank < size - 1){
                MPI_Send(&my_u_temp[last_index_to_compute], 1, MPI_FLOAT, rank + 1,
                         2, comm);
            }
        }
        
        // Swap the pointers
        double *temp = my_u;
        my_u = my_u_temp;
        my_u_temp = temp;
    }
    
    // 5. Collect the partial results
    if (rank == 0){
        for (int i = 0; i < mine_end; i++){
            U[i] = my_u[i];
        }
        
        for (int ranked = 1; ranked < size - 1; ranked++){
            MPI_Recv(U + si_begin[ranked] + 1, J - 2,
                     MPI_FLOAT, ranked, 3, comm, MPI_STATUS_IGNORE);
        }
        
        if (size > 1){
            MPI_Recv(U + si_begin[size - 1] + 1, J - 1,
                     MPI_FLOAT, size - 1, 3, comm, MPI_STATUS_IGNORE);
        }
    } else {
        int count = J - 2;
        if (rank == size - 1)
            count = J - 1;
        MPI_Send(&my_u[1], count, MPI_FLOAT, 0, 3, comm);
    }
    
    
    // 6. Print the result out
    if (rank == 0) {
        
        cout << "\nTrue and numerical values at M="<<M<<" space points, N=";
        cout << N << " time points, at time T="<<T<<":"<<endl;
        cout << "\nTrue values           Numerical solutions\n"<<endl;
        double error = 0.0;
        double Usol;
        //int aa = 10000;
        for (unsigned int i = 0; i <= M; i++){
            Usol = exp(-4*M_PI*M_PI*T) * sin(2*M_PI*i*dx) + 2*exp(-25*M_PI*M_PI*T) * sin(5*M_PI*i*dx) + 3*exp(-400*M_PI*M_PI*T) * sin(20*M_PI*i*dx);
            /*Usol = exp(-4*M_PI*M_PI*aa*dt) * sin(2*M_PI*i*dx) + 2*exp(-25*M_PI*M_PI*aa*dt) * sin(5*M_PI*i*dx) + 3*exp(-400*M_PI*M_PI*aa*dt) * sin(20*M_PI*i*dx);*/
            
            error += pow(abs(Usol - U[i]), 2.0);
            
            if (i == 0){
                cout << "Printing first and last 5" << endl;
            }
            if ((i < 5) or (i > M-5)){
                cout << Usol << "            " << U[i] << endl;
            }
        }
        cout << "MAE: " << error/(M+1) << endl;
    }
    
    
    
    
    
    

    // Wrap up
    MPI_Finalize();
    
    return 0;
}


