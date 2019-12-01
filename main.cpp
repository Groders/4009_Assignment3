
# include <cstdlib>
# include "mpi.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;


void print_2darray(bool* grid, int n);
void print_2darray_1(bool* grid, int n, int k);

//Helper function for debugging
void print_2darray(bool* grid, int n){
    for(int z = 0 ; z < n ; z++) {
        for(int y = 0 ; y < n ; y++){
			cout << grid[z*n + y];
        }
        cout << endl;
    }
}

void print_2darray_1(bool* grid, int n, int k){
    for(int z = 0 ; z < n ; z++) {
        for(int y = 0 ; y < n ; y++){
			if (z*n + y >= k){
				return;
			}
			cout << grid[z*n + y];
        }
        cout << endl;
    }
}

//Calculate changes from generation to generation
void calc_changes (bool* grida, bool* gridb, int n, int h, bool* bottem, bool* top){
	for(int z = 0 ; z < h ; z++) {
        for(int y = 0 ; y < n ; y++){
			int count = 0;
			if (h == 1){//edge case
				if(y == 0){
					count = top[y] + top[y+1] + 
								bottem[y] + bottem[y+1];
				}else if(y == n -1){
					count = top[y] + top[y-1] + 
								bottem[y] + bottem[y-1];
				}else {
					count = top[y] + top[y-1] + top[y+1] + 
								bottem[y] + bottem[y-1] + bottem[y+1];
				}
			}else{
				if (h == 0){//is this possible?
					return;
				}else if (y == 0 && z == 0){//Top Left Corner
					count = grida [z*n + y + 1] + top [0] + grida [(z + 1)*n + y] + 
						top [1] + grida [(z +1)*n + y + 1];
				} else if (y == 0 && z == h - 1){//Bottem Left Corner
					count = grida [z*n + y + 1] + grida [(z -1)*n + y] + bottem [y] + 
						grida [(z -1)*n + y + 1] + bottem [y + 1];
				} else if (y == n - 1 && z == 0){//Top Right Corner
					count = grida [z*n + y - 1] + top [y] + grida [(z + 1)*n + y] + 
						top [y -1] + grida [(z + 1)*n + y -1];
				} else if (y == n - 1 && z == h - 1){//Bottem Right Corner
					count = grida [z*n + y - 1] + grida [(z -1)*n + y] + bottem [y] + 
						grida [(z-1)*n + y -1] + bottem [y -1];
				}else if (y == 0) { //Left Edge
					count = grida [z*n + y + 1] + grida [(z -1)*n + y] + grida [(z + 1)*n + y] + 
						grida [(z -1)*n + y + 1] + grida [(z +1)*n + y + 1];
				} else if (y == n -1){ //Right Edge
					count = grida [z*n + y - 1] + grida [(z -1)*n + y] + grida [(z + 1)*n + y] + 
						grida [(z-1)*n + y -1] + grida [(z + 1)*n + y -1];
				}else if (z == 0){//Top Edge
					count = grida [z*n + y - 1] + grida [z*n + y + 1] + top [y] + grida [(z + 1)*n + y] + 
						top [y + 1] + top [y -1] + grida [(z +1)*n + y + 1] + grida [(z + 1)*n + y -1];
				}else if (z == h -1){//Bottem Edge
					count = grida [z*n + y - 1] + grida [z*n + y + 1] + grida [(z -1)*n + y] + bottem [y] + 
						grida [(z -1)*n + y + 1] + grida [(z-1)*n + y -1] + bottem [y + 1] + bottem [y -1];
				}else {//Normal
					count = grida [z*n + y - 1] + grida [z*n + y + 1] + grida [(z -1)*n + y] + grida [(z + 1)*n + y] + 
						grida [(z -1)*n + y + 1] + grida [(z-1)*n + y -1] + grida [(z +1)*n + y + 1] + grida [(z + 1)*n + y -1];
				}
			}
			if (grida[z*n + y] == true && count >= 2 && count <= 3){
				gridb[z*n + y] = true;
			} else if (grida[z*n + y] == false && count == 3) {
				gridb[z*n + y] = true;	
			} else {
				gridb[z*n + y] = false;				
			}
        }
    }
}

int main(int argc, char *argv[])
{
	int id;
	int p;
  	int n, k, m;
	bool *send_data, *gather_data;

	MPI::Init(argc, argv); //  Initialize MPI.
	p = MPI::COMM_WORLD.Get_size(); //  Get the number of processes.
	id = MPI::COMM_WORLD.Get_rank(); //  Get the individual process ID.

	double wtime [1];

	//  Process 0 prints an introductory message.
	if (id == 0) 
	{

		std::string ffn;
		std::string line;

		cout << "\n";
		cout << "--- HELLO_MPI. I am Master Process 0:\n";
		cout << "--- The number of processes is " << p << "\n";
		cout << "\n";
		
		cout << "--- Enter n: \n";
		cin >> n;
		cout << "--- Enter k: \n";
		cin >> k;
		cout << "--- Enter m: \n";
		cin >> m;
		cout << "--- Enter FileName: \n";
		cin >> ffn;
		//n= 5;
		//k = 100;
		//m = 5;
		//ffn = "test-files/test1N100.txt";
		
		std::ifstream  fn(ffn);

		send_data= new bool[n*n];
        for(int i =0; i<n; i++){
			std::getline(fn, line);
			for(int j = 0; j<n;j++){
				send_data [i*n + j] = line[j] == '1' ? true : false;
			}
		}
		
		//debug
		//print_2darray(send_data, n);
		
	} 
	
	int send_buf [3];
	send_buf[0] = n;
	send_buf[1] = k;
	send_buf[2] = m;
	
	MPI_Barrier(MPI_COMM_WORLD);
	wtime[0] = MPI::Wtime();
	MPI_Bcast(&send_buf, 3, MPI_INT, 0, MPI_COMM_WORLD);
	
	//int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)

	int *displs = new int[p + 1];
	int *sendcounts = new int[p];
	int N = send_buf[0];
	int K = send_buf[1];
	int M = send_buf[2];
	
	displs[0] = 0;

	//calculated displacement and sendcounts for each processor
	for (int i = 0; i < p; i++){
		displs[i + 1] = (N*(i+1)/p)*N;
		sendcounts[i] = displs[i + 1] - displs[i];
	}
	
	bool *rvc_buf = new bool[sendcounts[id]];
	
	/*
	for (int i = 0; i < 3; i++) {
		cout << send_buf[i];
		cout << "\n";
    }*/	
	//debug
	/*
	if (id == 0){
		for (int i = 0; i < p; i++)
			cout << sendcounts[i] << " " << displs[i] << " " << "\n";
	}*/
	

	//add the padding to each processors data set (n*2 needs to be added per processor)
	
	MPI_Scatterv(send_data, sendcounts, displs, MPI_C_BOOL, rvc_buf, sendcounts[id], MPI_C_BOOL, 0, MPI_COMM_WORLD);
	
	if (id == 0){
	 free(send_data);
	}


	MPI_Barrier(MPI_COMM_WORLD);


	//send edges to everyone (first row and last row)
	bool *sendbuf = new bool[2*N];
	bool *recvbuf = new bool[2*N];
	int *scounts = new int[p];
	int *sdispls = new int[p];
	int *recvcounts = new int[p];
	int *rdispls = new int[p];
	char* out_file = new char [100];

	bool *grida = rvc_buf;
	bool *gridb = new bool[sendcounts[id]];
	bool *temp;
	if (id == 0){
		gather_data = new bool[n*n];
	}

	//debug
	/*
	for (int i = 0; i < p; i++){
		MPI_Barrier(MPI_COMM_WORLD);
		if (id == i){
			cout << "ID: " << id << "\n\n";
			
			print_2darray_1 (grida, N, sendcounts[id]);
			cout << "\n";
			cout << "\n";

		}
		MPI_Barrier(MPI_COMM_WORLD);
	}*/
	
	if (!scounts || !recvcounts || !rdispls || !sdispls || !gridb || !recvbuf || !recvbuf) {
        cout << "couldn't allocate args\n";
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    
	for (int i = 0; i < K; i++){
		int h = sendcounts[id] / N;
		int t = sendcounts[id];

		for (int j = 0; j < p; j++){
			scounts[j] =0;
			recvcounts[j] = 0;
			sdispls[j] = 0;
			rdispls[j] = 0;
		}
		
		//send buffers
		if (p == 1) {
			scounts[0] = 2*N;
			recvcounts[0] = 2*N;
			for (int l = 0; l < N; l++){
				sendbuf[l] = false;
				sendbuf[l + N] = false;
			}			
		}
		else {			
			if (id == 0){
				scounts[0] = N;
				scounts[1] = N;
				sdispls[1] = N;
				recvcounts[0] = N;
				recvcounts[1] = N;
				rdispls[1] = N;
				for (int l = 0; l < N; l++){
					sendbuf[l] = false;
					sendbuf[l + N] = grida[t - N + l];
				}			
			} else if (id == p - 1){ 
				scounts[p-1] = N;
				scounts[p-2] = N;
				sdispls[p-1] = N;
				recvcounts[p-1] = N;
				recvcounts[p-2] = N;
				rdispls[p-1] = N;
				for (int l = 0; l < N; l++){
					sendbuf[l] = grida [l];
					sendbuf[l + N] = false;
				}
			}else { 
				scounts[id+1] = N;
				scounts[id-1] = N;
				sdispls[id+1] = N;
				recvcounts[id-1] = N;
				recvcounts[id+1] = N;
				rdispls[id+1] = N;
				for (int l = 0; l < N; l++){
					sendbuf[l] = grida[l];
					sendbuf[l + N] = grida[t - N + l];
				}
			}
		}
		
		
		
		//debug
		/*
		for (int l = 1; l < p;  l++){			
			if(l == id){
				cout << id << "\n";
				
				for (int q = 0; q< p; q++) {
					cout << "Send Count " << scounts[q] << " ";
					cout << "Send Displacsment: " <<sdispls[q] << " ";
					cout << "Recive Displacement: " << rdispls[q] << " ";
					cout << "Recive Count: "<< recvcounts[q] << "\n";
				}
				print_2darray_1 (sendbuf, N, 2*N);

				
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		*/
		//send the edges
		MPI_Alltoallv(sendbuf, scounts, sdispls, MPI_C_BOOL, recvbuf, recvcounts, rdispls, MPI_C_BOOL, MPI_COMM_WORLD);
		
		MPI_Barrier(MPI_COMM_WORLD);
		
	
		calc_changes (grida, gridb, N, h, recvbuf + N, recvbuf);
		
		//debug
		/*
		for (int j = 0; j < p; j++){
			MPI_Barrier(MPI_COMM_WORLD);
			if (id == j){
				cout << "ID: " << id << " " << i << "\n\n";
				print_2darray_1 (recvbuf + N, N, N);
				print_2darray_1 (recvbuf, N, N);				
				cout << "\n";
				cout << "\n";
				//print_2darray_1 (grida, N, N*h);
				cout << "\n";
				print_2darray_1 (gridb, N, N*h);
				
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}*/
		//swap a and b 
		temp = grida;
		grida = gridb;
		gridb = temp;
		
		
		//collect and save to file
		MPI_Barrier(MPI_COMM_WORLD);
		if (M != 0 && ((i % M == 0 && i != 0) || i == K - 1)){

			MPI_Gatherv(grida, t, MPI_C_BOOL, gather_data, sendcounts, displs, MPI_C_BOOL, 0, MPI_COMM_WORLD);
			if (id == 0){
				ofstream results_output; 
				sprintf(out_file, "./output_step_%d", i + 1);
				results_output.open(out_file);
				for (int r = 0; r < N*N; r++){
					if (r != 0 && r % n == 0){
						results_output << endl;
					}
					results_output << gather_data[r];
				}
				results_output << endl;
				//print_2darray(gather_data, n);
				//cout << "\n\n";
			}
		}
	}
	double* wtime_arr;
	if (id == 0){
		 wtime_arr = new double[p];
	}
	
	wtime[0] = MPI::Wtime() - wtime[0];

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	
	MPI_Gather(wtime, 1, MPI_DOUBLE, wtime_arr,1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	if (id == 0)
	{
		double longest = 0;
		for (int i = 0; i < p; i ++){
			if (longest < wtime_arr[i]){
				longest = wtime_arr[i];
			}
			//cout << wtime_arr[i] << " ";
		}
		//cout << endl;
		cout << "  Elapsed wall clock time = " << longest << " seconds.\n";
	}
	
	free(displs);
	free(sendcounts);
	free(sendbuf);
	free(recvbuf);
	free(scounts);
	free(sdispls);
	free(recvcounts);
	free(rdispls);
	free(out_file);
	free(gridb);
	free(rvc_buf);

	if (id == 0){
		free(wtime_arr);
		free(gather_data);
	}


//  Terminate MPI.
	MPI::Finalize();
	return 0;
}

