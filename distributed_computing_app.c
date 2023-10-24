#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define NUM_CLUSTERS 4
#define ROOT 0

int rank, nProcesses, arraySize, commError, nWorkers, neigh[2];
int *worker, *parents, *workStart;

void read_workers() {
    if (rank >= NUM_CLUSTERS)
        return;
    
    char fileName[30];
    memset(fileName, 0, sizeof(fileName));
    switch (rank) 
    {
    case 0:
        strcpy(fileName, "cluster0.txt");
        break;
    case 1:
        strcpy(fileName, "cluster1.txt");
        break;
    case 2:
        strcpy(fileName, "cluster2.txt");
        break;
    case 3:
        strcpy(fileName, "cluster3.txt");
        break;
    default:
        fprintf(stderr, "Invalid cluster to read: %d.\n", rank);
        break;
    }

    FILE *in = fopen(fileName, "r");
    fscanf(in, "%d", &nWorkers);

    worker = (int *) calloc(nWorkers, sizeof(int));
    if (worker == NULL)
        fprintf(stderr, "Could not allocate vector of %d workers for cluster %d.\n", nWorkers, rank);

    for (int i = 0; i < nWorkers; i++)
        fscanf(in, "%d", &worker[i]);
}


void create_cluster_neighs() {
    if (rank >= NUM_CLUSTERS)
        return;
    
    if (commError == 0) {
        switch (rank)
        {
        case 0:
            neigh[0] = 3;
            neigh[1] = 1;
            break;
        case 1:
            neigh[0] = 0;
            neigh[1] = 2;
            break;
        case 2:
            neigh[0] = 1;
            neigh[1] = 3;
            break;
        case 3:
            neigh[0] = 2;
            neigh[1] = 0;
            break;
        default:
            fprintf(stderr, "Invalid cluster in topology(commError = %d): %d.\n", commError, rank);
            break;
        }
    } else if (commError == 1) {
        switch (rank)
        {
        case 0:
            neigh[0] = 3;
            neigh[1] = -1;
            break;
        case 1:
            neigh[0] = -1;
            neigh[1] = 2;
            break;
        case 2:
            neigh[0] = 1;
            neigh[1] = 3;
            break;
        case 3:
            neigh[0] = 2;
            neigh[1] = 0;
            break;
        default:
            fprintf(stderr, "Invalid cluster in topology(commError = %d): %d.\n", commError, rank);
            break;
        }
    } else if (commError == 2) {
        switch (rank)
        {
        case 0:
            neigh[0] = 3;
            neigh[1] = -1;
            break;
        case 1:
            neigh[0] = -1;
            neigh[1] = -1;
            break;
        case 2:
            neigh[0] = -1;
            neigh[1] = 3;
            break;
        case 3:
            neigh[0] = 2;
            neigh[1] = 0;
            break;
        default:
            fprintf(stderr, "Invalid cluster in topology(commError = %d): %d.\n", commError, rank);
            break;
        }
    } else
        fprintf(stderr, "Invalid commError( %d ) for cluster %d.\n", commError, rank);
}

void init_topology() {
    // Vector of fathers for all nodes in topology
    parents = (int *) calloc(nProcesses, sizeof(nProcesses));
    if (parents == NULL)
        fprintf(stderr, "Could not allocate the parents vector for cluster %d.\n", rank);
    if (rank < NUM_CLUSTERS) {
        // Set all to -1 as if every worker is a cluster
        // Only the main 4 cluster will be interconacted
        for (int i = 0; i < nProcesses; i++)
            parents[i] = -1;

        for (int i = 0; i < nWorkers; i++)
            parents[ worker[i] ] = rank;
    }
}

void sync_topology() {
    if (rank < NUM_CLUSTERS) {
        int send_vec[nProcesses], recv_vec[nProcesses];
        if (commError == 0) {
            // inel topology so we use the inel algorithm to communicate
            // we need to run a full comunication 2 times in order to all clusters know complete topology
            // ROOT(0) will initiate the comunication
            for (int t = 0; t < 2; t++) {
                if (rank == ROOT) {
                    for (int j = 0; j < nProcesses; j++)
                        send_vec[j] = parents[j];
                    MPI_Send(send_vec, nProcesses, MPI_INT, neigh[1], 0, MPI_COMM_WORLD);
                    fprintf(stdout, "M(%d,%d)\n", rank, neigh[1]);

                    MPI_Status status;
                    MPI_Recv(recv_vec, nProcesses, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);
                    for (int i = 0; i < nProcesses; i++) {
                        if (parents[i] == -1 && recv_vec[i] != -1)
                            parents[i] = recv_vec[i];
                    }
                } else {
                    // receive info and update my topology parent list
                    MPI_Status status;
                    MPI_Recv(recv_vec, nProcesses, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);
                    for (int i = 0; i < nProcesses; i++) {
                        if (parents[i] == -1 && recv_vec[i] != -1)
                            parents[i] = recv_vec[i];
                    }
                    // send the topology i know to next neighbour on inel
                    for (int j = 0; j < nProcesses; j++)
                        send_vec[j] = parents[j];
                    MPI_Send(send_vec, nProcesses, MPI_INT, neigh[1], 0, MPI_COMM_WORLD);
                    fprintf(stdout, "M(%d,%d)\n", rank, neigh[1]);
                }
            }
        } else {
            // There is no cycle and we can use a dfs approch where we send to each neigh only if that neighbour
            // isn't the one we received information from

            // sync info about cluster
            // make a dfs send from cluster ROOT(0)
            if (rank == ROOT) {
                int neighbour = -1;
                for (int i = 0; i < 2; i++)
                    if (neigh[i] != -1)
                        neighbour = neigh[i];
                for (int j = 0; j < nProcesses; j++)
                    send_vec[j] = parents[j];
                
                if (neighbour == -1)
                    fprintf(stderr, "Rank %d try to sent to neighbour -1 in sync_topology\n", rank);

                // sent curent topology to only neighbour
                MPI_Send(send_vec, nProcesses, MPI_INT, neighbour, neighbour, MPI_COMM_WORLD);
                fprintf(stdout, "M(%d,%d)\n", rank, neighbour);

                // await the complete answer from neighbour after dfs
                MPI_Status status;
                MPI_Recv(recv_vec, nProcesses, MPI_INT, neighbour, rank, MPI_COMM_WORLD, &status);

                //update my known topology
                for (int i = 0; i < nProcesses; i++) {
                    if (parents[i] == -1 && recv_vec[i] != -1)
                        parents[i] = recv_vec[i];
                }
            } else {
                if (!(commError == 2 && rank == 1)) {
                    // await topology from parent in dfs from ROOT
                    MPI_Status status;
                    MPI_Recv(recv_vec, nProcesses, MPI_INT, MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);
                    int father = status.MPI_SOURCE;

                    //update my known topology
                    for (int i = 0; i < nProcesses; i++) {
                        if (parents[i] == -1 && recv_vec[i] != -1)
                            parents[i] = recv_vec[i];
                    }
                    
                    // send topology to next unvisited neighbour
                    int neighbour = -1;
                    for (int i = 0; i < 2; i++) {
                        if (neigh[i] == -1 || neigh[i] == status.MPI_SOURCE)
                            continue;
                        for (int j = 0; j < nProcesses; j++)
                            send_vec[j] = parents[j];
                        neighbour = neigh[i];
                        MPI_Send(send_vec, nProcesses, MPI_INT, neigh[i], neigh[i], MPI_COMM_WORLD);
                        fprintf(stdout, "M(%d,%d)\n", rank, neigh[i]);
                    }

                    // await answer with full topology from neighbour if it exists
                    if (neighbour != -1)
                        MPI_Recv(recv_vec, nProcesses, MPI_INT, neighbour, rank, MPI_COMM_WORLD, &status);

                    //update my topology fathers array
                    for (int i = 0; i < nProcesses; i++) {
                        if (parents[i] == -1 && recv_vec[i] != -1)
                            parents[i] = recv_vec[i];
                    }
                    
                    // send complete topology back to father
                    for (int j = 0; j < nProcesses; j++)
                        send_vec[j] = parents[j];
                    MPI_Send(send_vec, nProcesses, MPI_INT, father, father, MPI_COMM_WORLD);
                    fprintf(stdout, "M(%d,%d)\n", rank, father);
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void print_topology() {
    fprintf(stdout, "%d -> ", rank);
    for (int i = 0; i < NUM_CLUSTERS; i++) {
        int nr = 0;
        for (int j = 0; j < nProcesses; j++)
            if (parents[j] == i)
                nr++;
        if (nr != 0) {
            fprintf(stdout, "%d:", i);
            for (int j = 0; j < nProcesses; j++) {
                if (parents[j] == i) {
                    nr--;
                    fprintf(stdout, "%d", j);
                    if (nr)
                        fprintf(stdout, ",");
                    else
                        fprintf(stdout, " ");
                }
            }
        }
    }
    fprintf(stdout, "\n");
}

void inform_workers_and_show_topology() {
    if (rank < NUM_CLUSTERS) {
        int send_vec[nProcesses];
        // print topology
        print_topology();

        // inform workers about toology
        for (int j = 0; j < nProcesses; j++)
            send_vec[j] = parents[j];
       
        for (int i = 0; i < nWorkers; i++) {
            MPI_Send(send_vec, nProcesses, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);
        }
    } else {
        int recv_vec[nProcesses];
        MPI_Status status;
        MPI_Recv(recv_vec, nProcesses, MPI_INT, MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);

        // set my topology fathers array
        for (int i = 0; i < nProcesses; i++)
                parents[i] = recv_vec[i];
        
        // print topology
        print_topology();
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

int get_number_of_workers() {
    int nr = 0;
    for (int i = 0; i < nProcesses; i++) {
        if (parents[i] != -1)
            nr++;
        if (parents[i] == 1 && commError == 2)
            nr--;
    }
    return nr;
}

void compute_array_no_error() {
    int length, numWorkers;
    int *arr;
    if (rank < NUM_CLUSTERS) {
        // algorithm inel for communication

        // First phase the vector is sent beetween clusters to be proccessed by workers, as well as start positions
        if (rank == ROOT) {
            // ROOT allocate the vector and initiate the communication in inel
            int *array = (int *) malloc(arraySize * sizeof(int));
            if (array == NULL)
                fprintf(stderr, "Error allocating array to be porcessed.\n");
            for (int i = 0; i < arraySize; i++)
                array[i] = arraySize - i - 1;

            length = arraySize;
            MPI_Send(&length, 1, MPI_INT, neigh[1], 0, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neigh[1]);
            MPI_Send(array, length, MPI_INT, neigh[1], 0, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neigh[1]);

            MPI_Status status;
            MPI_Recv(&length, 1, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);
            arr = (int *) malloc(length * sizeof(int));
            MPI_Recv(arr, length, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);
        } else {
            MPI_Status status;
            MPI_Recv(&length, 1, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);
            arr = (int *) malloc(length * sizeof(int));
            MPI_Recv(arr, length, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);

            MPI_Send(&length, 1, MPI_INT, neigh[1], 0, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neigh[1]);
            MPI_Send(arr, length, MPI_INT, neigh[1], 0, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neigh[1]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
        
    if (rank < NUM_CLUSTERS) {
        // Second phase assing workers start position
        // fprintf(stdout, "\nSecond phase start\n");
        if (rank == ROOT) {
            // ROOT allocate the vector and initiate the communication in inel
            //assign work to each worker (start position in array)
            numWorkers = get_number_of_workers();
            int workPayloadSize = arraySize / numWorkers;
            workStart = (int *)calloc(numWorkers, sizeof(int));
            for (int i = 1; i < numWorkers; i++)
                workStart[i] = workStart[i - 1] + workPayloadSize;

            MPI_Send(&numWorkers, 1, MPI_INT, neigh[1], 0, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neigh[1]);
            MPI_Send(workStart, numWorkers, MPI_INT, neigh[1], 0, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neigh[1]);

            MPI_Status status;
            MPI_Recv(&numWorkers, 1, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);
            MPI_Recv(workStart, numWorkers, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);
        } else {
            MPI_Status status;
            MPI_Recv(&numWorkers, 1, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);
            workStart = (int *)calloc(numWorkers, sizeof(int));
            MPI_Recv(workStart, numWorkers, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);

            MPI_Send(&numWorkers, 1, MPI_INT, neigh[1], 0, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neigh[1]);
            MPI_Send(workStart, numWorkers, MPI_INT, neigh[1], 0, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neigh[1]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank < NUM_CLUSTERS) {
        // Third phase each cluster send his part of the array to be porecessed by his workers
        // fprintf(stdout, "\nThird phase start\n");
        int pos = 0;
        for (int i = 0; i < nProcesses; i++) {
            if (parents[i] < rank && parents[i] != -1)
                pos++;
            if (parents[i] == 1 && commError == 2)
                pos--;
        }

        int send_arr[length], recv_arr[length];
        for (int i = 0; i < nWorkers; i++) {
            // send array to workers and where to process
            for (int j = 0; j < length; j++)
                send_arr[j] = arr[j];
            int start, end;
            MPI_Send(&length, 1, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);
            MPI_Send(send_arr, length, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);
            
            start = workStart[pos + i];
            MPI_Send(&start, 1, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);

            end = (pos + i + 1 < numWorkers) ? workStart[pos + i + 1] : length;
            MPI_Send(&end, 1, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);

            MPI_Status status;
            MPI_Recv(recv_arr, length, MPI_INT, worker[i], rank, MPI_COMM_WORLD, &status);
            for (int j = start; j < end; j++)
                arr[j] = recv_arr[j];
        }
    } else {
        int start, end;
        MPI_Status status;
        MPI_Recv(&length, 1, MPI_INT, parents[rank], rank, MPI_COMM_WORLD, &status);
        arr = (int *) malloc(length * sizeof(int));
        MPI_Recv(arr, length, MPI_INT, parents[rank], rank, MPI_COMM_WORLD, &status);
        MPI_Recv(&start, 1, MPI_INT, parents[rank], rank, MPI_COMM_WORLD, &status);
        MPI_Recv(&end, 1, MPI_INT, parents[rank], rank, MPI_COMM_WORLD, &status);

        for (int i = start; i < end; i++)
            arr[i] *= 5;
        MPI_Send(arr, length, MPI_INT, parents[rank], parents[rank], MPI_COMM_WORLD);
        fprintf(stdout, "M(%d,%d)\n", rank, parents[rank]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank < NUM_CLUSTERS) {
        // Fourth phase reassemble the proccessed vector by clusters
        int send_arr[length], recv_arr[length];
        int pos = 0;
        for (int i = 0; i < nProcesses; i++) {
            if (parents[i] < rank && parents[i] != -1)
                pos++;
            if (parents[i] == 1 && commError == 2)
                pos--;
        }

        int start = workStart[pos];
        int end = (pos + nWorkers < numWorkers) ? workStart[pos + nWorkers] : length;

        // fprintf(stdout, "%d -> ", rank);
        // for (int i = start; i < end; i++)
        //         fprintf(stdout, "%d ", arr[i]);
        // fprintf(stdout, "\n");
        if (rank == ROOT) {
            // ROOT allocate the vector and initiate the communication in inel
            for (int i = 0; i < length; i++)
                send_arr[i] = -1;
            for (int i = start; i < end; i++)
                send_arr[i] = arr[i];
            MPI_Send(send_arr, length, MPI_INT, neigh[1], 0, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neigh[1]);

            MPI_Status status;
            MPI_Recv(recv_arr, length, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);
            for (int i = 0; i < length; i++)
                arr[i] = recv_arr[i];

            //print answer
            fprintf(stdout, "Rezultat: ");
            for (int i = 0; i < length; i++)
                fprintf(stdout, "%d ", arr[i]);
            fprintf(stdout, "\n");
        } else {
            for (int i = 0; i < length; i++)
                send_arr[i] = -1;
            for (int i = start; i < end; i++)
                send_arr[i] = arr[i];

            MPI_Status status;
            MPI_Recv(recv_arr, length, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);
            for (int i = 0; i < length; i++)
                arr[i] = recv_arr[i];

            for (int i = 0; i < length; i++)
                if (send_arr[i] == -1 && arr[i] != -1)
                    send_arr[i] = arr[i];
            MPI_Send(send_arr, length, MPI_INT, neigh[1], 0, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neigh[1]);
        }

    }
    
    MPI_Barrier(MPI_COMM_WORLD);
}

void print_array(int *arr, int len) {
    fprintf(stdout, "%d -> ", rank);
    for (int i = 0; i < len; i++)
        fprintf(stdout, "%d ", arr[i]);
    fprintf(stdout, "\n");
}

void compute_array_error() {
    int length, numWorkers;
    int *arr;
    if (rank < NUM_CLUSTERS) {
        if (rank == ROOT) {
            // ROOT allocate the vector and initiate the communication in dfs
            int *arr = (int *) malloc(arraySize * sizeof(int));
            for (int i = 0; i < arraySize; i++)
                arr[i] = arraySize - i - 1;
            length = arraySize;

            int neighbour = -1;
            for (int i = 0; i < 2; i++)
                if (neigh[i] != -1)
                    neighbour = neigh[i];
            
            if (neighbour == -1)
                fprintf(stderr, "Rank %d try to sent to neighbour -1 in sync_topology\n", rank);

            // sent length and the array which need to be processed
            MPI_Send(&length, 1, MPI_INT, neighbour, neighbour, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neighbour);
            MPI_Send(arr, length, MPI_INT, neighbour, neighbour, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neighbour);

            // assign work for all workers
            numWorkers = get_number_of_workers();
            int workPayloadSize = arraySize / numWorkers;
            workStart = (int *)calloc(numWorkers, sizeof(int));
            for (int i = 1; i < numWorkers; i++)
                workStart[i] = workStart[i - 1] + workPayloadSize;

            MPI_Send(&numWorkers, 1, MPI_INT, neighbour, neighbour, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neighbour);
            MPI_Send(workStart, numWorkers, MPI_INT, neighbour, neighbour, MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, neighbour);

            // send computation orders to workers
            int pos = 0;
            for (int i = 0; i < nProcesses; i++) {
                if (parents[i] < rank && parents[i] != -1)
                    pos++;
                if (parents[i] == 1 && commError == 2)
                    pos--;
            }

            int send_arr[length], recv_arr[length], st, dr;
            for (int i = 0; i < nWorkers; i++) {
                // send array to workers and where to process
                for (int j = 0; j < length; j++)
                    send_arr[j] = arr[j];
                int start, end;
                MPI_Send(&length, 1, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
                fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);
                MPI_Send(send_arr, length, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
                fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);
                
                start = workStart[pos + i];
                MPI_Send(&start, 1, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
                fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);

                end = (pos + i + 1 < numWorkers) ? workStart[pos + i + 1] : length;
                MPI_Send(&end, 1, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
                fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);

                if (i == 0)
                    st = start;
                if (i == nWorkers - 1)
                    dr = end;

                // await answer from worker
                MPI_Status status;
                MPI_Recv(recv_arr, length, MPI_INT, worker[i], rank, MPI_COMM_WORLD, &status);
                for (int j = start; j < end; j++)
                    arr[j] = recv_arr[j];
            }

            // await the complete answer from neighbour after dfs
            MPI_Status status;
            MPI_Recv(recv_arr, length, MPI_INT, neighbour, rank, MPI_COMM_WORLD, &status);

            //fprintf(stdout, "rank(%d) -> st: %d, dr: %d\n", rank, st, dr);

            // update my answer
            for (int i = 0; i < length; i++) {
                if (i < st || i >= dr)
                    arr[i] = recv_arr[i];
            }

            // print answer
            fprintf(stdout, "Rezultat: ");
            for (int i = 0; i < length; i++)
                fprintf(stdout, "%d ", arr[i]);
            fprintf(stdout, "\n");
        } else {
            if (!(commError == 2 && rank == 1)) {
                // await topology from parent in dfs from ROOT
                MPI_Status status;
                MPI_Recv(&length, 1, MPI_INT, MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);
                int father = status.MPI_SOURCE;
                arr = (int *) malloc(length * sizeof(int));
                MPI_Recv(arr, length, MPI_INT, MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);

                MPI_Recv(&numWorkers, 1, MPI_INT, MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);
                workStart = (int *)calloc(numWorkers, sizeof(int));
                MPI_Recv(workStart, numWorkers, MPI_INT, MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);

                
                // send the array which need to be proccessed to the next neighbour
                int neighbour = -1;
                for (int i = 0; i < 2; i++) {
                    if (neigh[i] == -1 || neigh[i] == father)
                        continue;
                    neighbour = neigh[i];
                }

                if (neighbour != -1) {
                    // send length and the array which need to be processed
                    MPI_Send(&length, 1, MPI_INT, neighbour, neighbour, MPI_COMM_WORLD);
                    fprintf(stdout, "M(%d,%d)\n", rank, neighbour);
                    MPI_Send(arr, length, MPI_INT, neighbour, neighbour, MPI_COMM_WORLD);
                    fprintf(stdout, "M(%d,%d)\n", rank, neighbour);

                    // send number of total workers and work distribution to next cluster
                    MPI_Send(&numWorkers, 1, MPI_INT, neighbour, neighbour, MPI_COMM_WORLD);
                    fprintf(stdout, "M(%d,%d)\n", rank, neighbour);
                    MPI_Send(workStart, numWorkers, MPI_INT, neighbour, neighbour, MPI_COMM_WORLD);
                    fprintf(stdout, "M(%d,%d)\n", rank, neighbour);
                }

                // send the computation orders to workers
                int pos = 0;
                for (int i = 0; i < nProcesses; i++) {
                    if (parents[i] < rank && parents[i] != -1)
                        pos++;
                    if (parents[i] == 1 && commError == 2)
                        pos--;
                }

                int send_arr[length], recv_arr[length], st = -1, dr = -1;
                for (int i = 0; i < nWorkers; i++) {
                    // send array to workers and where to process
                    for (int j = 0; j < length; j++)
                        send_arr[j] = arr[j];
                    int start, end;
                    MPI_Send(&length, 1, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
                    fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);
                    MPI_Send(send_arr, length, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
                    fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);
                    
                    start = workStart[pos + i];
                    MPI_Send(&start, 1, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
                    fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);

                    end = (pos + i + 1 < numWorkers) ? workStart[pos + i + 1] : length;
                    MPI_Send(&end, 1, MPI_INT, worker[i], worker[i], MPI_COMM_WORLD);
                    fprintf(stdout, "M(%d,%d)\n", rank, worker[i]);

                    if (i == 0)
                        st = start;
                    
                    if (i == nWorkers - 1)
                        dr = end;

                    // await answer from worker
                    MPI_Status status;
                    MPI_Recv(recv_arr, length, MPI_INT, worker[i], rank, MPI_COMM_WORLD, &status);
                    for (int j = start; j < end; j++)
                        arr[j] = recv_arr[j];

                }

                // await answer from neighbour if it exists
                if (neighbour != -1)
                    MPI_Recv(recv_arr, length, MPI_INT, neighbour, rank, MPI_COMM_WORLD, &status);

                //fprintf(stdout, "rank(%d) -> st: %d, dr: %d\n", rank, st, dr);

                // update my answer
                if (neighbour != -1) {
                    for (int i = 0; i < length; i++) {
                        if (i < st || i >= dr)
                            arr[i] = recv_arr[i];
                    }
                }
                
                // send answer back to father
                for (int j = 0; j < length; j++)
                    send_arr[j] = arr[j];
                MPI_Send(send_arr, length, MPI_INT, father, father, MPI_COMM_WORLD);
                fprintf(stdout, "M(%d,%d)\n", rank, father);
            }
        }
    } else {
        if (!(parents[rank] == 1 && commError == 2)) {
            int start, end;
            MPI_Status status;
            MPI_Recv(&length, 1, MPI_INT, parents[rank], rank, MPI_COMM_WORLD, &status);
            arr = (int *) malloc(length * sizeof(int));
            MPI_Recv(arr, length, MPI_INT, parents[rank], rank, MPI_COMM_WORLD, &status);
            MPI_Recv(&start, 1, MPI_INT, parents[rank], rank, MPI_COMM_WORLD, &status);
            MPI_Recv(&end, 1, MPI_INT, parents[rank], rank, MPI_COMM_WORLD, &status);

            //print_array(arr, length);
            //fprintf(stdout, "rank(%d) -> st: %d, dr: %d\n", rank, start, end);
            for (int i = start; i < end; i++)
                arr[i] *= 5;
            //print_array(arr, length);

            MPI_Send(arr, length, MPI_INT, parents[rank], parents[rank], MPI_COMM_WORLD);
            fprintf(stdout, "M(%d,%d)\n", rank, parents[rank]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

	if (argc != 3) {
		printf("please run with: mpirun --oversubscribe -np <nr_procese> %s <dim_vec> <eroare_comunicare>\n", argv[0]);
		MPI_Finalize();	
		exit(0);
	}

    if (rank == ROOT)
        arraySize = atoi(argv[1]);
    commError = atoi(argv[2]);

    // Read each main cluster's workers
    read_workers();

    // Create the topology between the main 4 cluster
    create_cluster_neighs();

    // All clusters initiate the parents vector only with the info of their workers
    init_topology();

    // Sync the topology between clusters with info from each clusters workers subordonated
    sync_topology();

    // Now all clusters know the topology and they can inform their workers an each can print the topology they know
    inform_workers_and_show_topology();
    
    // Compute the results on array
    if (commError == 0)
        compute_array_no_error();
    else
        compute_array_error();
    
    MPI_Finalize();
}
