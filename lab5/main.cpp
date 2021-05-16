#include <iostream>
#include <pthread.h>
#include <mpi.h>
#include <cmath>

#define NUM_THREADS 2
#define NUM_LISTS 8
//this parameter we choose by yourself
#define L 11000
#define NUM_TASKS 3600

#define TAG_FOR_LIST 2
#define TAG_FOR_REQUEST 1
#define TAG_FOR_NUM 0

#define STOP_RECV 0

double globSum = 0;
int *taskList = nullptr;
int numTasks;
int currTask;
pthread_mutex_t mutex;

void printConclusion(const std::string& message, int rank, double value) {
     std::cout << message << rank << ": " << value << std::endl;
}
void print (double timeOneIter, int numProc, int rankProc, int iterCounter, int sumNumTaskPerIter) {
    double maxTime = 0;
    double minTime = 0;
    MPI_Allreduce(&timeOneIter, &minTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&timeOneIter, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    double delta = maxTime - minTime;
    double percentDelta = delta * 100 / maxTime;
    for (int i = 0; i < numProc; ++i) {
        if (rankProc == i) {
            if (rankProc == 0) {
                std::cout << "                         ITER " << iterCounter << std::endl;
            }
            std::cout << "                         RANK " << rankProc << std::endl;
            printConclusion("Time per iter on proc with rank number ", rankProc, timeOneIter);
            printConclusion("Num of tasks for all iter on proc with rank number ", rankProc, sumNumTaskPerIter);
            printConclusion("Global sum on proc with rank number ", rankProc, globSum);
            printConclusion("Disbalance time on proc with rank number ", rankProc, delta);
            printConclusion("Percentage of disbalance on proc with rank number ", rankProc, percentDelta);
            if (rankProc == numProc - 1) {
                std::cout << "-------------------------------------------------------------" << std::endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
void makeTaskList (int sizeList, int numProc, int rankProc, int iterCounter) {
    for (int i = 0; i < sizeList; ++i) {
        taskList[i] = abs((50 - i % 100) * (rankProc - (iterCounter % numProc))) * L;
    }
}

void doTask(int repeatNum) {
    for (int i = 0; i < repeatNum; ++i) {
        globSum += sin(i);
    }
}

//function for 1st thread
void *doWork(void *args) {
    int numProc;
    int rankProc;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankProc);
    //flag for get additional job for curr proccess
    bool myWorkDone = true;
    for (int iterCounter = 0; iterCounter < NUM_LISTS; ++iterCounter) {
        double start = MPI_Wtime();
        int sumNumTaskPerIter = 0;
        numTasks = NUM_TASKS / numProc;
        int recvTask = 0;
        taskList = new int[numTasks];
        //fill list of tasks
        makeTaskList(numTasks, numProc, rankProc, iterCounter);
        //while for additional work
        while (true) {
            sumNumTaskPerIter += numTasks;
            //do task list
            for (currTask = 0; currTask < numTasks; ++currTask) {
                //mutex for get safety data
                pthread_mutex_lock(&mutex);
                int currWeight = taskList[currTask];
                pthread_mutex_unlock(&mutex);
                doTask(currWeight);
            }
            //done all work in curr list
            //request for additional job
            for (int i = 0; i < numProc; ++i) {
                if (i != rankProc) {
                    //send flag for additional job
                    MPI_Send(&myWorkDone, 1, MPI_CXX_BOOL, i, TAG_FOR_REQUEST, MPI_COMM_WORLD);
                    //get new num of task
                    MPI_Recv(&recvTask, 1, MPI_INT, i, TAG_FOR_NUM, MPI_COMM_WORLD, &status);
                    //if the proc have nothing to give
                    if (recvTask == 0) {
                        continue;
                    }
                    //get new tasks
                    MPI_Recv(taskList, recvTask, MPI_INT, i, TAG_FOR_LIST, MPI_COMM_WORLD, &status);
                }
            }
            //if no one get job to process, end the loop and curr list of task
            if (recvTask == 0) {
                break;
            }
            numTasks = recvTask;
        }
        delete[] taskList;
        double end = MPI_Wtime();
        double timeOneIter = end - start;
        //TODO FUNC FOR PRINT
        print(timeOneIter, numProc, rankProc, iterCounter, sumNumTaskPerIter);
        //sync all process
        MPI_Barrier(MPI_COMM_WORLD);
    }
    //flag for killing thread in wait
    myWorkDone = STOP_RECV;
    //send flag to kill thread in wait
    MPI_Send(&myWorkDone, 1, MPI_CXX_BOOL, rankProc, TAG_FOR_REQUEST, MPI_COMM_WORLD);
    pthread_exit(nullptr);
}


void *waitForRequest(void *args) {
    bool giveMeWork;
    int numShareTask;
    MPI_Status status;
    //while true for waiting requests
    while (true) {
        //start waiting flag from helper
        MPI_Recv(&giveMeWork, 1, MPI_CXX_BOOL, MPI_ANY_SOURCE, TAG_FOR_REQUEST, MPI_COMM_WORLD, &status);
        //if thread get STOP_RECV kill itself
        if (giveMeWork == STOP_RECV) {
            pthread_exit(nullptr);
        }
        //mutex to lock global values bc we will change its
        pthread_mutex_lock(&mutex);
        int restNumTasks = numTasks - currTask;
        numShareTask = restNumTasks / 2;
        //send the num of task for helper
        MPI_Send(&numShareTask, 1, MPI_INT, status.MPI_SOURCE, TAG_FOR_NUM, MPI_COMM_WORLD);
        //if process can't give tasks unlock mutex and start waiting new request
        if (numShareTask == 0) {
            pthread_mutex_unlock(&mutex);
            continue;
        }
        numTasks -= numShareTask;
        //send new tasks for helper
        MPI_Send(&taskList[numTasks], numShareTask, MPI_INT, status.MPI_SOURCE, TAG_FOR_LIST, MPI_COMM_WORLD);
        pthread_mutex_unlock(&mutex);
    }
}

int execute() {
    pthread_attr_t attr;
    pthread_t threads[NUM_THREADS];
    //init attributes
    if(pthread_attr_init(&attr)){
        perror("Cannot initialize attributes");
        return -1;
    }
    //do thread joinable
    if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE) !=0) {
        perror("lab5 : pthread_attr_setdetachstate() failed");
        return -1;
    }
    //create 1st thread
    if (pthread_create(&threads[0], &attr, doWork, nullptr) != 0) {
        perror("lab5 : pthread_create() failed for 1st thread");
        return -1;
    }
    //create 2nd thread
    if (pthread_create(&threads[1], &attr, waitForRequest, nullptr) != 0) {
        perror("lab5 : pthread_create() failed for 2nd thread");
        return  -1;
    }
    //main thread wait the 1st thread
    if (pthread_join(threads[0], nullptr) != 0) {
        perror("lab5 : pthread_join() failed for 1st thread");
        return  -1;
    }
    //main thread wait the 2nd thread
    if (pthread_join(threads[1], nullptr) != 0) {
        perror("lab5 : pthread_join() failed for 2nd thread");
        return -1;
    }
    //clear resources
    pthread_attr_destroy(&attr);
    return 0;
}

int main(int argc, char* argv[]) {
    int provided = 0;
    //init thread
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE) {
        perror("lab5 : Not required lvl.\n");
        //TODO check man for Finalize
        return  -1;
    }
    double start = MPI_Wtime();
    execute();
    double end = MPI_Wtime();
    std::cout << "Program works:" << end -start << std::endl;
    MPI_Finalize();
    return 0;
}
