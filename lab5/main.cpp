#include <iostream>
#include <pthread.h>
#include <mpi.h>
#include <cmath>

#define NUM_THREADS 2
#define NUM_LISTS 3
//this parameter we choose by yourself
#define L 12
#define NUM_TASKS 3600

#define TAG_FOR_LIST 2
#define TAG_FOR_REQUEST 1
#define TAG_FOR_NUM 0

#define STOP_RECV 0

int globSum = 0;
int *taskList = nullptr;
int numTasks;
int currTask;
pthread_mutex_t mutex;

void debug(const std::string& debugMessage, int rank) {
    std::cout << "NUM " << rank << ":" << debugMessage << std::endl;
}

void makeTaskList (int sizeList, int numProc, int rankProc, int iterCounter) {
    for (int i = 0; i < sizeList; ++i) {
        taskList[i] = abs((50 - i % 100) * rankProc - (iterCounter % numProc)) * L;
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
        numTasks = NUM_TASKS / numProc;
        taskList = new int[numTasks];
        //fill list of tasks
        makeTaskList(numTasks, numProc, rankProc, iterCounter);
        while (true) {
            for (currTask = 0; currTask < numTasks; ++currTask) {
                pthread_mutex_lock(&mutex);
                int currWait = taskList[currTask];
                pthread_mutex_unlock(&mutex);
                doTask(currWait);
                //delete task

            }
            //done all work in curr list
            numTasks = 0;
            //request for additional job
            for (int i = 0; i < numProc; ++i) {
                if (i != rankProc) {
                    MPI_Send(&myWorkDone, 1, MPI_CXX_BOOL, i, TAG_FOR_REQUEST, MPI_COMM_WORLD);
                    MPI_Recv(&numTasks, 1, MPI_INT, i, TAG_FOR_NUM, MPI_COMM_WORLD, &status);
                    //if the proc done
                    if (numTasks == 0) {
                        continue;
                    }
                    MPI_Recv(taskList, numTasks, MPI_INT, i, TAG_FOR_LIST, MPI_COMM_WORLD, &status);
                }
            }
            //if no one get job to process, end the loop and curr list of task
            if (numTasks == 0) {
                break;
            }
        }
        delete[] taskList;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    myWorkDone = STOP_RECV;
    MPI_Send(&myWorkDone, 1, MPI_CXX_BOOL, rankProc, TAG_FOR_REQUEST, MPI_COMM_WORLD);
    pthread_exit(nullptr);
}


void *waitForRequest(void *args) {
    bool giveMeWork;
    int numShareTask;
    ////////////debug///////////
    int numProc;
    int rankProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankProc);
    ////////////////////////////
    MPI_Status status;
    while (true) {
        MPI_Recv(&giveMeWork, 1, MPI_CXX_BOOL, MPI_ANY_SOURCE, TAG_FOR_REQUEST, MPI_COMM_WORLD, &status);
        if (giveMeWork == STOP_RECV) {
            pthread_exit(nullptr);
        }
        pthread_mutex_lock(&mutex);
        int restNumTasks = numTasks - currTask + 1;
        numShareTask = restNumTasks / 2;
        MPI_Send(&numShareTask, 1, MPI_INT, status.MPI_SOURCE, TAG_FOR_NUM, MPI_COMM_WORLD);
        if (numShareTask == 0) {
            continue;
        }
        numTasks -= numShareTask;
        MPI_Send(&taskList[numTasks], numShareTask, MPI_INT, status.MPI_SOURCE, TAG_FOR_LIST, MPI_COMM_WORLD);
        pthread_mutex_unlock(&mutex);
    }
}

int execute() {
    ////////////debug///////////
    int numProc;
    int rankProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankProc);
    ////////////////////////////
    pthread_attr_t attr;
    pthread_t threads[NUM_THREADS];
    //init attributes
    debug("EXEC :: INIT ATTRIBUTES", rankProc);
    if(pthread_attr_init(&attr)){
        perror("Cannot initialize attributes");
        return -1;
    }
    debug("EXEC :: INIT ATTRIBUTES END", rankProc);
    //do thread joinable
    debug("EXEC :: DO THREAD JOINABLE", rankProc);
    if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE) !=0) {
        perror("lab5 : pthread_attr_setdetachstate() failed");
        return -1;
    }
    debug("EXEC :: DO THREAD JOINABLE", rankProc);
    //create first thread
    debug("EXEC :: CREATE 1ST THREAD", rankProc);
    if (pthread_create(&threads[0], &attr, doWork, nullptr) != 0) {
        perror("lab5 : pthread_create() failed for 1st thread");
        return -1;
    }
    debug("EXEC :: CREATE 1ST THREAD END", rankProc);
    //create 2nd thread
    debug("EXEC :: CREATE 2ND THREAD", rankProc);
    if (pthread_create(&threads[1], &attr, waitForRequest, nullptr) != 0) {
        perror("lab5 : pthread_create() failed for 2nd thread");
        return  -1;
    }
    debug("EXEC :: CREATE 2ND THREAD END", rankProc);
    //
    if (pthread_join(threads[0], nullptr) != 0) {
        perror("lab5 : pthread_join() failed for 1st thread");
        return  -1;
    }
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
    debug("INIT THREADS");
    if (provided != MPI_THREAD_MULTIPLE) {
        perror("lab5 : Not required lvl.\n");
        //TODO check man for Finalize
        return  -1;
    }
    debug("MAIN :: START EXECUTE");
    execute();
    debug("MAIN :: END EXECUTE");
    MPI_Finalize();
    return 0;
}
