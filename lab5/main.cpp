#include <iostream>
#include <pthread.h>
#include <mpi.h>
#include <cmath>

#define NUM_THREADS 2
#define NUM_LISTS 3
//this parameter we choose by yourself
#define L 12
#define NUM_TASKS 56

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
    for (int i = 0; i < 8; ++i) {
        if (rank == i) {
            std::cout << "NUM " << rank << ":" << debugMessage << std::endl;
        }
        //MPI_Barrier(MPI_COMM_WORLD);
    }
}

void printTasks(int rank) {
    for (int i = 0; i < 8; ++i) {
        if (rank == i) {

            for (int j = 0; j < numTasks; ++j) {
                std::cout << taskList[j] << " ";
            }
            std::cout << std::endl;
        }
        //MPI_Barrier(MPI_COMM_WORLD);

    }
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
    //debug("DOWORK :: STARTING THE LOOP FOR LISTS", rankProc);
    for (int iterCounter = 0; iterCounter < NUM_LISTS; ++iterCounter) {
        double start = MPI_Wtime();
        int sumNumTaskPerIter = 0;
        numTasks = NUM_TASKS / numProc;
        int recvTask = 0;
        taskList = new int[numTasks];
        //fill list of tasks
        //debug("DOWORK :: MAKE LIST ", rankProc);
        makeTaskList(numTasks, numProc, rankProc, iterCounter);
        printTasks(rankProc);
        //debug("DOWORK :: MAKE LIST END", rankProc);
        //debug("DOWORK :: START LOOP WHILE TRUE", rankProc);

        while (true) {
            for (currTask = 0; currTask < numTasks; ++currTask) {
                //debug("DOWORK :: MUTEX LOCK", rankProc);
                pthread_mutex_lock(&mutex);
                //debug("DOWORK :: MUTEX LOCK END", rankProc);
                int currWeight = taskList[currTask];
               // debug("DOWORK :: MUTEX UNLOCK", rankProc);
                pthread_mutex_unlock(&mutex);
                //debug("DOWORK :: MUTEX UNLOCK END", rankProc);
               // debug("DOWORK :: DO TASK START", rankProc);
                doTask(currWeight);
                sumNumTaskPerIter += currWeight;
             //   debug("DOWORK ::DO TASK END", rankProc);
                //delete task

            }
            //done all work in curr list
            //debug("DOWORK :: DONE ALL WORK IN YOURSELF LIST  AND REQUEST FOR JOB", rankProc);
            //numTasks = 0;
            //request for additional job

            for (int i = 0; i < numProc; ++i) {
                if (i != rankProc) {
                    //debug("DOWORK :: START SEND REQUEST", rankProc);
                    MPI_Send(&myWorkDone, 1, MPI_CXX_BOOL, i, TAG_FOR_REQUEST, MPI_COMM_WORLD);
                    //debug("DOWORK :: END SEND REQUEST", rankProc);
                    //debug("DOWORK :: START RECV NUM TASK", rankProc);
                    MPI_Recv(&recvTask, 1, MPI_INT, i, TAG_FOR_NUM, MPI_COMM_WORLD, &status);
                  //  debug("DOWORK :: END RECV NUM TASK", rankProc);
                    //if the proc have nothing to give
                //    debug("DOWORK :: CHECK IF OTHER PROC HAVE JOB", rankProc);
                    if (recvTask == 0) {
                        continue;
                    }
              //      debug("DOWORK :: START RECV TASKS", rankProc);
                    MPI_Recv(taskList, recvTask, MPI_INT, i, TAG_FOR_LIST, MPI_COMM_WORLD, &status);
            //        debug("DOWORK :: END RECV TASKS", rankProc);
          //          debug("DOWORK :: PRINT", rankProc);
                    printTasks(rankProc);
                }
            }
            //if no one get job to process, end the loop and curr list of task
            if (recvTask == 0) {
                break;
            }
            numTasks = recvTask;
        }
        //debug("DOWORK :: END LOOP WHILE TRUE", rankProc);
        delete[] taskList;
        //debug("DOWORK :: START SUNC", rankProc);
        double end = MPI_Wtime();
        double timeOneIter = end - start;
        std::cout << "Time per iter on proc with rank number " << rankProc << ":" << timeOneIter <<std::endl;
        std::cout << "Num of tasks for all iter on proc with rank number " << rankProc << ":" << sumNumTaskPerIter <<std::endl;
        std::cout << "Global sum on proc with rank number " << rankProc << ":" << globSum <<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        //debug("DOWORK :: END SYNC", rankProc);
    }
    //debug("DOWORK :: END THE LOOP FOR LISTS", rankProc);
    myWorkDone = STOP_RECV;
    //debug("DOWORK :: START SEND TO KILL THREAD", rankProc);
    MPI_Send(&myWorkDone, 1, MPI_CXX_BOOL, rankProc, TAG_FOR_REQUEST, MPI_COMM_WORLD);
    //debug("DOWORK :: END SEND TO KILL THREAD", rankProc);

    //debug("DOWORK :: KILL THREAD", rankProc);
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

    //debug("WAIT :: WHILE TRUE LOOP", rankProc);
    while (true) {
        //debug("WAIT :: START RECIEVE FLAG FOR JOB", rankProc);
        MPI_Recv(&giveMeWork, 1, MPI_CXX_BOOL, MPI_ANY_SOURCE, TAG_FOR_REQUEST, MPI_COMM_WORLD, &status);
        //debug("WAIT :: END RECIEVE FLAG FOR JOB", rankProc);
        if (giveMeWork == STOP_RECV) {
        //    debug("WAIT :: KILL THREAD", rankProc);
            pthread_exit(nullptr);
        }
        //debug("WAIT :: START MUTEX LOCK", rankProc);
        pthread_mutex_lock(&mutex);
        //debug("WAIT :: END MUTEX LOCK", rankProc);
        int restNumTasks = numTasks - currTask;
        //debug("WAIT :: PRINT", rankProc);
        numShareTask = restNumTasks / 2;
        //std::cout << "WAIT :: REST NUM " << restNumTasks << " numSHARE " << numShareTask << " RANK: " <<rankProc << std::endl;
        //debug("WAIT :: START SEND NEW NUM OF TASK", rankProc);
        MPI_Send(&numShareTask, 1, MPI_INT, status.MPI_SOURCE, TAG_FOR_NUM, MPI_COMM_WORLD);
        //debug("WAIT :: END SEND NEW NUM OF TASK", rankProc);
        if (numShareTask == 0) {
            pthread_mutex_unlock(&mutex);
            continue;
        }
        numTasks -= numShareTask;
        //debug("WAIT :: START SEND NEW TASKS", rankProc);
        MPI_Send(&taskList[numTasks], numShareTask, MPI_INT, status.MPI_SOURCE, TAG_FOR_LIST, MPI_COMM_WORLD);
        //debug("WAIT :: END SEND NEW TASKS", rankProc);
        //debug("WAIT :: START UNLOCK MUTEX", rankProc);
        pthread_mutex_unlock(&mutex);
        //debug("WAIT :: END UNLOCK MUTEX", rankProc);
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
    //debug("EXEC :: INIT ATTRIBUTES", rankProc);
    if(pthread_attr_init(&attr)){
        perror("Cannot initialize attributes");
        return -1;
    }
   // debug("EXEC :: INIT ATTRIBUTES END", rankProc);
    //do thread joinable
    //debug("EXEC :: DO THREAD JOINABLE", rankProc);
    if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE) !=0) {
        perror("lab5 : pthread_attr_setdetachstate() failed");
        return -1;
    }
    //debug("EXEC :: DO THREAD JOINABLE END", rankProc);
    //create 1st thread
    //debug("EXEC :: CREATE 1ST THREAD", rankProc);
    if (pthread_create(&threads[0], &attr, doWork, nullptr) != 0) {
        perror("lab5 : pthread_create() failed for 1st thread");
        return -1;
    }
    //debug("EXEC :: CREATE 1ST THREAD END", rankProc);
    //create 2nd thread
    //debug("EXEC :: CREATE 2ND THREAD", rankProc);
    if (pthread_create(&threads[1], &attr, waitForRequest, nullptr) != 0) {
        perror("lab5 : pthread_create() failed for 2nd thread");
        return  -1;
    }
    //debug("EXEC :: CREATE 2ND THREAD END", rankProc);
    //main wait the 1st thread

    if (pthread_join(threads[0], nullptr) != 0) {
        perror("lab5 : pthread_join() failed for 1st thread");
        return  -1;
    }
    //main wait the 2nd thread
    if (pthread_join(threads[1], nullptr) != 0) {
        perror("lab5 : pthread_join() failed for 2nd thread");
        return -1;
    }
    //debug("EXEC :: BOTH THREAD DIED", rankProc);
    //clear resources
    pthread_attr_destroy(&attr);

    return 0;
}

int main(int argc, char* argv[]) {
    int provided = 0;
    //init thread
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    ////////////debug///////////
    int numProc;
    int rankProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankProc);
    ////////////////////////////
    debug("INIT THREADS", rankProc);
    if (provided != MPI_THREAD_MULTIPLE) {
        perror("lab5 : Not required lvl.\n");
        //TODO check man for Finalize
        return  -1;
    }
    debug("MAIN :: START EXECUTE", rankProc);
    double start = MPI_Wtime();
    execute();
    double end = MPI_Wtime();
    std::cout << "Program works:" << end -start << std::endl;
    debug("MAIN :: END EXECUTE", rankProc);
    MPI_Finalize();
    return 0;
}
