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

int globSum = 0;
int *taskList = nullptr;
int numTasks;
int currTask;


void debug(char* debugMessage) {
    std::cout << debugMessage << std::endl;
}

void makeTaskList (int* taskList, int sizeList, int numProc, int rankProc, int iterCounter) {
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
void* doWork(void* args) {
    int numProc;
    int rankProc;
    MPI_Status* status;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankProc);
    //flag for get additional job for curr proccess
    bool myWorkDone = true;
    for (int iterCounter = 0; iterCounter < NUM_LISTS; ++iterCounter) {
        numTasks = NUM_TASKS / numProc;
        //fill list of tasks
        makeTaskList(taskList, numTasks, numProc, rankProc, iterCounter);
        while (true) {
            for (currTask = 0; currTask < numTasks; ++currTask) {
                doTask(taskList[currTask]);
                //delete task
                taskList[currTask] = 0;
            }
            //done all work in curr list
            numTasks = 0;
            //request for additional job
            for (int i = 0; i < numProc; ++i) {
                MPI_Send(&myWorkDone, 1, MPI_CXX_BOOL, i, TAG_FOR_REQUEST, MPI_COMM_WORLD);
                MPI_Recv(&numTasks, 1, MPI_INT, i, TAG_FOR_NUM, MPI_COMM_WORLD, status);
                MPI_Recv(taskList, numTasks, MPI_INT, i, TAG_FOR_LIST, MPI_COMM_WORLD, status);
            }
            //if no one get job to process, end the loop and curr list of task
            if (numTasks == 0) {
                break;
            }
        }
        delete[] taskList;
    }

    return nullptr;
}


void* waitForRequest(void* args) {
    int restNumTasks = numTasks - currTask + 1;

    return nullptr;
}

int execute() {
    pthread_attr_t attr;
    pthread_t threads[NUM_THREADS];
    if(pthread_attr_init(&attr)){
        perror("Cannot initialize attributes");
        return -1;
    }
    if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE) !=0) {
        perror("lab5 : pthread_attr_setdetachstate() failed");
        return -1;
    }
    if (pthread_create(&threads[0], &attr, doWork, nullptr) != 0) {
        perror("lab5 : pthread_create() failed for 1st thread");
        return -1;
    }
    if (pthread_create(&threads[1], &attr, waitForRequest, nullptr) != 0) {
        perror("lab5 : pthread_create() failed for 2nd thread");
        return  -1;
    }

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
    if (provided != MPI_THREAD_MULTIPLE) {
        perror("lab5 : Not required lvl.\n");
        return  -1;
    }

    execute();

    //TODO Функция MPI_Finalize() должна вызываться только тем потоком,
    // который инициализировал MPI, и только после того, как все потоки
    // завершили выполнение MPI-операций. Поэтому все потоки должны быть
    // присоединяемыми (joinable) и присоединены к этому главному потоку
    // перед вызовом MPI_Finalize().
    MPI_Finalize();
    return 0;
}
