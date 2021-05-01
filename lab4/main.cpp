#include <iostream>
#include <mpi.h>
#include <cmath>

//initial coords
#define X0 -1
#define Y0 -1
#define Z0 -1
//size of zoned
#define Dx 2
#define Dy 2
#define Dz 2

#define Nx 5
#define Ny 5
#define Nz 5

#define a 100000

MPI_Request sendTop;
MPI_Request sendBottom;

MPI_Request recvTop;
MPI_Request recvBottom;

const double hx = (double)Dx / (Nx - 1);
const double hy = (double)Dy / (Ny - 1);
const double hz = (double)Dz / (Nz - 1);
#define epsilon 0,000001
//debug
void printMatrix(double* array) {
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            printf("%12lf", array[i * Nx + j]);

        }
        std::cout << std::endl;
    }
}

void printCube(double* layers, int partZ, int procNum, int procRank) {
    for (int r = 0; r < procNum; ++r) {
        if (procRank == r) {
            std::cout << "Process number:" << procRank << std::endl;
            for (int i = 0; i < partZ; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    for (int k = 0; k < Nx; ++k) {
                        printf("%12lf", layers[i * Nx * Ny + j * Nx + k]);
                    }
                    std::cout << std::endl;
                }
                std::cout << std:: endl << std::endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

double valuePhi(double x, double y, double z) {
    return x*x + y*y + z*z;
}

double valueRo(double x, double y, double z) {
    return (double)6 - a * valuePhi(x, y, z);
}

double coordX(int x) {
    return X0 + hx * x;
}

double coordY(int y) {
    return Y0 + hy * y;
}

double coordZ(int z) {
    return Z0 + hz * z;
}

double initLayers(double* layers, int offset, int partZ) {
    for (int z = 0; z < partZ; ++z) {
        for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
                if (x == 0 || x == Nx - 1 || y == 0 || y == Ny - 1 || (z + offset) == 0 || (z + offset) == Nz - 1)
                    layers[z * Nx * Ny + y * Nx + x] = valuePhi(coordX(offset + x), coordY(y), coordZ(z));
                else
                    layers[z * Nx * Ny + y * Nx + x] = 0;
            }
        }
    }
}

void getNumParts(int* partSizeZ,int procNum) {
    int partSize = Nz / procNum;
    int numNode = Nz;
    for (int i = 0; i < procNum - 1; ++i) {
        partSizeZ[i] = partSize;
        numNode -=partSize;
    }
    partSizeZ[procNum - 1] = numNode;
}

void getOffset(int* offset, const int* partSizeZ, int procNum) {
    offset[0] = 0;
    for (int i = 1; i < procNum; ++i) {
        offset[i] = partSizeZ[i - 1] + offset[i - 1];
    }
}

void sendBound(double* layer, double* topSlice, double* bottomSlice, int procRank, int procNum, int partZ) {
    if (procRank != 0) {
        MPI_Isend(&(layer[0]), Nx * Ny, MPI_DOUBLE, procRank - 1, 0, MPI_COMM_WORLD, &sendTop);
        MPI_Irecv(topSlice, Nx * Ny, MPI_DOUBLE, procRank - 1, 1, MPI_COMM_WORLD, &recvTop);
    }

    if (procRank != procNum - 1) {
        MPI_Isend(&(layer[(partZ - 1) * Nx * Ny]), Nx * Ny, MPI_DOUBLE, procRank + 1, 1, MPI_COMM_WORLD, &sendBottom);
        MPI_Irecv(bottomSlice, Nx * Ny, MPI_DOUBLE, procRank + 1, 0, MPI_COMM_WORLD, &recvBottom);
    }
}

void waitBounds(int procRank, int procNum) {
    if (procRank != 0) {
        MPI_Wait(&recvTop, MPI_STATUS_IGNORE);
        MPI_Wait(&sendTop, MPI_STATUS_IGNORE);
    }
    if (procRank != procNum - 1) {
        MPI_Wait(&recvBottom, MPI_STATUS_IGNORE);
        MPI_Wait(&sendBottom, MPI_STATUS_IGNORE);
    }
}

void calcCenter(const double* layer, double* buffer, int partZ, int offset) {
    auto factor = (double)1/((double)2 / (hx * hx) + 2 / (hy * hy) + 2 / (hz * hz) + a);
    for (int z = 1; z < partZ - 1; ++z) {
        for (int y = 1; y < Ny - 1; ++y) {
            for (int x = 1; x < Nx - 1; ++x) {
                double fz = (layer[(z - 1) * Nx * Ny + y * Nx + x] + layer[(z + 1) * Nx * Ny + y * Nx + x]) / (hx * hx);
                double fy = (layer[z * Nx * Ny + (y - 1) * Nx + x] + layer[z * Nx * Ny + (y + 1) * Nx + x]) / (hy * hy);
                double fx = (layer[z * Nx * Ny + y * Nx + x - 1] + layer[z * Nx * Ny + y * Nx + x + 1]) / (hx * hx);
                buffer[z * Nx * Ny + y * Nx + x] = (fz +fx + fy - valueRo(coordX(x), coordY(y), coordZ(z + offset))) * factor;
            }
        }
    }
}

void calcEdges(const double* layer, double* buffer, const double* topSlice, double* bottomSlice, int partZ, int offset, int procRank, int procNum) {
    auto factor = (double)1/((double)2 / (hx * hx) + 2 / (hy * hy) + 2 / (hz * hz) + a);
    for (int y = 1; y < Ny - 1; ++y) {
        for (int x = 1; x < Nx; ++x) {
            int z = 0;
            if (procRank !=0) {
                double fz = (topSlice[y * Nx + x] + layer[(z + 1) * Nx * Ny + y * Nx + x]) / (hx * hx);
                double fy = (layer[z * Nx * Ny + (y - 1) * Nx + x] + layer[z * Nx * Ny + (y + 1) * Nx + x]) / (hy * hy);
                double fx = (layer[z * Nx * Ny + y * Nx + x - 1] + layer[z * Nx * Ny + y * Nx + x + 1]) / (hx * hx);
                buffer[z * Nx * Ny + y * Nx + x] = (fz +fx + fy - valueRo(coordX(x), coordY(y), coordZ(z + offset))) * factor;
            }
            z = partZ - 1;
            if (procRank != procNum - 1) {
                double fz = (layer[(z - 1) * Nx * Ny + y * Nx + x] + bottomSlice[y * Nx + x]) / (hx * hx);
                double fy = (layer[z * Nx * Ny + (y - 1) * Nx + x] + layer[z * Nx * Ny + (y + 1) * Nx + x]) / (hy * hy);
                double fx = (layer[z * Nx * Ny + y * Nx + x - 1] + layer[z * Nx * Ny + y * Nx + x + 1]) / (hx * hx);
                buffer[z * Nx * Ny + y * Nx + x] = (fz +fx + fy - valueRo(coordX(x), coordY(y), coordZ(z + offset))) * factor;
            }
        }
    }
}

void methodJacobi() {
    //steps for next node
    int procNum;
    int procRank;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    //size for part of Z line
    int* partSizeZ = new int[procNum];
    int* offset = new int[procNum];
    getNumParts(partSizeZ, procNum);
    getOffset(offset, partSizeZ, procNum);
    auto* layer = new double[partSizeZ[procRank] * Nx * Ny];
    auto* bufferLayer = new double[partSizeZ[procRank] * Nx * Ny];
    auto* topSlice = new double [Nx * Ny];
    auto* bottomSlice = new double [Nx * Ny];
    initLayers(layer, offset[procRank], partSizeZ[procRank]);
    printCube(layer, partSizeZ[procRank], procNum, procRank);
    double maxDiff;
    do {
        sendBound(layer, topSlice, bottomSlice, procRank, procNum, partSizeZ[procRank]);
        calcCenter(layer, bufferLayer, partSizeZ[procRank], offset[procRank]);
        waitBounds(procRank, procNum);
        calcEdges(layer, bufferLayer, topSlice, bottomSlice, partSizeZ[procRank], offset[procRank], procRank, procNum);
        //debug
//        if (procRank == 3) {
//            std::cout << "topSlice" << "\n";
//            printMatrix(topSlice);
//            std::cout << "bottomSlice" << "\n";
//            printMatrix(bottomSlice);
//        }
    }
    while (maxDiff < epsilon);
    delete[] partSizeZ;
    delete[]layer;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    methodJacobi();
    MPI_Finalize();
    return 0;
}
