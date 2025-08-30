#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpi.h>

#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1
#define BS 64
#define COORDINADOR 0

void inicializar_matriz(double *mat, int n, double val, int orden, int stripSize);
void multiplicar_bloques(double *bloque_a, double *bloque_b, double *bloque_c, int n, int bs);

int main(int argc, char *argv[]) {
    double *A, *B, *BT, *C, *R, *resultMatriz;
    int N, bs, i, j, fila, desplI, posA, posB, desplJ, T, k;
    int numProcs, rank, stripSize, provided;
    double local_min[2] = {INT_MAX, INT_MAX}, local_max[2] = {INT_MIN, INT_MIN}, min[2] = {INT_MAX, INT_MAX}, max[2] = {INT_MIN, INT_MIN};
    double escalar = 0.0;
    MPI_Status status;
    double local_prom[2] = {0, 0}, prom[2] = {0, 0};
    double commTimes[8], maxCommTimes[8], minCommTimes[8], commTime, totalTime;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    N = strtol(argv[1], NULL, 10);
    T = strtol(argv[2], NULL, 10);

    if (argc != 3) {
        if (rank == COORDINADOR) {
            printf("Error: Se esperaban 2 argumentos: N y T (cant hilos por proceso).\n");
        }
        MPI_Finalize();
        return 1;
    }

    if (N % numProcs != 0) {
        if (rank == COORDINADOR) {
            printf("Error: N debe ser multiplo del numero de procesos.\n");
        }
        MPI_Finalize();
        return 1;
    }

    double size = N * N;
    stripSize = N / numProcs;
    double sizeWorker = N * stripSize;
    bs = (N / (numProcs * T) < BS ? N / (numProcs * T) : BS);

    if (rank == COORDINADOR) {
        A = (double *)malloc(size * sizeof(double));
        C = (double *)malloc(size * sizeof(double));
        R = (double *)malloc(size * sizeof(double));
        resultMatriz = (double *)malloc(size * sizeof(double));
    } else {
        A = (double *)malloc(sizeWorker * sizeof(double));
        C = (double *)malloc(sizeWorker * sizeof(double));
        R = (double *)malloc(sizeWorker * sizeof(double));
        resultMatriz = (double *)malloc(sizeWorker * sizeof(double));
    }

    B = (double *)malloc(size * sizeof(double));
    BT = (double *)malloc(size * sizeof(double));

    if (!A || !B || !BT || !C || !R || !resultMatriz) {
        printf("Error: Fallo en la asignacion de memoria.\n");
        MPI_Finalize();
        return 1;
    }

    if (rank == COORDINADOR) {
        inicializar_matriz(A, N, 1.0, ORDENXFILAS, N);
        inicializar_matriz(B, N, 1.0, ORDENXCOLUMNAS, N);
        inicializar_matriz(C, N, 1.0, ORDENXFILAS, N);
        inicializar_matriz(resultMatriz, N, 0.0, ORDENXFILAS, N);
        inicializar_matriz(R, N, 0.0, ORDENXFILAS, N);
    } else {
        inicializar_matriz(R, N, 0.0, ORDENXFILAS, stripSize);
        inicializar_matriz(resultMatriz, N, 0.0, ORDENXFILAS, stripSize);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    commTimes[0] = MPI_Wtime();

    MPI_Scatter(A, sizeWorker, MPI_DOUBLE, A, sizeWorker, MPI_DOUBLE, COORDINADOR, MPI_COMM_WORLD);
    MPI_Scatter(C, sizeWorker, MPI_DOUBLE, C, sizeWorker, MPI_DOUBLE, COORDINADOR, MPI_COMM_WORLD);

    commTimes[1] = MPI_Wtime();

    if (rank == COORDINADOR) {
        for (i = 0; i < N; i++) {
            desplI = i * N;
            for (j = 0; j < N; j++) {
                BT[j * N + i] = B[desplI + j];
            }
        }
    }

    MPI_Bcast(B, size, MPI_DOUBLE, COORDINADOR, MPI_COMM_WORLD);
    MPI_Bcast(BT, size, MPI_DOUBLE, COORDINADOR, MPI_COMM_WORLD);

    #pragma omp parallel num_threads(T) private(i, j, k, fila, desplI, desplJ, posA, posB)
    {
        #pragma omp for reduction(+:local_prom[0]) reduction(max: local_max[0]) reduction(min: local_min[0]) nowait schedule(static)
        for (i = 0; i < stripSize; i++) {
            fila = i * N;
            for (j = 0; j < N; j++) {
                posA = A[fila + j];
                local_min[0] = (posA < local_min[0]) ? posA : local_min[0];
                local_max[0] = (posA > local_max[0]) ? posA : local_max[0];
                local_prom[0] += posA;
            }
        }

        int inicio = stripSize * rank;
        int fin = inicio + stripSize;
        #pragma omp for reduction(+:local_prom[1]) reduction(max: local_max[1]) reduction(min: local_min[1]) schedule(static)
        for (i = inicio; i < fin; i++) {
            fila = i * N;
            for (j = 0; j < N; j++) {
                posB = B[fila + j];
                local_min[1] = (posB < local_min[1]) ? posB : local_min[1];
                local_max[1] = (posB > local_max[1]) ? posB : local_max[1];
                local_prom[1] += posB;
            }
        }

        commTimes[2] = MPI_Wtime();
    }

    MPI_Reduce(&local_min, &min, 2, MPI_DOUBLE, MPI_MIN, COORDINADOR, MPI_COMM_WORLD);
    MPI_Reduce(&local_max, &max, 2, MPI_DOUBLE, MPI_MAX, COORDINADOR, MPI_COMM_WORLD);
    MPI_Reduce(&local_prom, &prom, 2, MPI_DOUBLE, MPI_SUM, COORDINADOR, MPI_COMM_WORLD);

    commTimes[3] = MPI_Wtime();

    if (rank == COORDINADOR) {
        prom[0] = prom[0] / size;
        prom[1] = prom[1] / size;
        escalar = ((max[0] * max[1]) - (min[0] * min[1])) / (prom[0] * prom[1]);
    }

    #pragma omp barrier

    #pragma omp parallel num_threads(T) private(i, j, k, desplI, desplJ)
    {
        #pragma omp for nowait schedule(static)
        for (i = 0; i < stripSize; i += bs) {
            desplI = i * N;
            for (j = 0; j < N; j += bs) {
                desplJ = j * N;
                for (k = 0; k < N; k += bs) {
                    multiplicar_bloques(&A[desplI + k], &B[desplJ + k], &resultMatriz[desplI + j], N, bs);
                }
            }
        }

        #pragma omp for nowait schedule(static)
        for (i = 0; i < stripSize; i += bs) {
            desplI = i * N;
            for (j = 0; j < N; j += BS) {
                desplJ = j * N;
                for (k = 0; k < N; k += bs) {
                    multiplicar_bloques(&C[desplI + k], &BT[desplJ + k], &R[desplI + j], N, bs);
                }
            }
        }

        commTimes[4] = MPI_Wtime();
    }

    MPI_Bcast(&escalar, 1, MPI_DOUBLE, COORDINADOR, MPI_COMM_WORLD);
    commTimes[5] = MPI_Wtime();

    #pragma omp parallel num_threads(T) private(i, j, desplI)
    {
        #pragma omp for nowait schedule(static)
        for (i = 0; i < stripSize; i++) {
            desplI = i * N;
            for (j = 0; j < N; j++) {
                R[desplI + j] += resultMatriz[desplI + j] * escalar;
            }
        }
    }

    commTimes[6] = MPI_Wtime();
    #pragma omp barrier

    MPI_Gather(R, sizeWorker, MPI_DOUBLE, R, sizeWorker, MPI_DOUBLE, COORDINADOR, MPI_COMM_WORLD);
    commTimes[7] = MPI_Wtime();

    MPI_Reduce(commTimes, minCommTimes, 8, MPI_DOUBLE, MPI_MIN, COORDINADOR, MPI_COMM_WORLD);
    MPI_Reduce(commTimes, maxCommTimes, 8, MPI_DOUBLE, MPI_MAX, COORDINADOR, MPI_COMM_WORLD);

    MPI_Finalize();

    if (rank == COORDINADOR) {
        totalTime = maxCommTimes[7] - minCommTimes[0];
        commTime = (maxCommTimes[1] - minCommTimes[0]) +
                   (maxCommTimes[3] - minCommTimes[2]) +
                   (maxCommTimes[5] - minCommTimes[4]) +
                   (maxCommTimes[7] - minCommTimes[6]);
        printf("N=%d  Tiempo total=%lf  Tiempo comunicacion=%lf\n", N, totalTime, commTime);
    }

    free(A);
    free(B);
    free(BT);
    free(C);
    free(R);
    free(resultMatriz);

    return 0;
}

void inicializar_matriz(double *mat, int n, double val, int orden, int stripSize) {
    int i, j, fila;
    if (orden == ORDENXFILAS) {
        for (i = 0; i < stripSize; i++) {
            fila = i * n;
            for (j = 0; j < n; j++) {
                mat[fila + j] = val;
            }
        }
    } else {
        for (i = 0; i < stripSize; i++) {
            for (j = 0; j < n; j++) {
                mat[j * n + i] = val;
            }
        }
    }
}

void multiplicar_matrices_bloques(double *a, double *bt, double *c, int n, int bs, int stripSize) {
    int i, j, k, actI, actJ;
    for (i = 0; i < stripSize; i += bs) {
        actI = i * n;
        for (j = 0; j < n; j += bs) {
            actJ = j * n;
            for (k = 0; k < n; k += bs) {
                multiplicar_bloques(&a[actI + k], &bt[actJ + k], &c[actI + j], n, bs);
            }
        }
    }
}

void multiplicar_bloques(double *bloque_a, double *bloque_b, double *bloque_c, int n, int bs) {
    int i, j, k, actI, actJ;
    double suma;
    for (i = 0; i < bs; i++) {
        actI = i * n;
        for (j = 0; j < bs; j++) {
            actJ = j * n;
            suma = 0.0;
            for (k = 0; k < bs; k++) {
                suma += bloque_a[actI + k] * bloque_b[actJ + k];
            }
            bloque_c[actI + j] += suma;
        }
    }
}
