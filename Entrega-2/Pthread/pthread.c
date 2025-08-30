#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <limits.h>
#include <float.h>

#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1
#define BS 64

double *A, *B, *BT, *C, *R, *resultMatriz, *minA, *minB, *maxA, *maxB, *totalA, *totalB;
int N, T, posA, posB;  
int size;
double MinB =DBL_MAX, MaxB= DBL_MIN, promB = 0.0; 
double MinA= DBL_MAX, MaxA= DBL_MIN, promA = 0.0;
double escalar;
pthread_barrier_t barrier; 

void inicializar_matriz(double *mat, int n, double val, int orden);
void multiplicar_matrices_bloques(double *a, double *bt, double *c, int n, int bs, int filas_por_bloque, int id);
void multiplicar_bloques(double *bloque_a, double *bloque_b, double *bloque_c, int n, int bs);
void transpuesta(double *mat, double *matT, int n, int inicio, int fin);
double dwalltime();

double dwalltime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

void *calculo(void *arg) {
    double actA, actB;
    int desplI;
    int i, j, actI;
    int id = *((int *)arg);
    int chunk = (N + T -1 ) / T;   // Cantidad de filas por hilo (para N filas y T hilos)
    int inicio = id * chunk;
    int fin = (inicio + chunk < N) ? (inicio + chunk) : N;

    for (i = inicio; i < fin; i++){
	actI = i * N;
	for (j = 0; j < N; j++) {
        	if (actA < minA[id]) minA[id] = actA; 
        	if (actA > maxA[id]) maxA[id] = actA;
        	totalA[id] += actA;
 	}
    }

    for (i = inicio; i < fin; i++){
	actI = i * N;
	for (j = 0; j < N; j++) {
        	if (actB < minB[id]) minB[id] = actB; 
        	if (actB > maxB[id]) maxB[id] = actB;
        	totalB[id] += actB;
	}
    }

    pthread_barrier_wait(&barrier);

    // Hilo 0 calcula los valores globales y el escalar
    if (id == 0) {
        for (int i = 0; i < T; i++) {
            if (minA[i] < MinA) MinA = minA[i];
            if (maxA[i] > MaxA) MaxA = maxA[i];
            if (minB[i] < MinB) MinB = minB[i];
            if (maxB[i] > MaxB) MaxB = maxB[i];
            promA += totalA[i];
            promB += totalB[i];
        }
        promA /= size;
        promB /= size;
        escalar = (MaxA * MaxB - MinA * MinB) / (promA * promB);
    }

    // Es necesario esperar si o si porque necestio escalar, por mas que haya muchas operaciones entre medio.
    // No tengo garantizado que escalar se va a calcular antes de que se lo necesite sin la barrera.
    // Creeria que hacerlo aca o al finalizar todas las demas operaciones sobre la matriz es lo mismo, ya que siempre dependo del Hilo 0.
    pthread_barrier_wait(&barrier);

    transpuesta(B, BT, N, inicio, fin);
 
    // Multiplicaciones de matrices
    multiplicar_matrices_bloques(A, B, resultMatriz, N, BS, chunk, id);
    multiplicar_matrices_bloques(C, BT, R, N, BS, chunk, id);
	
    // Aplicar escalar en base a rango lineal
    for (i = inicio; i < fin; i++) {
        desplI = i * N;
        for (j = 0; j < N; j++) {
            R[desplI + j] += resultMatriz[desplI + j] * escalar;
        }
    }   
    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    if (argc != 3 || (N = atoi(argv[1])) <= 0 || (T = atoi(argv[2])) <= 0) {
        printf("Error: Se esperaban 2 argumentos: N (tamaño matriz) y T (cantidad de hilos).\n");
        exit(1);
    }
    N = atoi(argv[1]);
    T = atoi(argv[2]);

    int i;
    double inicioTiempo;
    int ids[T];
    pthread_t threads[T];
    size = N * N;

    A = (double *)malloc(size * sizeof(double));
    B = (double *)malloc(size * sizeof(double));
    BT = (double *)malloc(size * sizeof(double));
    C = (double *)malloc(size * sizeof(double));
    R = (double *)malloc(size * sizeof(double));
    resultMatriz = (double *)malloc(size * sizeof(double));
    minA = (double *)malloc(T * sizeof(double));
    minB = (double *)malloc(T * sizeof(double));
    maxA = (double *)malloc(T * sizeof(double));
    maxB = (double *)malloc(T * sizeof(double));
    totalA = (double *)malloc(T * sizeof(double));
    totalB = (double *)malloc(T * sizeof(double));

    inicializar_matriz(A, N, 1.0, ORDENXFILAS);
    inicializar_matriz(B, N, 1.0, ORDENXCOLUMNAS);
    inicializar_matriz(C, N, 1.0, ORDENXFILAS);
    inicializar_matriz(resultMatriz, N, 0.0, ORDENXFILAS);
    inicializar_matriz(R, N, 0.0, ORDENXFILAS);	    

    for (i = 0; i < T; i++) {
        minA[i] = DBL_MAX;  
        maxA[i] = -DBL_MAX;
        minB[i] = DBL_MAX;
        maxB[i] = -DBL_MAX;
        totalA[i] = 0.0;
        totalB[i] = 0.0;
        ids[i] = i;
    }

    pthread_barrier_init(&barrier, NULL, T);

    inicioTiempo = dwalltime();

    for (i = 0; i < T; i++) {
        pthread_create(&threads[i], NULL, calculo, &ids[i]);
    }

    for (i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
    }

    double totalTiempo = dwalltime() - inicioTiempo;

    printf("Tiempo ejecucion: %f segundos\n", totalTiempo);

    free(A);
    free(B);
    free(BT);
    free(C);
    free(R);
    free(minA);
    free(minB);
    free(maxA);
    free(maxB);
    free(totalA);
    free(totalB);
    free(resultMatriz);
    pthread_barrier_destroy(&barrier);
    
    return 0;
}

void inicializar_matriz(double *mat, int n, double val, int orden) {
    int i, j;

    if (orden == ORDENXFILAS) {
        for (i = 0; i < n; i++) 
            for (j = 0; j < n; j++) 
                mat[i * n + j] = val;
    } else {
        for (i = 0; i < n; i++) 
            for (j = 0; j < n; j++) 
                mat[j * n + i] = val;
    }
}

void transpuesta(double *mat, double *matT, int n, int inicio, int fin) {
    int i, j, actI;

    for (i = inicio; i < fin; i++) {
	actI = i * n;
        for (j = 0; j < n; j++) {
            matT[j * n + i] = mat[actI + j];
	}
    }
}

void multiplicar_matrices_bloques(double *a, double *bt, double *c, int n, int bs, int filas_por_bloque, int id) {
    int i, j, k, comienzo, final, actI, actJ;

    comienzo = id * filas_por_bloque;
    final = comienzo + filas_por_bloque;

    for (i = comienzo; i < final; i += bs) {
	actI = i * n;
        for (j = 0; j < n; j += bs) {
	    actJ = j * n;
            for (k = 0; k < n; k += bs) {
                multiplicar_bloques(&a[actI + k], &bt[actJ + k], &c[actI  + j], n, bs);
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
            suma = 0.0;
	    actJ = j * n;
            for (k = 0; k < bs; k++) {
                suma += bloque_a[actI  + k] * bloque_b[actJ  + k];
            }
            bloque_c[actI  + j] += suma;
        }
    }
}
