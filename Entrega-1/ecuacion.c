#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1

// Prototipos de funciones
void inicializar_matriz(double *mat, int n, double val, int orden);
void multiplicar_matrices_bloques(double *a, double *bt, double *c, int n, int bs);
void multiplicar_bloques(double *bloque_a, double *bloque_b, double *bloque_c, int n, int bs);
void transpuesta(double *mat, double *matT, int n);
double dwalltime();

//para calcular tiempo
double dwalltime() {
	double sec;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}

int main(int argc, char *argv[]) {
	double *A, *B, *BT, *C, *R, *resultMatriz;
	int N, BS, i, j, fila;
	double inicioTiempo;
	double minB, maxB, promB = 0.0;
	double minA, maxA, promA = 0.0;
	double posA, posB;

	//validacion de parametros
	if (argc != 3) {
        printf("Error: Se esperaban 2 argumentos, pero se recibieron %d. \n", argc - 1);
    	printf("Se espera un valor N para el tamaÃ±o de la matriz y un valor BS para el tamaÃ±o de los bloques.\n", argv[0]);
    	exit(1);
	}

	N = strtol(argv[1], NULL, 10);
	BS = strtol(argv[2], NULL, 10);

	if ((N <= 0) || (BS <= 0) || (N % BS != 0)) {
    	printf("Parametros mal ingresados. Ejemplo de valores esperados 512 128\n");
    	exit(1);
	}

	//asignacion de memoria
	double size = N * N;
	A = (double *)malloc(size * sizeof(double));
	B = (double *)malloc(size * sizeof(double));
	BT = (double *)malloc(size * sizeof(double));
	C = (double *)malloc(size * sizeof(double));
	R = (double *)malloc(size * sizeof(double));
	resultMatriz = (double *)malloc(size * sizeof(double));

	inicializar_matriz(A, N, 1.0, ORDENXFILAS);
	inicializar_matriz(B, N, 1.0, ORDENXCOLUMNAS);
	inicializar_matriz(C, N, 1.0, ORDENXFILAS);
	inicializar_matriz(resultMatriz, N, 0.0, ORDENXFILAS);
	transpuesta(B, BT, N);

	inicioTiempo = dwalltime();

	// Calcular min, max y promedios de las matrices A y B
	for (i = 0; i < N; i++) {
    	fila = i * N;
    	for (j = 0; j < N; j++) {
            posA = A[fila + j];
            posB = B[fila + j];

            minA = (posA < minA) ? posA : minA;
            maxA = (posA > maxA) ? posA : maxA;
            promA += posA;

            minB = (posB < minB) ? posB : minB;
            maxB = (posB > maxB) ? posB : maxB;
            promB += posB;
    	}
	}

	promA /= size;
	promB /= size;

	if (promA == 0 || promB == 0) {
    	printf("Error: Division por cero \n");
    	exit(1);
	}

	double escalar = ((maxA * maxB - minA * minB) / (promA * promB));

    multiplicar_matrices_bloques(A, B, resultMatriz, N, BS);
    multiplicar_matrices_bloques(C, BT, R, N, BS);

    for (i = 0; i < N * N; i++) {
        R[i] += resultMatriz[i] * escalar;
    }

	double totalTiempo = dwalltime() - inicioTiempo;

	printf("TamaÃ±o de matriz: %d, Bloque: %d, Tiempo: %f segundos\n", N, BS, totalTiempo);

	free(A);
	free(B);
	free(BT);
	free(C);
	free(R);
	free(resultMatriz);

	return 0;
}

// FunciÃ³n para inicializar una matriz
void inicializar_matriz(double *mat, int n, double val, int orden) {
	int i, j;
	if (orden == ORDENXFILAS) {
    	for (i = 0; i < n; i++) {
        	for (j = 0; j < n; j++) {
            	mat[i * n + j] = val;
        	}
    	}
	} else {
    	for (i = 0; i < n; i++) {
        	for (j = 0; j < n; j++) {
            	mat[j * n + i] = val;
        	}
    	}
	}
}

void transpuesta(double *mat, double *matT, int n) {
	int i, j;
	for (i = 0; i < n; i++) {
    	for (j = 0; j < n; j++) {
        	matT[j * n + i] = mat[i * n + j];
    	}
	}
}

// MultiplicaciÃ³n de matrices en bloques
void multiplicar_matrices_bloques(double *a, double *bt, double *c, int n, int bs) {
	int i, j, k;
	for (i = 0; i < n; i += bs) {
    	for (j = 0; j < n; j += bs) {
        	for (k = 0; k < n; k += bs) {
            	multiplicar_bloques(&a[i * n + k], &bt[j * n + k], &c[i * n + j], n, bs);
        	}
    	}
	}
}
// MultiplicaciÃ³n de bloques de matrices
void multiplicar_bloques(double *bloque_a, double *bloque_b, double *bloque_c, int n, int bs) {
    int i, j, k;
    for (i = 0; i < bs; i++) {
        for (j = 0; j < bs; j++) {// J itera sobre las columnas del bloque b 
            for (k = 0; k < bs; k++) { // K itera sobre las filas de BT
                bloque_c[i * n + j] += bloque_a[i * n + k] * bloque_b[j * n + k];
            }
        }
    }
}
