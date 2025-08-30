#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1
#define BS 64

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
	int N, i, j, fila, desplI, posA, posB;
	double inicioTiempo;
	double minB, maxB, promB = 0.0;
	double minA, maxA, promA = 0.0;

	//validacion de parametros
	if (argc != 2) {
        	printf("Error: Se esperaban 2 argumentos, pero se recibieron %d. \n", argc - 1);
    		printf("Se espera un valor N para el tamaño de la matriz y un valor BS para el tamaño de los bloques.\n", argv[0]);
    		exit(1);
	}

	N = strtol(argv[1], NULL, 10);

	if ((N <= 0) || (BS <= 0) || (N % BS != 0)) {
    		printf("Parametros mal ingresados. Ejemplo de valores esperados N: 512, BS: 128\n");
    		exit(1);
	}

	// Asignacion de memoria
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
	inicializar_matriz(R, N, 0.0, ORDENXFILAS);

	inicioTiempo = dwalltime();

	// Calculo transpuesta
	transpuesta(B, BT, N);

	// Variables de A
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

    	for (i = 0; i < size; i++) {
		R[i] += resultMatriz[i] * escalar;	
    	}

	double totalTiempo = dwalltime() - inicioTiempo;

	printf("Tamaño de matriz: %d, Bloque: %d, Tiempo: %f segundos\n", N, BS, totalTiempo);

	free(A);
	free(B);
	free(BT);
	free(C);
	free(R);
	free(resultMatriz);

	return 0;
}

// Función para inicializar una matriz
void inicializar_matriz(double *mat, int n, double val, int orden) {
	int i, j;
	if (orden == ORDENXFILAS) {
    		for (i = 0; i < n; i++) {
        		for (j = 0; j < n; j++) {
            			mat[i * n + j] = val;
        		}
    		}
	} 
	else {
    		for (i = 0; i < n; i++) {
        		for (j = 0; j < n; j++) {
            			mat[j * n + i] = val;
        		}
    		}
	}
}

void transpuesta(double *mat, double *matT, int n) {
	int i, j, actI;

	for (i = 0; i < n; i++) {
		actI = i * n;
    		for (j = 0; j < n; j++) {
        		matT[j * n + i] = mat[actI + j];
    		}
	}
}

// Multiplicación de matrices en bloques
void multiplicar_matrices_bloques(double *a, double *bt, double *c, int n, int bs) {
	int i, j, k, actI, actJ;

	for (i = 0; i < n; i += bs) {
		actI = i * n;
    		for (j = 0; j < n; j += bs) {
			actJ = j * n;
        		for (k = 0; k < n; k += bs) {
            			multiplicar_bloques(&a[actI + k], &bt[actJ + k], &c[actI + j], n, bs);
        		}
    		}
	}
}

// Multiplicación de bloques de matrices
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