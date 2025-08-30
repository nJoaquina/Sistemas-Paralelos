#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <float.h>
#include <math.h>

#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1
#define BS 64

// Prototipos de funciones
void inicializar_matriz(double *mat, int n, double val, int orden);
//void multiplicar_matrices_bloques(double *a, double *bt, double *c, int n, int bs,  int chunk, int id);
void multiplicar_bloques(double *bloque_a, double *bloque_b, double *bloque_c, int n, int bs);
void transpuesta(double *mat, double *matT, int n);
double dwalltime();

int main(int argc, char *argv[]) {
    double *A, *B, *BT, *C, *R, *resultMatriz;
    int N, T;
    int i, j,k ;
    double inicioTiempo, totalTiempo;
    double MinA , MaxA, promA, totalA;
    double MinB , MaxB, promB, totalB;
    double escalar;

    // ValidaciÃ³n de parÃ¡metros
    if (argc != 3) {
        printf("Error: Se esperaban 2 argumentos (N y T).\n");
        exit(1);
    }
 
    N = strtol(argv[1], NULL, 10);
    T = strtol(argv[2], NULL, 10);


    if (N <= 0 || T <= 0 || N % BS != 0) {
        printf("Parametros mal ingresados. Ejemplo de valores esperados: 1024 4\n");
        exit(1);
    }

    // Seteo de cantidad de hilos
     omp_set_num_threads(T);

    // Reserva de memoria
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

    inicioTiempo = dwalltime();

    transpuesta(B, BT, N);
    MinA = MinB =DBL_MAX;
    MaxA = MaxB =DBL_MIN;
    totalA = totalB = 0.0;

      #pragma omp parallel private(i,j,k)
    { 
        #pragma omp for  reduction(+:totalA) reduction(max: MaxA) reduction (min: MinA) schedule (static) //static divide en cant iteraciones/cant hilos
        for (i = 0; i < N; i++) {
            double valA = A[i * N];
            if (valA < MinA) MinA = valA;
            if (valA > MaxA) MaxA = valA;
            totalA += valA;
        }//barrera implicita

        // Calcula maximo y minimo local de B
        #pragma omp for  reduction(+:totalB) reduction(max: MaxB) reduction(min: MinB) schedule(static)
        for (i = 0; i < N; i++) {
            double valB = B[i * N];
            if (valB < MinB) MinB = valB;
            if (valB > MaxB) MaxB = valB;
            totalB += valB;
        }//barrera implicita sincroniza a todos los hilos que lo integran

        // Esperamos que todos los hilos terminen para poder calcular el resultado correcto
       // #pragma omp barrier // esta barrera estaria demas por que  el #pragma de arriba tiene barrera implicita


  //     #pragma omp single { //barrera implicta la final, los demas hilos esperan a que esta parte se ejecuten
       #pragma omp single
 {
            totalA /= size;
            totalB /= size;

            if (totalA == 0 || totalB == 0) {
                printf("Error: Division por cero\n");
                exit(1);
            }

            escalar = (MaxA * MaxB - MinA * MinB) / (totalA * totalB);
      }
        // Multiplicaciones matrices AxB //cuando un hilo termina no espera sigue para multiplicar CxBT ya uqe no hay datos compartidos
        #pragma omp for nowait schedule (static)
        for (i = 0; i < N; i += BS) {
          int actI = i * N;
          for (j = 0; j < N; j += BS) {
            int actJ = j * N;
            for (k = 0; k < N; k += BS) {
                multiplicar_bloques(&A[actI + k], &B[actJ + k], &C[actI + j], N, BS);
            }
        }
    }
	// Multiplicaciones matrices CxBT
        #pragma omp for schedule (static)
        for (i = 0; i < N; i += BS) {
          int actI = i * N;
          for (j = 0; j < N; j += BS) {
            int actJ = j * N;
            for (k = 0; k < N; k += BS) {
                multiplicar_bloques(&C[actI + k], &BT[actJ + k], &R[actI + j], N, BS);
            }
        }
    }

        // R = R + resultMatriz * escalar
        #pragma omp for nowait schedule (static)
       for (int i = 0; i < N; i++) {
          int desplI = i* N; 
        for(int j=0; j< N; j++){
           R[desplI + j] += resultMatriz[desplI + j] * escalar;
        }
      }
    }
    totalTiempo = dwalltime() - inicioTiempo;

    printf("Tamaï¿½o de matriz: %d, Bloque: %d, Hilos: %d, Tiempo: %f segundos\n", N, BS, T, totalTiempo);

    free(A);
    free(B);
    free(BT);
    free(C);
    free(R);
    free(resultMatriz);

    return 0;
}

// Calculo de tiempo
double dwalltime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

// Inicializar matriz
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

// Transponer matriz
void transpuesta(double *mat, double *matT, int n) {
    int i, j;
   #pragma omp for private (i,j) 
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            matT[j * n + i] = mat[i * n + j];
        }
    }
}
/*
void multiplicar_matrices_bloques(double *a, double *bt, double *c, int n, int bs, int filas_por_bloque, int id) {
    int i, j, k, comienzo, final;

    comienzo= id * filas_por_bloque;
    final= (comienzo + filas_por_bloque < n)? (comienzo + filas_por_bloque): n;

    #pragma omp for private(i,j,k)schedule (static)
    for (i = comienzo; i < final; i += bs) {
        int actI = i * n;
        for (j = 0; j < n; j += bs) {
            int actJ = j * n;
            for (k = 0; k < n; k += bs) {
                multiplicar_bloques(&a[actI + k], &bt[actJ + k], &c[actI + j], n, BS);
            }
        }
    }
} */

void multiplicar_bloques(double *bloque_a, double *bloque_b, double *bloque_c, int n, int bs) {
    int i, j, k;
    double suma;
    for (i = 0; i < bs; i++) {
        int actI = i * n;
        for (j = 0; j < bs; j++) {
            int actJ = j * n;
            suma = 0.0;
            for (k = 0; k < bs; k++) {
                suma += bloque_a[actI + k] * bloque_b[actJ + k];
            }
              bloque_c[actI + j] += suma; 
        }
         
    }   
}
