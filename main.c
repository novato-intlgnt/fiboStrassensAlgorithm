#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define STEP 1
#define MAX_SIZE 2
#define REPS 12
#define MIN_SIZE 5
#define ROWS MAX_SIZE
#define COLS MAX_SIZE

int myPow(int a, int n) {
  if (n == 0) {
    return 1;
  }
  if (n == 1) {
    return a;
  }

  int res = myPow(a, n / 2);
  return (res * res * myPow(a, n % 2));
}

static void multOrdMatrix(int n, const int *matrixA, const int *matrixB,
                          int *matrixC) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        // C[i][j] += A[i][k] * B[k][j]
        matrixC[i * n + j] += matrixA[i * n + k] * matrixB[k * n + j];
      }
    }
  }
}

// Without offsets, is simple but reserve more memory
static void multDCMatrix(int n, const int *matrixA, const int *matrixB,
                         int *matrixC) {
  if (n == 2) {
    multOrdMatrix(n, matrixA, matrixB, matrixC);
    return;
  }

  int size = n / 2;
  int bytes = size * size * sizeof(int);

  int *matA_A = NULL, *matA_B = NULL, *matA_C = NULL, *matA_D = NULL;
  int *matB_E = NULL, *matB_F = NULL, *matB_G = NULL, *matB_H = NULL;
  matA_A = malloc(bytes);
  matA_B = malloc(bytes);
  matA_C = malloc(bytes);
  matA_D = malloc(bytes);
  matB_E = malloc(bytes);
  matB_F = malloc(bytes);
  matB_G = malloc(bytes);
  matB_H = malloc(bytes);

  int *p1 = calloc(size * size, sizeof(int)),
      *p2 = calloc(size * size, sizeof(int));

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      matA_A[i * size + j] = matrixA[i * n + j];
      matA_B[i * size + j] = matrixA[i * n + (j + size)];
      matA_C[i * size + j] = matrixA[(i + size) * n + j];
      matA_D[i * size + j] = matrixA[(i + size) * n + (j + size)];

      matB_E[i * size + j] = matrixB[i * n + j];
      matB_F[i * size + j] = matrixB[i * n + (j + size)];
      matB_G[i * size + j] = matrixB[(i + size) * n + j];
      matB_H[i * size + j] = matrixB[(i + size) * n + (j + size)];
    }
  }

  // Cr = (Aa*Be) + (Ab*Bg)
  multDCMatrix(size, matA_A, matB_E, p1);
  multDCMatrix(size, matA_B, matB_G, p2);
  for (int i = 0; i < size * size; i++)
    matrixC[(i / size) * n + (i % size)] = p1[i] + p2[i];

  // Cs = (Aa*Bf) + (Ab*Bh)
  for (int i = 0; i < size * size; i++) {
    p1[i] = 0;
    p2[i] = 0;
  } // Reset temporales
  multDCMatrix(size, matA_A, matB_F, p1);
  multDCMatrix(size, matA_B, matB_H, p2);
  for (int i = 0; i < size * size; i++)
    matrixC[(i / size) * n + (i % size + size)] = p1[i] + p2[i];

  // Ct = (Ac*Be) + (Ad*Bg)
  for (int i = 0; i < size * size; i++) {
    p1[i] = 0;
    p2[i] = 0;
  }
  multDCMatrix(size, matA_C, matB_E, p1);
  multDCMatrix(size, matA_D, matB_G, p2);
  for (int i = 0; i < size * size; i++)
    matrixC[((i / size) + size) * n + (i % size)] = p1[i] + p2[i];

  // Cu = (Ac*Bf) + (Ad*Bh)
  for (int i = 0; i < size * size; i++) {
    p1[i] = 0;
    p2[i] = 0;
  }
  multDCMatrix(size, matA_C, matB_F, p1);
  multDCMatrix(size, matA_D, matB_H, p2);
  for (int i = 0; i < size * size; i++)
    matrixC[((i / size) + size) * n + (i % size + size)] = p1[i] + p2[i];

  free(matA_A);
  free(matA_B);
  free(matA_C);
  free(matA_D);
  free(matB_E);
  free(matB_F);
  free(matB_G);
  free(matB_H);
  free(p1);
  free(p2);
}

static void genRandomSquareMatrix(int *matrix, int n, int maxValue) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      matrix[i * n + j] = rand() % (maxValue + 1);
    }
  }
}

static long long elapsedNanoseconds(struct timespec start,
                                    struct timespec end) {
  long long seconds = (long long)(end.tv_sec - start.tv_sec);
  long long nanoseconds = (long long)(end.tv_nsec - start.tv_nsec);
  return seconds * 1000000000LL + nanoseconds;
}

static double measureOrdMulti(const int *matrixA, const int *matrixB,
                              int *matrixC, int n) {
  struct timespec start;
  struct timespec end;

  clock_gettime(CLOCK_MONOTONIC, &start);
  multOrdMatrix(n, matrixA, matrixB, matrixC);
  clock_gettime(CLOCK_MONOTONIC, &end);

  return (double)elapsedNanoseconds(start, end) / 1000.0;
}

static double measuremultDCMulti(const int *matrixA, const int *matrixB,
                                 int *matrixC, int n) {
  struct timespec start;
  struct timespec end;

  clock_gettime(CLOCK_MONOTONIC, &start);
  multDCMatrix(n, matrixA, matrixB, matrixC);
  clock_gettime(CLOCK_MONOTONIC, &end);

  return (double)elapsedNanoseconds(start, end) / 1000.0;
}

int main(void) {
  srand((unsigned int)time(NULL));

  printf("# n ord_us dc_us\n");

  for (int n = 1; n <= REPS; n++) {
    int pow2 = myPow(2, n);
    size_t totalSize = (size_t)pow2 * pow2 * sizeof(int);

    // baseMatrix se aloca con el tamaño correcto en cada iteración
    int *baseMatrix = malloc(totalSize);
    int *matrixA = malloc(totalSize);
    int *matrixB = malloc(totalSize);
    int *matrixC = calloc((size_t)pow2 * pow2, sizeof(int));
    int *matrixD = calloc((size_t)pow2 * pow2, sizeof(int));

    if (baseMatrix == NULL || matrixA == NULL || matrixB == NULL ||
        matrixC == NULL || matrixD == NULL) {
      fprintf(stderr, "Error: no se pudo reservar memoria para n=%d\n", pow2);
      free(baseMatrix);
      free(matrixA);
      free(matrixB);
      free(matrixC);
      free(matrixD);
      return EXIT_FAILURE;
    }

    genRandomSquareMatrix(baseMatrix, pow2, 100);
    memcpy(matrixA, baseMatrix, totalSize);
    memcpy(matrixB, baseMatrix, totalSize);
    double ordMultTemp = measureOrdMulti(matrixA, matrixB, matrixC, pow2);
    double multDCTemp = measuremultDCMulti(matrixA, matrixB, matrixD, pow2);

    printf("%d %.3f %.3f\n", pow2, ordMultTemp, multDCTemp);
    free(baseMatrix);
    free(matrixA);
    free(matrixB);
    free(matrixC);
    free(matrixD);
  }
  return EXIT_SUCCESS;
}
