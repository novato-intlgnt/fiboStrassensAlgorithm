#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define STEP 10
#define MAX_SIZE 1000
#define REPS 40.0
#define MIN_SIZE 5
#define ROWS 2
#define COLS 2

static void multOrdMatrix(const int *matrixA, const int *matrixB, int *matrixC,
                          int n) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; i++) {
      matrixC[i * n + j] = 0;
      for (int k = 0; k < n; k++) {
        matrixC[i * n + j] += matrixA[i * n + j] * matrixB[j * n + k];
      }
    }
  }
}

static void multDCMatrix(const int *matrixA, const int *matrixB, int *matrixC,
                         int n) {
  if (n == 1) {
    matrixC[0] = matrixA[0] * matrixB[0];
    return;
  }

  int size = n / 2;
  int bytes = size * size * sizeof(int);

  int *matA_A = malloc(bytes), *matA_B = malloc(bytes), *matA_C = malloc(bytes),
      *matA_D = malloc(bytes);
  int *matB_E = malloc(bytes), *matB_F = malloc(bytes), *matB_G = malloc(bytes),
      *matB_H = malloc(bytes);

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
}

static void merge(int arr[], int left, int mid, int right) {
  int leftSize = mid - left + 1;
  int rightSize = right - mid;

  int *leftArr = (int *)malloc(leftSize * sizeof(int));
  int *rightArr = (int *)malloc(rightSize * sizeof(int));

  if (leftArr == NULL || rightArr == NULL) {
    fprintf(stderr, "Error: no se pudo reservar memoria en merge()\n");
    free(leftArr);
    free(rightArr);
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < leftSize; i++) {
    leftArr[i] = arr[left + i];
  }

  for (int j = 0; j < rightSize; j++) {
    rightArr[j] = arr[mid + 1 + j];
  }

  int i = 0;
  int j = 0;
  int k = left;

  while (i < leftSize && j < rightSize) {
    if (leftArr[i] <= rightArr[j]) {
      arr[k++] = leftArr[i++];
    } else {
      arr[k++] = rightArr[j++];
    }
  }

  while (i < leftSize) {
    arr[k++] = leftArr[i++];
  }

  while (j < rightSize) {
    arr[k++] = rightArr[j++];
  }

  free(leftArr);
  free(rightArr);
}

static void mergeSortRecursive(int arr[], int left, int right) {
  if (left >= right) {
    return;
  }

  int mid = left + (right - left) / 2;

  mergeSortRecursive(arr, left, mid);
  mergeSortRecursive(arr, mid + 1, right);
  merge(arr, left, mid, right);
}

static void mergeSort(int arr[], int n) {
  if (n > 1) {
    mergeSortRecursive(arr, 0, n - 1);
  }
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
  multOrdMatrix(matrixA, matrixB, matrixC, n);
  clock_gettime(CLOCK_MONOTONIC, &end);

  return (double)elapsedNanoseconds(start, end) / 1000.0;
}

static double measureMergeSort(int baseArr[], int n) {
  struct timespec start;
  struct timespec end;

  clock_gettime(CLOCK_MONOTONIC, &start);
  mergeSort(baseArr, n);
  clock_gettime(CLOCK_MONOTONIC, &end);

  return (double)elapsedNanoseconds(start, end) / 1000.0;
}

int main(void) {
  srand((unsigned int)time(NULL));

  printf("# n insertion_us insertion_error merge_us merge_error\n");

  size_t totalSize = ROWS * COLS * sizeof(int);
  int *baseMatrix = (int *)malloc(totalSize);

  for (int n = 2; n <= REPS; n += STEP) {
    int *matrixA = (int *)malloc(totalSize);
    int *matrixB = (int *)malloc(totalSize);
    int *matrixC = (int *)malloc(totalSize);
    int *matrixD = (int *)malloc(totalSize);

    if (matrixA == NULL || matrixB == NULL || matrixC == NULL) {
      fprintf(stderr, "Error: no se pudo reservar memoria para n=%d\n", n);
      return EXIT_FAILURE;
    }

    genRandomSquareMatrix(baseMatrix, n, 100);
    memcpy(baseMatrix, matrixA, totalSize);
    memcpy(baseMatrix, matrixB, totalSize);
    double ordMultTemp = measureOrdMulti(matrixA, matrixB, matrixC, n);
    double mergeTime = measureMergeSort(baseMatrix, n);

    printf("%d %.3f %.3f\n", n, ordMultTemp, mergeTime);
    free(baseMatrix);
    free(matrixA);
    free(matrixB);
    free(matrixC);
  }
  return EXIT_SUCCESS;
}
