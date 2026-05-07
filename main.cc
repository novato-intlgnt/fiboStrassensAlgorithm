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
#define THRESHOLD 32

struct Matrix {
  int _n;
  int *_data;

  Matrix(int size) : _n(size), _data(new int[size * size]) {}

  ~Matrix() { delete[] _data; }
  int &operator()(int i, int j) { return _data[i * _n + j]; }

  int *getData() const { return _data; }
  const int &operator()(int i, int j) const { return _data[i * _n + j]; }
};

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

static void multOrdMatrix(int n, const Matrix &matrixA, const Matrix &matrixB,
                          Matrix &matrixC) {
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < n; k++) {
      int temp = matrixA(i, k);
      for (int j = 0; j < n; j++) {
        // C[i][j] += A[i][k] * B[k][j]
        matrixC(i, j) +=
            temp * matrixB(k, j); // remind the exchange of j and k bucles
      }
    }
  }
}

static void multiply_iter(int n, const Matrix &matrixA, const Matrix &matrixB,
                          Matrix &matrixC, int rowA, int colA, int rowB,
                          int colB, int rowC, int colC) {
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < n; k++) {
      int temp = matrixA((rowA + i), (colA + k));
      for (int j = 0; j < n; j++) {
        // C[i][j] += A[i][k] * B[k][j]
        matrixC((rowC + i), (colC + j)) +=
            temp * matrixB((rowB + k), (colB + j));
      }
    }
  }
}

// Without offsets, is simple but reserve more memory
static void multDCMatrix(int n, const Matrix &matrixA, int rowA, int colA,
                         const Matrix &matrixB, int rowB, int colB,
                         Matrix &matrixC, int rowC, int colC) {
  if (n <= THRESHOLD) {
    multiply_iter(n, matrixA, matrixB, matrixC, rowA, colA, rowB, colB, rowC,
                  colC);
    return;
  }

  int size = n / 2;

  // Cr = (Aa*Be) + (Ab*Bg)
  multDCMatrix(size, matrixA, rowA, colA, matrixB, rowB, colB, matrixC, rowC,
               colC); // (Aa*Be)
  multDCMatrix(size, matrixA, rowA, colA + size, matrixB, rowB + size, colB,
               matrixC, rowC,
               colC); // (Ab*Bg)

  // Cs = (Aa*Bf) + (Ab*Bh)
  multDCMatrix(size, matrixA, rowA, colA, matrixB, rowB, colB + size, matrixC,
               rowC,
               colC + size); // (Aa*Bf)
  multDCMatrix(size, matrixA, rowA, colA + size, matrixB, rowB + size,
               colB + size, matrixC, rowC,
               colC + size); // (Ab*Bh)

  // Ct = (Ac*Be) + (Ad*Bg)
  multDCMatrix(size, matrixA, rowA + size, colA, matrixB, rowB, colB, matrixC,
               rowC + size,
               colC); // (Ac*Be)
  multDCMatrix(size, matrixA, rowA + size, colA + size, matrixB, rowB + size,
               colB, matrixC, rowC + size,
               colC); // (Ad*Bg)

  // Cu = (Ac*Bf) + (Ad*Bh)
  multDCMatrix(size, matrixA, rowA + size, colA, matrixB, rowB, colB + size,
               matrixC, rowC + size,
               colC + size); // (Ac*Bf)
  multDCMatrix(size, matrixA, rowA + size, colA + size, matrixB, rowB + size,
               colB + size, matrixC, rowC + size,
               colC + size); // (Ad*Bh)
}

static void genRandomSquareMatrix(Matrix &matrix, int n, int maxValue) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      matrix(i, j) = rand() % (maxValue + 1);
    }
  }
}

static long long elapsedNanoseconds(struct timespec start,
                                    struct timespec end) {
  long long seconds = (long long)(end.tv_sec - start.tv_sec);
  long long nanoseconds = (long long)(end.tv_nsec - start.tv_nsec);
  return seconds * 1000000000LL + nanoseconds;
}

static double measureOrdMulti(int n, const Matrix &matrixA,
                              const Matrix &matrixB, Matrix &matrixC) {
  struct timespec start;
  struct timespec end;

  clock_gettime(CLOCK_MONOTONIC, &start);
  multOrdMatrix(n, matrixA, matrixB, matrixC);
  clock_gettime(CLOCK_MONOTONIC, &end);

  return (double)elapsedNanoseconds(start, end) / 1000.0;
}

static double measuremultDCMulti(int n, const Matrix &matrixA, int rowA,
                                 int colA, const Matrix &matrixB, int rowB,
                                 int colB, Matrix &matrixC, int rowC,
                                 int colC) {
  struct timespec start;
  struct timespec end;

  clock_gettime(CLOCK_MONOTONIC, &start);
  multDCMatrix(n, matrixA, rowA, colA, matrixB, rowB, colB, matrixC, rowC,
               colC);
  clock_gettime(CLOCK_MONOTONIC, &end);

  return (double)elapsedNanoseconds(start, end) / 1000.0;
}

int main(void) {
  srand((unsigned int)time(NULL));

  printf("# n ord_us dc_us\n");

  for (int n = 1; n <= REPS; n++) {
    int size = myPow(2, n);
    int totalSize = size * size;

    Matrix baseMatrix(size), matrixA(size), matrixB(size), matrixC(size),
        matrixD(size);

    genRandomSquareMatrix(baseMatrix, size, 100);
    memcpy(matrixA.getData(), baseMatrix.getData(), totalSize * sizeof(int));
    memcpy(matrixB.getData(), baseMatrix.getData(), totalSize * sizeof(int));
    memset(matrixC.getData(), 0, totalSize * sizeof(int));
    memset(matrixD.getData(), 0, totalSize * sizeof(int));

    double ordMultTemp = measureOrdMulti(size, matrixA, matrixB, matrixC);
    double multDCTemp =
        measuremultDCMulti(size, matrixA, 0, 0, matrixB, 0, 0, matrixD, 0, 0);

    printf("%d %.3f %.3f\n", size, ordMultTemp, multDCTemp);
  }
  return EXIT_SUCCESS;
}
