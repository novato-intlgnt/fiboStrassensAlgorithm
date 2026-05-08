#define _POSIX_C_SOURCE 200809L
#include <algorithm>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define STEP 1
#define MAX_SIZE 12
#define REPS 10
#define MIN_SIZE 5
#define THRESHOLD 64

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

void addMat(int n, const Matrix &matrixA, int rowA, int colA,
            const Matrix &matrixB, int rowB, int colB, Matrix &dest) {
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      dest(i, j) = matrixA(rowA + i, colA + j) + matrixB(rowB + i, colB + j);
}
void subMat(int n, const Matrix &matrixA, int rowA, int colA,
            const Matrix &matrixB, int rowB, int colB, Matrix &dest) {
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      dest(i, j) = matrixA(rowA + i, colA + j) - matrixB(rowB + i, colB + j);
}

static void multStrassenMat(int n, const Matrix &matrixA, int rowA, int colA,
                            const Matrix &matrixB, int rowB, int colB,
                            Matrix &matrixC, int rowC, int colC) {

  if (n <= THRESHOLD) {
    multiply_iter(n, matrixA, matrixB, matrixC, rowA, colA, rowB, colB, rowC,
                  colC);
    return;
  }
  int size = n / 2;
  Matrix P1(size), P2(size), P3(size), P4(size), P5(size), P6(size), P7(size),
      auxA(size), auxB(size);

  // P1 = Aa * (Bf - Bh)
  subMat(size, matrixB, rowB, colB + size, matrixB, rowB + size, colB + size,
         auxB);
  multStrassenMat(size, matrixA, rowA, colA, auxB, 0, 0, P1, 0, 0);

  // P2 = (Aa + Ab) * Bh
  addMat(size, matrixA, rowA, colA, matrixA, rowA, colA + size, auxA);
  multStrassenMat(size, auxA, 0, 0, matrixB, rowB + size, colB + size, P2, 0,
                  0);

  // P3 = (Ac + Ad) * Be
  addMat(size, matrixA, rowA + size, colA, matrixA, rowA + size, colA + size,
         auxA);
  multStrassenMat(size, auxA, 0, 0, matrixB, rowB, colB, P3, 0, 0);

  // P4 = Ad * (Bg - Be)
  subMat(size, matrixB, rowB + size, colB, matrixB, rowB, colB, auxB);
  multStrassenMat(size, matrixA, rowA + size, colA + size, auxB, 0, 0, P4, 0,
                  0);

  // P5 = (Aa + Ad) * (Be + Bh)
  addMat(size, matrixA, rowA, colA, matrixA, rowA + size, colA + size, auxA);
  addMat(size, matrixB, rowB + size, colB, matrixB, rowB, colB, auxB);
  multStrassenMat(size, auxA, 0, 0, auxB, 0, 0, P5, 0, 0);

  // P6 = (Ab - Ad) * (Bg + Bh)
  subMat(size, matrixA, rowA, colA + size, matrixA, rowA + size, colA + size,
         auxA);
  addMat(size, matrixB, rowB + size, colB, matrixB, rowB + size, colB + size,
         auxB);
  multStrassenMat(size, auxA, 0, 0, auxB, 0, 0, P6, 0, 0);

  // P7 = (Aa - Ac) * (Be + Bf)
  subMat(size, matrixA, rowA, colA, matrixA, rowA + size, colA, auxA);
  addMat(size, matrixB, rowB, colB, matrixB, rowB, colB + size, auxB);
  multStrassenMat(size, auxA, 0, 0, auxB, 0, 0, P7, 0, 0);

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      // r = p5 + p4 - p2 + p6
      matrixC(rowC + i, colC + j) = P5(i, j) + P4(i, j) - P2(i, j) + P6(i, j);

      // s = p1 + p2
      matrixC(rowC + i, colC + j + size) = P1(i, j) + P2(i, j);

      // t = p3 + p4
      matrixC(rowC + i + size, colC + j) = P3(i, j) + P4(i, j);

      // u = p5 + p1 - p3 - p7
      matrixC(rowC + i + size, colC + j + size) =
          P5(i, j) + P1(i, j) - P3(i, j) + P7(i, j);
    }
  }
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

static double measuremultStraMulti(int n, const Matrix &matrixA, int rowA,
                                   int colA, const Matrix &matrixB, int rowB,
                                   int colB, Matrix &matrixC, int rowC,
                                   int colC) {
  struct timespec start;
  struct timespec end;

  clock_gettime(CLOCK_MONOTONIC, &start);
  multStrassenMat(n, matrixA, rowA, colA, matrixB, rowB, colB, matrixC, rowC,
                  colC);
  clock_gettime(CLOCK_MONOTONIC, &end);

  return (double)elapsedNanoseconds(start, end) / 1000.0;
}

double calculateMedian(double *arr, int n) {
  std::sort(arr, arr + n);
  if (n % 2 == 0) {
    return (arr[n / 2 - 1] + arr[n / 2]) / 2.0;
  } else {
    return arr[n / 2];
  }
}
int main(void) {
  srand((unsigned int)time(NULL));

  printf("# n ord_us dc_us st_us\n");

  for (int n = 1; n <= MAX_SIZE; n++) {
    int size = myPow(2, n);
    size_t totalElements = (size_t)size * size;
    size_t totalSizeBytes = totalElements * sizeof(int);

    Matrix baseMatrix(size);
    genRandomSquareMatrix(baseMatrix, size, 100);

    double tmpORD[REPS], tmpDC[REPS], tmpSTRA[REPS];

    for (int i = 0; i < REPS; i++) {
      Matrix matrixA(size), matrixB(size), matrixC(size);

      memcpy(matrixA.getData(), baseMatrix.getData(), totalSizeBytes);
      memcpy(matrixB.getData(), baseMatrix.getData(), totalSizeBytes);

      memset(matrixC.getData(), 0, totalSizeBytes);
      tmpORD[i] = measureOrdMulti(size, matrixA, matrixB, matrixC);

      memset(matrixC.getData(), 0, totalSizeBytes);
      tmpDC[i] =
          measuremultDCMulti(size, matrixA, 0, 0, matrixB, 0, 0, matrixC, 0, 0);

      memset(matrixC.getData(), 0, totalSizeBytes);
      tmpSTRA[i] = measuremultStraMulti(size, matrixA, 0, 0, matrixB, 0, 0,
                                        matrixC, 0, 0);
    }

    double medianORD = calculateMedian(tmpORD, REPS);
    double medianDC = calculateMedian(tmpDC, REPS);
    double medianSTRA = calculateMedian(tmpSTRA, REPS);

    printf("%d %.3f %.3f %.3f\n", size, medianORD, medianDC, medianSTRA);
  }
  return EXIT_SUCCESS;
}
