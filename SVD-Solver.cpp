#include <iostream>
#include <vector>
#include <chrono>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <cmath>
#pragma once
using namespace std;
using namespace chrono;
const double EPS = 1e-10;
const int MAX_ITER = 100;
const double SVD_THRESHOLD = 1e-10;
const double SVD_EPS = 1e-10;

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;
Matrix createMatrix(int N);
Matrix transpose(const Matrix& A);
Matrix multiply(const Matrix& A, const Matrix& B);
Vector createRightHandSide(const Matrix& A);
double computeError(const Vector& x, const Vector& x_exact);
double computeConditionNumber(const Matrix& A);
void computeEigenvalues(const Matrix& A, Vector& eigenvalues, Matrix& eigenvectors);
void svdDecomposition(const Matrix& A, Matrix& U, Vector& S, Matrix& V);
Vector solveSVD(const Matrix& U, const Vector& S, const Matrix& V, const Vector& f);

Matrix createMatrix(int N) {
    Matrix A(N, Vector(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = 1.0 / (2.0 + (i + 1) + 3.0 * (j + 1));
        }
    }
    return A;
}

Vector createRightHandSide(const Matrix& A) {
    int N = A.size();
    Vector f(N, 0.0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            f[i] += A[i][j];
        }
    }
    return f;
}

Matrix transpose(const Matrix& A) {
    int m = A.size(), n = A[0].size();
    Matrix At(n, Vector(m));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            At[j][i] = A[i][j];
    return At;
}

Matrix multiply(const Matrix& A, const Matrix& B) {
    int m = A.size(), n = B[0].size(), p = B.size();
    Matrix C(m, Vector(n, 0));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < p; ++k)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

double computeError(const Vector& x, const Vector& x_exact) {
    double diff_norm = 0.0;
    double exact_norm = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        diff_norm += (x[i] - x_exact[i]) * (x[i] - x_exact[i]);
        exact_norm += x_exact[i] * x_exact[i];
    }
    return std::sqrt(diff_norm) / std::sqrt(exact_norm);
}


double computeConditionNumber(const Matrix& A) {
    Matrix U, V;
    Vector S;
    svdDecomposition(A, U, S, V);

    if (S.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    double sigma_max = S.front();

    double sigma_min = 0.0;
    const double epsilon = 1e-15;

    for (auto it = S.rbegin(); it != S.rend(); ++it) {
        if (std::abs(*it) > epsilon) {
            sigma_min = *it;
            break;
        }
    }

    if (sigma_min == 0.0) {
        return std::numeric_limits<double>::infinity();
    }

    return sigma_max / sigma_min;
}

void householderQR(const Matrix& A, Matrix& Q, Matrix& R) {
    int n = A.size();
    Q = Matrix(n, Vector(n, 0.0));
    for (int i = 0; i < n; ++i) Q[i][i] = 1.0;
    R = A;

    for (int k = 0; k < n - 1; ++k) {
        double norm = 0.0;
        for (int i = k; i < n; ++i)
            norm += R[i][k] * R[i][k];
        norm = sqrt(norm);

        if (fabs(norm) < 1e-15) continue;

        double alpha = -copysign(norm, R[k][k]);
        Vector v(n, 0.0);
        for (int i = k; i < n; ++i)
            if (i == k) {
                v[i] = R[i][k] - alpha;
            }
            else {
                v[i] = R[i][k];
            }

        double beta = 0.0;
        for (int i = k; i < n; ++i)
            beta += v[i] * v[i];
        beta = 2.0 / beta;

        // R update
        for (int j = k; j < n; ++j) {
            double dot = 0.0;
            for (int i = k; i < n; ++i)
                dot += v[i] * R[i][j];
            for (int i = k; i < n; ++i)
                R[i][j] -= beta * v[i] * dot;
        }

        // Q update
        for (int j = 0; j < n; ++j) {
            double dot = 0.0;
            for (int i = k; i < n; ++i)
                dot += Q[j][i] * v[i];
            for (int i = k; i < n; ++i)
                Q[j][i] -= beta * v[i] * dot;
        }
    }
}


Vector solveQR(const Matrix& Q, const Matrix& R, const Vector& f) {
    int N = Q.size();

    // Compute Q^T * f  (Q is orthogonal, so Q^T = Q^H = Q.transpose())
    Vector y(N, 0.0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            y[i] += Q[j][i] * f[j];
        }
    }

    // Back substitution for Rx = y
    Vector x(N, 0.0);
    for (int i = N - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < N; ++j) {
            x[i] -= R[i][j] * x[j];
        }
        x[i] /= R[i][i];
    }

    return x;
}

void computeEigenvalues(const Matrix& A, Vector& eigenvalues, Matrix& eigenvectors) {
    int n = A.size();
    Matrix Ak = A;
    Matrix Q, R;
    eigenvectors = Matrix(n, Vector(n, 0.0));

    for (int i = 0; i < n; ++i)
        eigenvectors[i][i] = 1.0;

    for (int iter = 0; iter < MAX_ITER; ++iter) {

        householderQR(Ak, Q, R);

        Ak = multiply(R, Q);

        Matrix temp = multiply(eigenvectors, Q);
        eigenvectors = temp;
    }

    eigenvalues.resize(n);
    for (int i = 0; i < n; ++i) {
        eigenvalues[i] = Ak[i][i];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (eigenvalues[i] < eigenvalues[j]) {
                std::swap(eigenvalues[i], eigenvalues[j]);
                for (int k = 0; k < n; ++k)
                    std::swap(eigenvectors[k][i], eigenvectors[k][j]);
            }
        }
    }
}

void svdDecomposition(const Matrix& A, Matrix& U, Vector& S, Matrix& Vt) {
    int m = A.size();
    if (m == 0) return;
    int n = A[0].size();
    int k = std::min(m, n);

    // 1. Calculating AtA
    Matrix At = transpose(A);
    Matrix AtA = multiply(At, A);

    // 2. Calculating AtA's eigen values using QR decomposition
    Vector eigenvalues;
    Matrix V;
    computeEigenvalues(AtA, eigenvalues, V);

    // 3. Singular number - root of eigen number
    S.resize(k);
    for (int i = 0; i < k; ++i) {
        S[i] = sqrt(fabs(eigenvalues[i]));
    }

    // 4. Calculating Vt
    Vt = transpose(V);

    // 5. Calculating U like A*V*diag(S)^(-1)
    U = Matrix(m, Vector(k, 0.0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < k; ++j) {
            if (S[j] > SVD_EPS) {
                for (int l = 0; l < n; ++l) {
                    U[i][j] += A[i][l] * Vt[j][l] / S[j];
                }
            }
        }
    }

    // 6. New U's columns normalization
    for (int j = 0; j < k; ++j) {

        // first normalization
        double norm = 0.0;
        for (int i = 0; i < m; ++i) {
            norm += U[i][j] * U[i][j];
        }
        norm = sqrt(norm);
        if (norm > SVD_EPS) {
            for (int i = 0; i < m; ++i) {
                U[i][j] /= norm;
            }
        }

        // second normalization
        for (int p = 0; p < j; ++p) {
            double dot = 0.0;
            for (int i = 0; i < m; ++i) {
                dot += U[i][p] * U[i][j];
            }
            for (int i = 0; i < m; ++i) {
                U[i][j] -= dot * U[i][p];
            }
        }
    }
}


Vector solveSVD(const Matrix& U, const Vector& S, const Matrix& Vt, const Vector& f) {
    int m = U.size();
    int n = Vt[0].size();

    // 1. Ut * f
    Vector y(m, 0.0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            y[i] += U[j][i] * f[j];
        }
    }

    // 2. diag(S)^(-1) * (Ut * f)
    double max_s = *std::max_element(S.begin(), S.end());
    double threshold = max_s * std::max(m, n) * SVD_THRESHOLD;
    for (int i = 0; i < S.size(); ++i) {
        if (S[i] > threshold) {
            y[i] /= S[i];
        }
        else {
            y[i] = 0.0;
        }
    }

    // 3. V * (diag(S)^(-1) * Ut * f)
    Vector x(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            x[i] += Vt[j][i] * (j < y.size() ? y[j] : 0.0);
        }
    }

    return x;
}

int main() {
    std::vector<int> sizes = { 5, 10, 20 };
    for (int N : sizes) {
        Matrix A = createMatrix(N);
        Vector f = createRightHandSide(A);
        Vector x_exact(N, 1.0);
        double cond = computeConditionNumber(A);
        double final_error = 0.0;
        auto start = steady_clock::now();
        Matrix U, V;
        Vector S;
        svdDecomposition(A, U, S, V);
        Vector x = solveSVD(U, S, V, f);
        auto stop = steady_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        for (int i = 0; i < N; i++) {
            cout << S[i] << " ";
        }
        cout << endl;
        final_error = computeError(x, x_exact);
        std::cout << "Time (s): " << duration.count() << endl << "Error: " << final_error * 1e-3 << endl << "Conditional: " << cond << endl;
        cout << "___________" << endl;
    }

    return 0;
}
