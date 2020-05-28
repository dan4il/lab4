#include <iostream>
#include <limits>
#include <cmath>

using namespace std;

const int SIZE = 100;

double* gauss(double a[SIZE][SIZE], double f[SIZE]);
double* pvr(const double a[SIZE][SIZE], const double f[SIZE]);
double* nev(double a[SIZE][SIZE], double f[SIZE], double* x);
double det(const double a[SIZE][SIZE], int n = SIZE, int iexcl = -1, int jexc = -1);
double norm(const double v[SIZE]);
double norm(const double a[SIZE][SIZE]);
void copy(const double* source, double* dest, int n = SIZE);
void inverseMatrix(const double a[SIZE][SIZE], double aInv[SIZE][SIZE]);
bool converge(double xk[SIZE], double xkp[SIZE], int SIZE);

int main()
{
    double a[SIZE][SIZE];   // матрица коэффициентов
    double f[SIZE];         //  свободные эл-ты

    // инициализация 
    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < SIZE; j++)
        {
            if (i == j) 
                a[i][j] = 10;
            else if (j <= i + 1) 
                a[i][j] = 1. / (i + 1);
            else 
                a[i][j] = 0;
        }
        f[i] = i + 1;
    }


    double a0[SIZE][SIZE];
    double f0[SIZE];
    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < SIZE; j++)
            a0[i][j] = a[i][j];
    }
    copy(f, f0);

    double* x1 = gauss(a, f);
    double* x2 = pvr(a, f);


    double* gaussNev = nev(a0, f0, x1);
    double* pvrNev = nev(a0, f0, x2);

    cout << "Norma nev Gauss:" << norm(gaussNev) << endl;
    cout << "Norma nev PVR:" << norm(pvrNev) << endl;

    double inv[SIZE][SIZE];
    inverseMatrix(a0, inv);

    cout << "Chislo obuslov: " << norm(a0) * norm(inv) << endl;

    delete[] x1;
    delete[] x2;
    delete[] gaussNev;
    delete[] pvrNev;
    return 0;
}

// Решение Гауссом
double* gauss(double a[SIZE][SIZE], double f[SIZE])
    {
    double* x = new double[SIZE];
    for (int i = 0; i < SIZE; i++)
    {
        for (int j = i + 1; j < SIZE; j++)
        {
            double c = a[j][i] / a[i][i];
            for (int k = i; k < SIZE; k++) 
                a[j][k] -= a[i][k] * c;
            f[j] -= f[i] * c;
        }
    }
    for (int i = SIZE - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int j = SIZE - 1; j > i; j--)
            sum += a[i][j] * x[j];
        x[i] = (f[i] - sum) / a[i][i];
    }
    return x;
}

double* pvr(const double a[SIZE][SIZE], const double f[SIZE]) 
{
    int iter = 0;
    double* p = new double[SIZE];
    double* x = new double[SIZE];
    for (int i = 0; i < SIZE; i++)
        x[i] = 1;
    do
    {
        for (int i = 0; i < SIZE; i++)
            p[i] = x[i];
        for (int i = 0; i < SIZE; i++)
        {
            double var = 0;
            for (int j = 0; j < i; j++)
                var += (a[i][j] * x[j]);
            for (int j = i + 1; j < SIZE; j++)
                var += (a[i][j] * p[j]);
            x[i] = (f[i] - var) / a[i][i];
        }
        iter++;
    } while (!converge(x, p, SIZE));

    cout<< "Iterations: " << iter << endl;
    return x;
 }

double* nev(double a[SIZE][SIZE], double f[SIZE], double* x)
{
    double* vnev = new double[SIZE];
    for (int i = 0; i < SIZE; i++) 
    {
        double sum = 0;
        for (int j = 0; j < SIZE; j++) 
            sum += a[i][j] * x[j];
        vnev[i] = f[i] - sum;
    }
    return vnev;
}
void copy(const double* source, double* dest, int n)
{
    for (int i = 0; i < n; i++)
        dest[i] = source[i];
}
double norm(const double v[SIZE])
{
    double sum = 0;
    for (int i = 0; i < SIZE; i++)
        sum += v[i] * v[i];
    return sqrt(sum);
}

double norm(const double a[SIZE][SIZE]) 
{
    double sum;
    double max = -1 * INFINITY;
    for (int i = 0; i < SIZE; i++) 
    {
        sum = 0;
        for (int j = 0; j < SIZE; j++) 
            sum += a[i][j];
        if (sum > max) max = sum;
    }
    return max;
}

double det(const double a[SIZE][SIZE], int n, int iexcl, int jexc)
{
    double an[SIZE][SIZE];
    for (int i = 0, i0 = 0; i < n; i++, i0++) 
    {
        if (i0 == iexcl) i0++;
        for (int j = 0, j0 = 0; j < n; j++, j0++) 
        {
            if (j0 == jexc) j0++;
            an[i][j] = a[i0][j0];
        }
    }
    for (int i = 0; i < n; i++) 
    {
        for (int j = i + 1; j < n; j++) 
        {
            double c = an[j][i] / an[i][i];
            for (int k = i; k < n; k++)
                an[j][k] -= an[i][k] * c;
        }
    }
    double deter = 1;
    for (int i = 0; i < n; i++) 
        deter *= an[i][i];
    return deter;
}

void inverseMatrix(const double a[SIZE][SIZE], double inv[SIZE][SIZE])
{
    double deter = det(a);
    for (int i = 0; i < SIZE; i++) 
    {
        for (int j = 0; j < SIZE; j++) 
            inv[j][i] = det(a, SIZE - 1, i, j) * pow(-1, i + j) / deter;
    }
}

bool converge(double xk[SIZE], double xkp[SIZE], int SIZE )
{
    double norm = 0;
    for (int i = 0; i < SIZE; i++)
        norm += (xk[i] - xkp[i]) * (xk[i] - xkp[i]);
    return (sqrt(norm) < numeric_limits<double>::epsilon());
}

