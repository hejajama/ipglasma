#ifndef Matrix_h
#define Matrix_h

#include <complex>
#include <string>
#include <vector>

#include "Spinor.h"

using std::ostream;

class Matrix {
  private:
    // Small-buffer optimization: matrices up to 3x3 (SU(2), SU(3)) live
    // entirely inline -- no heap allocation per Matrix. Larger dimensions
    // (not used anywhere in IP-Glasma at present) fall back to heap storage.
    static constexpr int kInlineCap = 9;
    int ndim = 0;
    int nn = 0;
    complex<double> sbuf[kInlineCap];  // zeroed explicitly where needed
    std::vector<complex<double>> hbuf;  // only used when nn > kInlineCap
    complex<double> *e = sbuf;

    void setupStorage() {
        if (nn <= kInlineCap) {
            e = sbuf;
        } else {
            hbuf.assign(nn, complex<double>(0.0, 0.0));
            e = hbuf.data();
        }
    }

  public:
    // Tag type for constructing a Matrix without zero-initializing its
    // entries. Used internally by operators whose result has every element
    // written before use -- profiling showed ~15% of total runtime was
    // spent zero-filling temporaries that were immediately overwritten.
    struct NoInitTag {};
    static constexpr NoInitTag noInit{};

    // constructor(s)
    Matrix() : ndim(0), nn(0) { e = sbuf; }
    Matrix(int n);
    Matrix(int n, double a);
    Matrix(int n, NoInitTag) : ndim(n), nn(n * n) { setupStorage(); }

    // raw access to the element array (row-major, e[j + ndim*i])
    complex<double> *data() { return e; }
    const complex<double> *data() const { return e; }

    Matrix(const Matrix &o) : ndim(o.ndim), nn(o.nn) {
        setupStorage();
        for (int i = 0; i < nn; i++) e[i] = o.e[i];
    }
    Matrix &operator=(const Matrix &o) {
        if (this != &o) {
            if (nn != o.nn) {
                ndim = o.ndim;
                nn = o.nn;
                setupStorage();
            } else {
                ndim = o.ndim;
            }
            for (int i = 0; i < nn; i++) e[i] = o.e[i];
        }
        return *this;
    }

    // destructor
    ~Matrix() {}

    Matrix &inv();
    Matrix &logm_pade(const int m);
    Matrix &sqrtm(const int scale = 1);
    double OneNorm();
    double FrobeniusNorm();

    Matrix &logm();

    void setRe(int i, double a) { e[i] = complex<double>(a, e[i].imag()); };
    void setRe(int i, int j, double a) {
        e[j + ndim * i] = complex<double>(a, e[j + ndim * i].imag());
    };
    void setIm(int i, double a) { e[i] = complex<double>(e[i].real(), a); };
    void setIm(int i, int j, double a) {
        e[j + ndim * i] = complex<double>(e[j + ndim * i].real(), a);
    };

    void set(int i, complex<double> a) { e[i] = a; };
    void set(int i, int j, complex<double> a) { e[j + ndim * i] = a; };

    complex<double> get(int i) { return e[i]; };
    complex<double> get(int i, int j) { return e[j + ndim * i]; };

    double getRe(int i) { return e[i].real(); };
    double getIm(int i) { return e[i].imag(); };

    int getNDim() const { return ndim; }
    int getNN() const { return nn; }

    std::string MatrixToString();

    Matrix &expm(double t = 1.0, const int p = 6);

    // Matrix exponential of traceless hermitian
    // matrix using coefficients of t^a as input
    std::vector<complex<double>> expmCoeff(std::vector<double> &Q, int n);

    complex<double> det();
    complex<double> trace();

    void reu() {
        if (ndim == 3) {
            Spinor e1(ndim);
            Spinor e2(ndim);
            Spinor e3(ndim);

            Spinor a1(ndim, e[0], e[1], e[2]);
            Spinor a2(ndim, e[3], e[4], e[5]);

            e1 = a1.normalize();
            e2 = a2.GramSchmidt(e1);
            e3 = (e1 % e2).normalize();

            e[0] = e1(0);
            e[1] = e1(1);
            e[2] = e1(2);
            e[3] = e2(0);
            e[4] = e2(1);
            e[5] = e2(2);
            e[6] = e3(0);
            e[7] = e3(1);
            e[8] = e3(2);
        } else if (ndim == 2) {
            Spinor e1(ndim);
            Spinor e2(ndim);
            Spinor a1(ndim, e[0], e[1]);
            Spinor a2(ndim, e[2], e[3]);

            // cout << "a1=" << a1 << endl << endl;
            // cout << "a2=" << a2 << endl << endl;

            e1 = a1.normalize();
            e2 = a2.GramSchmidt(e1);

            // cout << "e1=" << e1 << endl << endl;
            // cout << "e2=" << e2 << endl << endl;

            e[0] = e1(0);
            e[1] = e1(1);
            e[2] = e2(0);
            e[3] = e2(1);
        }
    }

    void reu2();

    // operators:

    //()
    std::complex<double> operator()(const int i) const { return e[i]; }
    std::complex<double> operator()(const int i, const int j) const {
        return e[j + ndim * i];
    }
    // std::complex<double> operator () (const int i, const int j) { return
    // e[j+ndim*i]; }

    //=
    // const Matrix& operator = (const Matrix& p)  {
    //    nn = p.getNN();
    //    if (&p != this ) {
    //        for (int i=0; i<nn; i++) {
    //            e[i] = p.e[i];
    //        }
    //    }
    //    return *this;
    //}

    //==
    bool operator==(const Matrix &p) const {
        for (int i = 0; i < nn; i++)
            if (e[i] != p.e[i]) return false;
        return true;
    }

    //!=
    bool operator!=(const Matrix &p) const {
        for (int i = 0; i < nn; i++)
            if (e[i] != p.e[i]) return true;
        return false;
    }

    //+=
    Matrix &operator+=(const Matrix &a) {
        for (int i = 0; i < nn; i++) e[i] += a.e[i];
        return *this;
    }

    //-=
    Matrix &operator-=(const Matrix &a) {
        for (int i = 0; i < nn; i++) e[i] -= a.e[i];
        return *this;
    }

    //*=
    Matrix &operator*=(const complex<double> a) {
        for (int i = 0; i < nn; i++) e[i] *= a;
        return *this;
    }

    // /=
    Matrix &operator/=(const complex<double> a) {
        for (int i = 0; i < nn; i++) e[i] /= a;
        return *this;
    }

    double square() const {
        double tr = 0.0;
        for (int i = 0; i < nn; i++) {
            tr += e[i].real() * e[i].real() + e[i].imag() * e[i].imag();
        }
        return 0.5 * tr;
    }

    Matrix &imag();

    Matrix &conjg();
    Matrix prodABconj(const Matrix &a, const Matrix &b);
    Matrix prodAconjB(const Matrix &a, const Matrix &b);

    complex<double> traceOfProdcutOfMatrix(Matrix &M1, Matrix &M2) const;

    //<<
    friend ostream &operator<<(ostream &os, const Matrix &p) {
        for (int i = 0; i < p.getNDim(); i++) {
            for (int j = 0; j < p.getNDim(); j++) os << p(i, j);
            if (i < p.getNDim() - 1) os << std::endl;
        }
        return os;
    }
};

Matrix operator+(const Matrix &a, const Matrix &b);
Matrix operator-(const Matrix &a, const Matrix &b);
Matrix operator-(const Matrix &a);
Matrix operator/(const Matrix &a, const Matrix &b);

Matrix operator*(const double a, const Matrix &b);
Matrix operator*(const std::complex<double> a, const Matrix &b);
Matrix operator*(const Matrix &a, const double b);
Matrix operator*(const Matrix &a, const Matrix &b);

#endif
