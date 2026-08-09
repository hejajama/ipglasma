// GaugeFix.cpp is part of the CYM evolution.
// Copyright (C) 2012 Bjoern Schenke.

#include "GaugeFix.h"

#include <complex>
#include <iostream>
#include <vector>

#include "Instrumentation.h"
#include "Matrix.h"

using std::cout;
using std::endl;

namespace {

// chi is a traceless Hermitian SU(3) algebra element,
// chi = q[a] t[a], with t[a] = lambda[a]/2.  Recover the eight real
// coefficients directly from the 3x3 matrix and evaluate exp(i chi) with
// the analytic SU(3) exponential used by the Wilson-line hot path.
inline void expGaugeRotationSU3(const Matrix &chi, Matrix &out) {
    double q[8];
    q[0] = 2.0 * chi.getRe(1);
    q[1] = -2.0 * chi.getIm(1);
    q[2] = chi.getRe(0) - chi.getRe(4);
    q[3] = 2.0 * chi.getRe(2);
    q[4] = -2.0 * chi.getIm(2);
    q[5] = 2.0 * chi.getRe(5);
    q[6] = -2.0 * chi.getIm(5);
    q[7] = std::sqrt(3.0) * (chi.getRe(0) + chi.getRe(4));

    complex<double> U[9];
    out.expmCoeff(q, U);

    const complex<double> I(0.0, 1.0);
    const double invSqrt3 = 1.0 / std::sqrt(3.0);

    out.set(0, 0, U[0] + 0.5 * U[3] + 0.5 * invSqrt3 * U[8]);
    out.set(1, 1, U[0] - 0.5 * U[3] + 0.5 * invSqrt3 * U[8]);
    out.set(2, 2, U[0] - invSqrt3 * U[8]);

    out.set(0, 1, 0.5 * (U[1] - I * U[2]));
    out.set(1, 0, 0.5 * (U[1] + I * U[2]));
    out.set(0, 2, 0.5 * (U[4] - I * U[5]));
    out.set(2, 0, 0.5 * (U[4] + I * U[5]));
    out.set(1, 2, 0.5 * (U[6] - I * U[7]));
    out.set(2, 1, 0.5 * (U[6] + I * U[7]));
}

}  // namespace

//**************************************************************************
// GaugeFix class.

void GaugeFix::FFTChi(
    FFT *fft, Lattice *lat, Group *group, Parameters *param, int steps) {
    const int N = param->getSize();
    int nn[2];
    nn[0] = N;
    nn[1] = N;
    int Nc = param->getNc();
    int Nc2m1 = Nc * Nc - 1;

    Matrix one(Nc, 1.);

    const int max_gfiter = steps;

    Matrix zero(Nc, 0.);
    double gresidual_prev = 10000.;
    double gresidual = 0.;
    std::vector<double> residualSite(static_cast<std::size_t>(N) * N);

    Matrix **chi;
    {
        IPG_PROFILE_SCOPE("observables.gluon_multiplicity.gauge_fix.setup");
        chi = new Matrix *[N * N];

        for (int i = 0; i < N * N; i++) {
            chi[i] = new Matrix(Nc, 0.);
        }
    }

    cout << "gauge fixing" << endl;

    for (int gfiter = 0; gfiter < max_gfiter; gfiter++) {
        gresidual = 0.;
        {
            IPG_PROFILE_SCOPE("observables.gluon_multiplicity.gauge_fix.divergence");
#pragma omp parallel
            {
                Matrix localg(Nc), localgdag(Nc), localDivA(Nc);
                Matrix localUx(Nc), localUy(Nc), localUxMx(Nc), localUyMy(Nc);

#pragma omp for collapse(2)
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        const int localpos = i * N + j;

                        // use periodic boundary conditions to have fast convergence
                        const int localposmX =
                            (i == 0) ? (N - 1) * N + j : (i - 1) * N + j;
                        const int localposmY =
                            (j == 0) ? i * N + N - 1 : i * N + j - 1;

                        localUx = lat->Ux[localpos];
                        localUy = lat->Uy[localpos];
                        localUxMx = lat->Ux[localposmX];
                        localUyMy = lat->Uy[localposmY];

                        localDivA = (localUx - localUxMx + localUy - localUyMy);

                        localg = zero;
                        for (int ig = 0; ig < Nc2m1; ig++) {
                            localg =
                                localg
                                + ((localDivA)*group->getT(ig)).trace().imag()
                                      * group->getT(ig);
                        }

                        *chi[localpos] = localgdag = localg;
                        localgdag.conjg();
                        residualSite[static_cast<std::size_t>(localpos)] =
                            ((localgdag * localg).trace()).real()
                            / static_cast<double>(Nc);
                    }
                }
            }

            // Preserve the original serial accumulation order so the
            // convergence test sees the same floating-point summation order.
            for (int localpos = 0; localpos < N * N; ++localpos) {
                gresidual += residualSite[static_cast<std::size_t>(localpos)];
            }
        }

        gresidual /= N * N;

        if (gfiter % 10 == 0) {
            cout << gfiter << " " << gresidual << endl;
            gresidual_prev = gresidual;
        }

        if (gresidual < 1e-9) {
            break;
        }

        if (gresidual > gresidual_prev && gresidual < 1e-6) {
            // make sure it is progressively converging
            // otherwise, break the loop with less accuracy
            break;
        }

        {
            IPG_PROFILE_SCOPE("observables.gluon_multiplicity.gauge_fix.fft_forward");
            fft->fftn(chi, chi, nn, 1);
        }

        {
            IPG_PROFILE_SCOPE("observables.gluon_multiplicity.gauge_fix.poisson");
#pragma omp parallel for
            for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double kx, ky, kt2;
                int localpos = i * N + j;
                kx = sin(
                    M_PI
                    * (-0.5 + static_cast<double>(i) / static_cast<double>(N)));
                ky = sin(
                    M_PI
                    * (-0.5 + static_cast<double>(j) / static_cast<double>(N)));
                kt2 = 4. * (kx * kx + ky * ky);  // lattice momentum squared
                *chi[localpos] = -1.5 * (1. / (kt2 + 1e-9)) * (*chi[localpos]);
                }
            }
        }

        {
            IPG_PROFILE_SCOPE("observables.gluon_multiplicity.gauge_fix.fft_backward");
            fft->fftn(chi, chi, nn, -1);
        }

        {
            IPG_PROFILE_SCOPE("observables.gluon_multiplicity.gauge_fix.exponentiate");
#pragma omp parallel
            {
                Matrix localg(Nc);
#pragma omp for
                for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    int localpos = i * N + j;
                    // chi is Hermitian and traceless here, so evaluate
                    // exp(i chi) directly in SU(3).  This avoids the generic
                    // Pade exponential and the subsequent reunitarization.
                    expGaugeRotationSU3(*chi[localpos], localg);

                    if (localg(2) != localg(2)) {
                        cout << "problem at " << i << " " << j
                             << " with g=" << localg << endl;
                        localg = one;
                    }

                    lat->Ux1[localpos] = (localg);
                    }
                }
            }
        }

        {
            IPG_PROFILE_SCOPE("observables.gluon_multiplicity.gauge_fix.transform");

            // Apply the complete gauge transformation directly to each field
            // element.  Ux1 stores the already-computed g(x) and is read-only
            // throughout this pass, so each lattice site can be updated by an
            // independent OpenMP iteration without neighboring threads writing
            // the same link.
#pragma omp parallel
            {
                Matrix g(Nc), gdag(Nc), gdagX(Nc), gdagY(Nc);

#pragma omp for collapse(2)
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        const int localpos = i * N + j;
                        const int localpospX =
                            (i == N - 1) ? j : (i + 1) * N + j;
                        const int localpospY =
                            (j == N - 1) ? i * N : i * N + j + 1;

                        g = lat->Ux1[localpos];
                        gdag = g;
                        gdag.conjg();

                        gdagX = lat->Ux1[localpospX];
                        gdagX.conjg();
                        gdagY = lat->Ux1[localpospY];
                        gdagY.conjg();

                        // The old serial sweep left-multiplied an outgoing link
                        // at x and right-multiplied it when visiting x+e_i.
                        // Preserve that operation order.  Across the periodic
                        // wrap the neighbor was visited first, so preserve the
                        // right-then-left order for those boundary links.
                        if (i == N - 1) {
                            lat->Ux[localpos] =
                                g * (lat->Ux[localpos] * gdagX);
                        } else {
                            lat->Ux[localpos] =
                                (g * lat->Ux[localpos]) * gdagX;
                        }
                        if (j == N - 1) {
                            lat->Uy[localpos] =
                                g * (lat->Uy[localpos] * gdagY);
                        } else {
                            lat->Uy[localpos] =
                                (g * lat->Uy[localpos]) * gdagY;
                        }

                        // Site-local adjoint fields are independent once g(x)
                        // is known and retain the original left-to-right matrix
                        // multiplication order.
                        lat->U[localpos] =
                            (g * lat->U[localpos]) * gdag;
                        lat->U2[localpos] =
                            (g * lat->U2[localpos]) * gdag;
                        lat->Uy2[localpos] =
                            (g * lat->Uy2[localpos]) * gdag;
                        lat->Ux2[localpos] =
                            (g * lat->Ux2[localpos]) * gdag;
                    }
                }
            }
        }
    }  // gfiter loop

    {
        IPG_PROFILE_SCOPE("observables.gluon_multiplicity.gauge_fix.cleanup");
        for (int i = 0; i < N * N; i++) {
            delete chi[i];
        }
        delete[] chi;
    }
}

void GaugeFix::gaugeTransform(Lattice *lat, Parameters *param, int i, int j) {
    int N = param->getSize();
    int pos, posmX, posmY;
    int Nc = param->getNc();
    Matrix g(Nc), gdag(Nc);

    pos = i * N + j;

    // use periodic boundary conditions to have fast convergence
    if (i == 0) {
        posmX = (N - 1) * N + j;
    } else {
        posmX = (i - 1) * N + j;
    }
    if (j == 0) {
        posmY = i * N + N - 1;
    } else {
        posmY = i * N + j - 1;
    }

    g = gdag = lat->Ux1[pos];
    gdag.conjg();

    // gauge transform Ux and Uy
    lat->Ux[pos] = (g * lat->Ux[pos]);
    lat->Uy[pos] = (g * lat->Uy[pos]);

    lat->Ux[posmX] = (lat->Ux[posmX] * gdag);
    lat->Uy[posmY] = (lat->Uy[posmY] * gdag);

    // gauge transform Ex and Ey
    lat->U[pos] = (g * lat->U[pos] * gdag);
    lat->U2[pos] = (g * lat->U2[pos] * gdag);

    // gauge transform phi and pi
    lat->Uy2[pos] = (g * lat->Uy2[pos] * gdag);
    lat->Ux2[pos] = (g * lat->Ux2[pos] * gdag);
}
