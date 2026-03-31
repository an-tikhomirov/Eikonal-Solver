#include "local_solver_anisotropic.hpp"

#include <cmath>
#include <limits>
#include <iostream>

using namespace INMOST;

static constexpr double INF_AT = 1.0e12;
static inline bool IsFiniteAT(double at) { return at < INF_AT * 0.5; }

static inline bool Invert3x3(const double A[9], double invA[9])
{
    const double a00 = A[0], a01 = A[1], a02 = A[2];
    const double a10 = A[3], a11 = A[4], a12 = A[5];
    const double a20 = A[6], a21 = A[7], a22 = A[8];

    const double c00 =  (a11*a22 - a12*a21);
    const double c01 = -(a10*a22 - a12*a20);
    const double c02 =  (a10*a21 - a11*a20);

    const double c10 = -(a01*a22 - a02*a21);
    const double c11 =  (a00*a22 - a02*a20);
    const double c12 = -(a00*a21 - a01*a20);

    const double c20 =  (a01*a12 - a02*a11);
    const double c21 = -(a00*a12 - a02*a10);
    const double c22 =  (a00*a11 - a01*a10);

    const double det = a00*c00 + a01*c01 + a02*c02;
    if (std::abs(det) < 1e-30) return false;

    const double invdet = 1.0 / det;

    invA[0] = c00 * invdet; invA[1] = c10 * invdet; invA[2] = c20 * invdet;
    invA[3] = c01 * invdet; invA[4] = c11 * invdet; invA[5] = c21 * invdet;
    invA[6] = c02 * invdet; invA[7] = c12 * invdet; invA[8] = c22 * invdet;

    return true;
}

static inline double QuadForm(const double M[9], double dx, double dy, double dz)
{
    const double v0 = M[0]*dx + M[1]*dy + M[2]*dz;
    const double v1 = M[3]*dx + M[4]*dy + M[5]*dz;
    const double v2 = M[6]*dx + M[7]*dy + M[8]*dz;
    return dx*v0 + dy*v1 + dz*v2;
}

double local_solver_anisotropic(Node v4, Tag at_tag, Tag rt_tag, Tag sigma_tag)
{
    double best_at = std::numeric_limits<double>::max();

    const Storage::real_array c4 = v4.Coords();
    const double x4[3] = {c4[0], c4[1], c4[2]};

    const double rt4 = v4.Real(rt_tag);

    ElementArray<Cell> cells = v4.getCells();

    // Плотность дискретизации минимума:
    const int N_EDGE = 40;
    const int N_TRI  = 25;

    for (int cell_idx = 0; cell_idx < cells.size(); ++cell_idx)
    {
        Cell cell = cells[cell_idx];
        ElementArray<Node> nodes = cell.getNodes();

        if (nodes.size() != 4) continue;

        // Находим локальный индекс v4 внутри тетраэдра
        int idx4 = -1;
        for (int i = 0; i < 4; ++i)
        {
            if (nodes[i].GetHandle() == v4.GetHandle())
            {
                idx4 = i;
                break;
            }
        }
        if (idx4 == -1) continue;

        // Считываем sigma_tensor
        Storage::real_array sig_arr = cell.RealArray(sigma_tag);
        if (sig_arr.size() < 9) continue;

        double sigma[9];
        for (int k = 0; k < 9; ++k) sigma[k] = sig_arr[k];

        // Инвертируем sigma
        double sigma_inv[9];
        if (!Invert3x3(sigma, sigma_inv))
            continue;

        // Остальные 3 вершины тетраэдра -> грань, противоположная v4
        Node v[3];
        int kk = 0;
        for (int i = 0; i < 4; ++i)
            if (i != idx4)
                v[kk++] = nodes[i];

        const Storage::real_array c1 = v[0].Coords();
        const Storage::real_array c2 = v[1].Coords();
        const Storage::real_array c3 = v[2].Coords();

        const double x1[3] = {c1[0], c1[1], c1[2]};
        const double x2[3] = {c2[0], c2[1], c2[2]};
        const double x3[3] = {c3[0], c3[1], c3[2]};

        const double at1 = v[0].Real(at_tag);
        const double at2 = v[1].Real(at_tag);
        const double at3 = v[2].Real(at_tag);

        const bool f1 = IsFiniteAT(at1);
        const bool f2 = IsFiniteAT(at2);
        const bool f3 = IsFiniteAT(at3);

        const int finite_count = (f1 ? 1 : 0) + (f2 ? 1 : 0) + (f3 ? 1 : 0);
        if (finite_count == 0) continue;

        auto metric_distance = [&](const double xs[3]) -> double
        {
            const double dx = x4[0] - xs[0];
            const double dy = x4[1] - xs[1];
            const double dz = x4[2] - xs[2];

            double q = QuadForm(sigma_inv, dx, dy, dz);

            // Для SPD-тензора q должно быть >= 0.
            // Небольшой минус возможен только из-за round-off.
            if (q < -1e-12)
            {
                return std::numeric_limits<double>::max();
            }

            if (q < 0.0) q = 0.0;
            return std::sqrt(q);
        };

        // ------------------------------
        // Случай: известна только 1 вершина
        // ------------------------------
        if (finite_count == 1)
        {
            const double* xs = nullptr;
            double at_s = 0.0;

            if (f1) { xs = x1; at_s = at1; }
            else if (f2) { xs = x2; at_s = at2; }
            else { xs = x3; at_s = at3; }

            const double dist = metric_distance(xs);
            if (dist == std::numeric_limits<double>::max()) continue;

            const double cand = at_s + dist;

            if (cand >= rt4 && cand < best_at)
                best_at = cand;

            continue;
        }

        // ------------------------------
        // Случай: известны 2 вершины
        // ------------------------------
        if (finite_count == 2)
        {
            const double* A = nullptr;
            const double* B = nullptr;
            double atA = 0.0;
            double atB = 0.0;

            if (f1 && f2)
            {
                A = x1; B = x2;
                atA = at1; atB = at2;
            }
            else if (f1 && f3)
            {
                A = x1; B = x3;
                atA = at1; atB = at3;
            }
            else
            {
                A = x2; B = x3;
                atA = at2; atB = at3;
            }

            double cell_best = std::numeric_limits<double>::max();

            for (int i = 0; i <= N_EDGE; ++i)
            {
                const double t = static_cast<double>(i) / static_cast<double>(N_EDGE);

                const double xs[3] = {
                    (1.0 - t)*A[0] + t*B[0],
                    (1.0 - t)*A[1] + t*B[1],
                    (1.0 - t)*A[2] + t*B[2]
                };

                const double at_s = (1.0 - t)*atA + t*atB;
                const double dist = metric_distance(xs);
                if (dist == std::numeric_limits<double>::max()) continue;

                const double cand = at_s + dist;
                if (cand < cell_best) cell_best = cand;
            }

            if (cell_best >= rt4 && cell_best < best_at)
                best_at = cell_best;

            continue;
        }

        // ------------------------------
        // Случай: известны все 3 вершины
        // ------------------------------
        double cell_best = std::numeric_limits<double>::max();

        for (int i = 0; i <= N_TRI; ++i)
        {
            for (int j = 0; j <= N_TRI - i; ++j)
            {
                const double l1 = static_cast<double>(i) / static_cast<double>(N_TRI);
                const double l2 = static_cast<double>(j) / static_cast<double>(N_TRI);
                const double l3 = 1.0 - l1 - l2;

                const double xs[3] = {
                    l1*x1[0] + l2*x2[0] + l3*x3[0],
                    l1*x1[1] + l2*x2[1] + l3*x3[1],
                    l1*x1[2] + l2*x2[2] + l3*x3[2]
                };

                const double at_s = l1*at1 + l2*at2 + l3*at3;
                const double dist = metric_distance(xs);
                if (dist == std::numeric_limits<double>::max()) continue;

                const double cand = at_s + dist;
                if (cand < cell_best) cell_best = cand;
            }
        }

        if (cell_best >= rt4 && cell_best < best_at)
            best_at = cell_best;
    }

    return (best_at == std::numeric_limits<double>::max()) ? INF_AT : best_at;
}