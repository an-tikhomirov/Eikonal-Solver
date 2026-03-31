#include "local_solver_isotropic.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

using namespace INMOST;

double local_solver_isotropic(INMOST::Node v4,
                            INMOST::Tag at_tag,
                            INMOST::Tag rt_tag)
{
    double best_phi = std::numeric_limits<double>::max();

    auto c4 = v4.Coords();
    double x4[3] = {c4[0], c4[1], c4[2]};

    ElementArray<Cell> cells = v4.getCells();

    for (int cell_idx = 0; cell_idx < cells.size(); ++cell_idx)
    {
        Cell cell = cells[cell_idx];
        ElementArray<Node> nodes = cell.getNodes();

        int idx = -1;
        for (int i = 0; i < 4; ++i)
        {
            if (nodes[i].GetHandle() == v4.GetHandle())
            {
                idx = i;
                break;
            }
        }
        if (idx == -1) continue;

        Node v1 = nodes[(idx + 1) % 4];
        Node v2 = nodes[(idx + 2) % 4];
        Node v3 = nodes[(idx + 3) % 4];

        const auto c1 = v1.Coords();
        const auto c2 = v2.Coords();
        const auto c3 = v3.Coords();

        double x1[3] = {c1[0], c1[1], c1[2]};   
        double x2[3] = {c2[0], c2[1], c2[2]};
        double x3[3] = {c3[0], c3[1], c3[2]};

        const double phi1 = v1.Real(at_tag);
        const double phi2 = v2.Real(at_tag);
        const double phi3 = v3.Real(at_tag);
        

        auto compute_phi = [&](double l1, double l2) -> double
        {
            const double l3 = 1.0 - l1 - l2;

            double xs[3];
            xs[0] = l1 * x1[0] + l2 * x2[0] + l3 * x3[0];
            xs[1] = l1 * x1[1] + l2 * x2[1] + l3 * x3[1];
            xs[2] = l1 * x1[2] + l2 * x2[2] + l3 * x3[2];

            const double dx = x4[0] - xs[0];
            const double dy = x4[1] - xs[1];
            const double dz = x4[2] - xs[2];

            const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

            return l1 * phi1 + l2 * phi2 + l3 * phi3 + dist;
        };

        double phi_min = std::numeric_limits<double>::max();
        const int N = 15;

        for (int i = 0; i <= N; ++i)
        {
            for (int j = 0; j <= N - i; ++j)
            {
                const double l1 = i / double(N);
                const double l2 = j / double(N);

                const double phi = compute_phi(l1, l2);
                if (phi < phi_min) phi_min = phi;
            }
        }

        if (phi_min < best_phi) best_phi = phi_min;
    }

    return (best_phi == std::numeric_limits<double>::max()) ? 1e30 : best_phi;
}