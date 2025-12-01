#include "local_solver_isotropic.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

using namespace INMOST;

double local_solver_isotropic(Node v4, Tag phi_tag) {
    double best_phi = std::numeric_limits<double>::max();
    auto c4 = v4.Coords();
    double x4[3] = {c4[0], c4[1], c4[2]};
    // получаем все смежные тетраэдры для которых v4 вершина
    ElementArray<Cell> cells = v4.getCells();
    
    for (int cell_idx = 0; cell_idx < cells.size(); ++cell_idx) {
        Cell cell = cells[cell_idx];
        ElementArray<Node> nodes = cell.getNodes();
        
        int idx = -1;
        for (int i = 0; i < 4; ++i) {
            if (nodes[i].GetHandle() == v4.GetHandle()) { // дескриптор типа Handle
                idx = i;
                break;
            }
        }
        if (idx == -1) continue;
        
        Node v1 = nodes[(idx + 1) % 4];
        Node v2 = nodes[(idx + 2) % 4];
        Node v3 = nodes[(idx + 3) % 4];
        
        auto c1 = v1.Coords(); //const добавить может
        auto c2 = v2.Coords();
        auto c3 = v3.Coords();
        
        double x1[3] = {c1[0], c1[1], c1[2]};
        double x2[3] = {c2[0], c2[1], c2[2]};
        double x3[3] = {c3[0], c3[1], c3[2]};
        
        double phi1 = v1.Real(phi_tag);
        double phi2 = v2.Real(phi_tag);
        double phi3 = v3.Real(phi_tag);
        
        auto compute_phi = [&](double l1, double l2) -> double {
            double xs[3];
            double l3 = 1.0 - l1 - l2;
            xs[0] = l1*x1[0] + l2*x2[0] + l3*x3[0];
            xs[1] = l1*x1[1] + l2*x2[1] + l3*x3[1];
            xs[2] = l1*x1[2] + l2*x2[2] + l3*x3[2];
            
            double dx = x4[0] - xs[0];
            double dy = x4[1] - xs[1];
            double dz = x4[2] - xs[2];
            
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            
            return l1*phi1 + l2*phi2 + l3*phi3 + dist;
        };
        
        double phi_min = std::numeric_limits<double>::max();
        const int N = 15;
        
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N - i; ++j) {
                double l1 = i / double(N);
                double l2 = j / double(N);
                double phi = compute_phi(l1, l2);
                if (phi < phi_min) phi_min = phi;
            }
        }
        
        // вдоль трёх рёбер, избыточно
        double dx1 = x4[0] - x1[0], dy1 = x4[1] - x1[1], dz1 = x4[2] - x1[2];
        double dx2 = x4[0] - x2[0], dy2 = x4[1] - x2[1], dz2 = x4[2] - x2[2];
        double dx3 = x4[0] - x3[0], dy3 = x4[1] - x3[1], dz3 = x4[2] - x3[2];
        
        double dist1 = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
        double dist2 = std::sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
        double dist3 = std::sqrt(dx3*dx3 + dy3*dy3 + dz3*dz3);
        
        double phi_edge1 = phi1 + dist1;
        double phi_edge2 = phi2 + dist2;
        double phi_edge3 = phi3 + dist3;
        
        double phi_fallback = std::min({phi_edge1, phi_edge2, phi_edge3});
        
        // выбрать минимум
        double cell_candidate = std::min(phi_min, phi_fallback);
        if (cell_candidate < best_phi)
            best_phi = cell_candidate;
    }
    
    return (best_phi == std::numeric_limits<double>::max()) ? 1e30 : best_phi;
}