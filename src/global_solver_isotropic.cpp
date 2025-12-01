#include "global_solver_isotropic.hpp"
#include "local_solver_isotropic.hpp"
#include <queue>
#include <set>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>

using namespace INMOST;

const double INF_PHI = 1.0e12;


void global_solver_isotropic(
    Mesh* mesh,
    Tag phi_tag,
    Tag active_tag,
    int max_iterations
) {
    const double TOLERANCE = 1e-6;
    
    std::queue<HandleType> active_list;
    std::set<HandleType> processed;
    
    // инициализировать очередь активных узлов и их соседей
    for (Mesh::iteratorNode it = mesh->BeginNode(); 
         it != mesh->EndNode(); ++it) {
        if (it->Real(phi_tag) < INF_PHI * 0.5) {
            it->Integer(active_tag) = 1;
            active_list.push(it->GetHandle());
            processed.insert(it->GetHandle());
            
            // добавить соседей источника сразу в очередь
            std::set<HandleType> neighbors_set;
            ElementArray<Cell> cells = it->getCells();
            for (int c = 0; c < cells.size(); ++c) {
                Cell cell = cells[c];
                ElementArray<Node> cell_nodes = cell.getNodes();
                for (int n = 0; n < cell_nodes.size(); ++n) {
                    neighbors_set.insert(cell_nodes[n].GetHandle());
                }
            }
            
            for (auto neighbor_handle : neighbors_set) {
                if (processed.find(neighbor_handle) == processed.end()) {
                    Node neighbor(mesh, neighbor_handle);
                    processed.insert(neighbor_handle);
                    active_list.push(neighbor_handle);
                }
            }
        } else {
            it->Integer(active_tag) = 0;
        }
    }
    
    std::cout << "\nИнициализация: " << active_list.size() 
              << " активных узлов\n";
    std::cout << "Запуск глобального солвера...\n\n";
    
    // Главный цикл
    int iteration = 0;
    int total_updates = 0;
    
    while (!active_list.empty() && iteration < max_iterations) {
        iteration++;
        
        int updates_this_iter = 0;
        int size_this_iter = active_list.size();
        
        // Обработать все узлы в текущей очереди
        for (int i = 0; i < size_this_iter; ++i) {
            HandleType node_handle = active_list.front();
            active_list.pop();
            
            Node node(mesh, node_handle);
            if (!node.isValid()) continue;
            
            double phi_old = node.Real(phi_tag);
            double phi_new = local_solver_isotropic(node, phi_tag);
            
            // Если улучшилось → обновить
            if (phi_new < phi_old - TOLERANCE) {
                node.Real(phi_tag) = phi_new;
                updates_this_iter++;
                total_updates++;
                
                // Получить соседних узлов через ячейки
                std::set<HandleType> neighbors_set;
                ElementArray<Cell> cells = node.getCells();
                for (int c = 0; c < cells.size(); ++c) {
                    Cell cell = cells[c];
                    ElementArray<Node> cell_nodes = cell.getNodes();
                    for (int n = 0; n < cell_nodes.size(); ++n) {
                        neighbors_set.insert(cell_nodes[n].GetHandle());
                    }
                }
                
                // исправить, добавить active_set
                for (auto neighbor_handle : neighbors_set) {
                    if (processed.find(neighbor_handle) == processed.end()) {
                        Node neighbor(mesh, neighbor_handle);
                        processed.insert(neighbor_handle);
                        active_list.push(neighbor_handle);
                    }
                }
            }
        }
        
        std::cout << "Итерация " << iteration << ": " << updates_this_iter 
                  << " обновлений, в очереди " << active_list.size() 
                  << " узлов\n";
        
        // Если нет обновлений → сходимость
        if (updates_this_iter == 0) {
            std::cout << "Сходимость достигнута!\n";
            break;
        }
    }
    
    // Статистика
    std::cout << "\n═══════════════════════════════════════════════════════\n";
    std::cout << "Глобальный солвер завершён\n";
    std::cout << "Итераций: " << iteration << "\n";
    std::cout << "Всего обновлений: " << total_updates << "\n";
    
    double phi_min = INF_PHI;
    double phi_max = 0;
    int computed_nodes = 0;
    
    for (Mesh::iteratorNode it = mesh->BeginNode(); 
         it != mesh->EndNode(); ++it) {
        double phi = it->Real(phi_tag);
        
        if (phi < INF_PHI * 0.9) {
            computed_nodes++;
            phi_min = std::min(phi_min, phi);
            phi_max = std::max(phi_max, phi);
        }
    }
    
    std::cout << "Вычисленных узлов: " << computed_nodes << "\n";
    std::cout << "Min φ: " << phi_min << ", Max φ: " << phi_max << "\n";
    std::cout << "═══════════════════════════════════════════════════════\n\n";
}