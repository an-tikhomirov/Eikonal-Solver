#include "global_solver_anisotropic.hpp"
#include "local_solver_anisotropic.hpp"

#include <queue>
#include <iostream>

using namespace INMOST;

static constexpr double INF_AT = 1.0e12;
static inline bool IsFiniteAT(double at) { return at < INF_AT * 0.5; }

void global_solver_anisotropic(
    Mesh* mesh,
    Tag at_tag,
    Tag rt_tag,
    Tag active_tag,
    Tag sigma_tag,
    int max_iterations
)
{
    const double TOL = 1e-6;

    std::queue<HandleType> active_queue;

    // Сбрасываем active_tag у всех узлов
    for (Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
        it->Integer(active_tag) = 0;

    // Инициализация очереди:
    // если в тетраэдре есть и конечные, и бесконечные AT,
    // то бесконечные узлы добавляем в очередь
    for (Mesh::iteratorCell icell = mesh->BeginCell(); icell != mesh->EndCell(); ++icell)
    {
        ElementArray<Node> nodes = icell->getNodes();

        int finite_count = 0;
        for (int i = 0; i < 4; ++i)
        {
            const double at = nodes[i].Real(at_tag);
            if (IsFiniteAT(at)) ++finite_count;
        }

        // если все узлы уже конечны или все бесконечны, этот тетраэдр не даёт стартовых кандидатов
        if (finite_count == 0 || finite_count == 4) continue;

        for (int i = 0; i < 4; ++i)
        {
            const double at = nodes[i].Real(at_tag);
            if (IsFiniteAT(at)) continue;

            if (nodes[i].Integer(active_tag) == 0)
            {
                nodes[i].Integer(active_tag) = 1;
                active_queue.push(nodes[i].GetHandle());
            }
        }
    }

    std::cout << "\n[ANISO] init queue: " << active_queue.size() << " nodes\n";

    int iteration = 0;
    int total_updates = 0;

    while (!active_queue.empty() && iteration < max_iterations)
    {
        ++iteration;

        int updates_this_iter = 0;
        const int batch_size = static_cast<int>(active_queue.size());

        for (int i = 0; i < batch_size; ++i)
        {
            const HandleType h = active_queue.front();
            active_queue.pop();

            Node node(mesh, h);
            if (!node.isValid()) continue;

            node.Integer(active_tag) = 0;

            const double at_old = node.Real(at_tag);
            const double at_new = local_solver_anisotropic(node, at_tag, rt_tag, sigma_tag);

            if (at_new < at_old - TOL)
            {
                node.Real(at_tag) = at_new;
                ++updates_this_iter;
                ++total_updates;

                ElementArray<Cell> cells = node.getCells();
                for (int c = 0; c < cells.size(); ++c)
                {
                    ElementArray<Node> cn = cells[c].getNodes();
                    for (int n = 0; n < cn.size(); ++n)
                    {
                        Node nb = cn[n];
                        if (!nb.isValid()) continue;

                        const double at_nb = nb.Real(at_tag);

                        if (!IsFiniteAT(at_nb) || at_nb > at_new)
                        {
                            if (nb.Integer(active_tag) == 0)
                            {
                                nb.Integer(active_tag) = 1;
                                active_queue.push(nb.GetHandle());
                            }
                        }
                    }
                }
            }
        }

        std::cout << "[ANISO] iter " << iteration
                  << ": updates " << updates_this_iter
                  << ", queue " << active_queue.size() << "\n";

        if (updates_this_iter == 0)
        {
            std::cout << "[ANISO] converged\n";
            break;
        }
    }

    std::cout << "[ANISO] done. iters=" << iteration
              << " total_updates=" << total_updates << "\n";
}