#include "Mesh.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>

const double INF_PHI = 1.0e12;

Mesh::Mesh()
{
    mesh = new INMOST::Mesh();
}

Mesh::~Mesh()
{
    delete mesh;
}

bool Mesh::Load(const std::string& filename)
{
    mesh->SetFileOption("VERBOSITY", "2");

    try
    {
        mesh->Load(filename);

        std::cout << "Mesh loaded: "
                  << mesh->NumberOfNodes() << " nodes, "
                  << mesh->NumberOfCells() << " cells"
                  << std::endl;

        return true;
    }
    catch (...)
    {
        std::cerr << "Failed to load mesh: " << filename << std::endl;
        return false;
    }
}

bool Mesh::Save(const std::string& filename)
{
    try
    {
        mesh->Save(filename);
        std::cout << "Mesh saved: " << filename << std::endl;
        return true;
    }
    catch (...)
    {
        std::cerr << "Failed to save mesh: " << filename << std::endl;
        return false;
    }
}

void Mesh::InitializeFields()
{
    at_tag = mesh->CreateTag("AT", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);
    rt_tag = mesh->CreateTag("RT", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);
    di_tag = mesh->CreateTag("DI", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);

    active_tag = mesh->CreateTag("FIM_active", INMOST::DATA_INTEGER, INMOST::NODE, INMOST::NONE, 1);
    cycle_tag  = mesh->CreateTag("FIM_cycle",  INMOST::DATA_INTEGER, INMOST::NODE, INMOST::NONE, 1);

    vm_tag = mesh->CreateTag("Vm", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);
    w_tag  = mesh->CreateTag("w",  INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);

    idiff_tag = mesh->CreateTag("Idiff", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);

    stim_tag = mesh->CreateTag("FIM_stim", INMOST::DATA_INTEGER, INMOST::NODE, INMOST::NONE, 1);

    excited_tag = mesh->CreateTag("FIM_excited", INMOST::DATA_INTEGER, INMOST::NODE, INMOST::NONE, 1);
    src_tag     = mesh->CreateTag("FIM_src",     INMOST::DATA_INTEGER, INMOST::NODE, INMOST::NONE, 1);

    for (INMOST::Mesh::iteratorNode inode = mesh->BeginNode();
         inode != mesh->EndNode(); ++inode)
    {
        inode->Real(at_tag) = INF_PHI;
        inode->Real(rt_tag) = -INF_PHI;
        inode->Real(di_tag) = INF_PHI;

        inode->Integer(active_tag) = 0;
        inode->Integer(cycle_tag)  = 0;

        inode->Real(vm_tag) = -80.0;
        inode->Real(w_tag)  = 0.0;
        inode->Real(idiff_tag) = 0.0;

        inode->Integer(stim_tag) = 0;

        inode->Integer(excited_tag) = 0;
        inode->Integer(src_tag) = 0;
    }

    std::cout << "Fields initialized: AT/RT/DI + Vm/w (+diagnostics)\n";
}

void Mesh::SetIsotropicSpeed(double speed)
{
    speed_tag = mesh->CreateTag("Speed",
                                INMOST::DATA_REAL,
                                INMOST::CELL,
                                INMOST::NONE,
                                1);

    for (INMOST::Mesh::iteratorCell icell = mesh->BeginCell();
         icell != mesh->EndCell(); ++icell)
    {
        icell->Real(speed_tag) = speed;
    }

    std::cout << "Isotropic speed set to: " << speed << std::endl;
}

static inline void Normalize3(double f[3])
{
    const double n = std::sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
    if (n > 0.0)
    {
        f[0] /= n;
        f[1] /= n;
        f[2] /= n;
    }
}

void Mesh::BuildConductivityTensor(double lambda_longitudinal,
                                   double lambda_transversal,
                                   bool force_fiber_x)
{
    sigma_tag = mesh->CreateTag("sigma_tensor",
                                INMOST::DATA_REAL,
                                INMOST::CELL,
                                INMOST::NONE,
                                9);

    INMOST::Tag fibers_tag;

    if (!force_fiber_x)
    {
        fibers_tag = mesh->GetTag("fibers");
    }

    const double a = lambda_longitudinal - lambda_transversal;

    for (INMOST::Mesh::iteratorCell it = mesh->BeginCell();
         it != mesh->EndCell(); ++it)
    {
        double f[3] = {1.0, 0.0, 0.0};

        if (!force_fiber_x)
        {
            INMOST::Storage::real_array fa = it->RealArray(fibers_tag);

            f[0] = fa[0];
            f[1] = fa[1];
            f[2] = fa[2];

            Normalize3(f);
        }

        INMOST::Storage::real_array sig = it->RealArray(sigma_tag);

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                sig[i*3 + j] = a * f[i] * f[j];

        sig[0] += lambda_transversal;
        sig[4] += lambda_transversal;
        sig[8] += lambda_transversal;
    }

    std::cout << "Conductivity tensor sigma built.\n"
              << "lambda_L = " << lambda_longitudinal << "\n"
              << "lambda_T = " << lambda_transversal << "\n"
              << "force_fiber_x = " << (force_fiber_x ? "true" : "false")
              << std::endl;
}

void Mesh::ApplyStimulusSphere(const std::array<double,3>& center,
                               double radius,
                               double t_stim)
{
    const double r2 = radius * radius;

    for (INMOST::Mesh::iteratorNode it = mesh->BeginNode();
         it != mesh->EndNode(); ++it)
    {
        auto c = it->Coords();

        double dx = c[0] - center[0];
        double dy = c[1] - center[1];
        double dz = c[2] - center[2];

        double d2 = dx*dx + dy*dy + dz*dz;

        if (d2 <= r2)
        {
            it->Integer(src_tag) = 1;
            it->Integer(stim_tag) = 1;
            it->Real(at_tag) = t_stim;
        }
        else
        {
            it->Integer(src_tag) = 0;
            it->Integer(stim_tag) = 0;
        }
    }
}

void Mesh::UpdateRefractory(double ERP)
{
    for (INMOST::Mesh::iteratorNode it = mesh->BeginNode();
         it != mesh->EndNode(); ++it)
    {
        const double at = it->Real(at_tag);

        if (at < INF_PHI*0.5)
            it->Real(rt_tag) = at + ERP;
    }
}