#include <iostream>
#include <string>
#include <array>

#include "Mesh.hpp"
#include "global_solver_isotropic.hpp"
#include "global_solver_anisotropic.hpp"

static bool HasArg(int argc, char* argv[], const std::string& s)
{
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == s) return true;
    return false;
}

static bool GetArgValue(int argc, char* argv[], const std::string& key, double& out)
{
    for (int i = 1; i + 1 < argc; ++i)
    {
        if (std::string(argv[i]) == key)
        {
            out = std::stod(argv[i + 1]);
            return true;
        }
    }
    return false;
}

static bool GetArgVec3(int argc, char* argv[], const std::string& key, std::array<double,3>& out)
{
    for (int i = 1; i + 3 < argc; ++i)
    {
        if (std::string(argv[i]) == key)
        {
            out[0] = std::stod(argv[i + 1]);
            out[1] = std::stod(argv[i + 2]);
            out[2] = std::stod(argv[i + 3]);
            return true;
        }
    }
    return false;
}

static void PrintUsage(const char* prog)
{
    std::cout
        << "Usage:\n  " << prog << " <input.vtu> <output.vtu> [--iso|--aniso|--aniso-x]\n"
        << "        [--erp <val>] [--s1 x y z] [--r1 <val>]\n\n"
        << "Defaults:\n"
        << "  mode: --aniso (if nothing specified)\n"
        << "  ERP=20, S1=(0,0,0) r1=1\n";
}

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        PrintUsage(argv[0]);
        return 1;
    }

    const bool run_iso     = HasArg(argc, argv, "--iso");
    const bool run_aniso_x = HasArg(argc, argv, "--aniso-x");
    const bool run_aniso   = HasArg(argc, argv, "--aniso") || (!run_iso && !run_aniso_x);
    (void)run_aniso;

    INMOST::Mesh::Initialize(&argc, &argv);

    try
    {
        Mesh mesh;

        if (!mesh.Load(argv[1]))
        {
            std::cerr << "Failed to load mesh: " << argv[1] << "\n";
            INMOST::Mesh::Finalize();
            return 1;
        }

        mesh.InitializeFields();

        double ERP = 20.0;
        (void)GetArgValue(argc, argv, "--erp", ERP);

        std::array<double,3> s1 = {28.7816, 17.8192, -5.42657};
        (void)GetArgVec3(argc, argv, "--s1", s1);

        double r1 = 0.1;
        (void)GetArgValue(argc, argv, "--r1", r1);

        const double lambda_L = 0.13175;
        const double lambda_T = 0.0176;

        std::cout << "ERP=" << ERP << "\n";
        std::cout << "S1 center=(" << s1[0] << "," << s1[1] << "," << s1[2] << ") r1=" << r1 << " t=0\n";

        mesh.ApplyStimulusSphere(s1, r1, 0.0);

        if (run_iso)
        {
            std::cout << "Mode: ISOTROPIC\n";
            mesh.SetIsotropicSpeed(1.0);

            global_solver_isotropic(mesh.GetMesh(),
                                    mesh.GetATTag(),
                                    mesh.GetRTTag(),
                                    mesh.GetActiveTag(),
                                    1000);
        }
        else
        {
            if (run_aniso_x)
            {
                std::cout << "Mode: ANISOTROPIC (test fiber = X)\n";
                mesh.BuildConductivityTensor(lambda_L, lambda_T, true);
            }
            else
            {
                std::cout << "Mode: ANISOTROPIC (from VTU fibers)\n";
                mesh.BuildConductivityTensor(lambda_L, lambda_T, false);
            }

            global_solver_anisotropic(mesh.GetMesh(),
                                      mesh.GetATTag(),
                                      mesh.GetRTTag(),
                                      mesh.GetActiveTag(),
                                      mesh.GetSigmaTag(),
                                      1000);
        }

        mesh.UpdateRefractory(ERP);

        if (!mesh.Save(argv[2]))
        {
            std::cerr << "Failed to save result: " << argv[2] << "\n";
            INMOST::Mesh::Finalize();
            return 1;
        }

        std::cout << "Result saved to: " << argv[2] << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << "\n";
        INMOST::Mesh::Finalize();
        return 1;
    }

    INMOST::Mesh::Finalize();
    return 0;
}
