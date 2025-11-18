#include <iostream>
#include "Mesh.hpp"

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <input.vtu> <output.vtu>" << std::endl;
        return 1;
    }

    INMOST::Mesh::Initialize(&argc, &argv);
    
    try {
        Mesh mesh;
        if (!mesh.Load(argv[1])) {
            std::cerr << "Failed to load mesh: " << argv[1] << std::endl;
            return 1;
        }
        

        mesh.InitializeEikonalData({0.5, 0.5, 0.5}, 0.1);
        mesh.SetIsotropicSpeed(1.0);
        
        std::cout << "Mesh initialized successfully!" << std::endl;
        
        if (!mesh.Save(argv[2])) {
            std::cerr << "Failed to save result: " << argv[2] << std::endl;
            return 1;
        }
        
        std::cout << "Result saved to: " << argv[2] << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }
    
    INMOST::Mesh::Finalize();
    return 0;
}