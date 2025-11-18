#include "Mesh.hpp"
#include <iostream>
#include <cmath>

Mesh::Mesh() {
    mesh = new INMOST::Mesh();
}

Mesh::~Mesh() {
    delete mesh;
}

bool Mesh::Load(const std::string& filename) {
    mesh->SetFileOption("VERBOSITY", "2");
    try {
        mesh->Load(filename);
        std::cout << "Mesh loaded: " << mesh->NumberOfNodes() << " nodes, " 
                  << mesh->NumberOfCells() << " cells" << std::endl;
        return true;
    } catch (...) {
        std::cerr << "Failed to load mesh: " << filename << std::endl;
        return false;
    }
}

bool Mesh::Save(const std::string& filename) {
    try {
        mesh->Save(filename);
        std::cout << "Mesh saved: " << filename << std::endl;
        return true;
    } catch (...) {
        std::cerr << "Failed to save mesh: " << filename << std::endl;
        return false;
    }
}

void Mesh::InitializeEikonalData(const std::array<double, 3>& source_point, double source_radius) {
    phi_tag = mesh->CreateTag("Phi", INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);
    active_tag = mesh->CreateTag("Active", INMOST::DATA_INTEGER, INMOST::NODE, INMOST::NONE, 1);

    for (INMOST::Mesh::iteratorNode inode = mesh->BeginNode(); inode != mesh->EndNode(); ++inode) {
        double coords[3];
        inode->Barycenter(coords);
        
        double dist_sq = 0.0;
        for (int i = 0; i < 3; i++) {
            double diff = coords[i] - source_point[i];
            dist_sq += diff * diff;
        }
        
        if (dist_sq < source_radius * source_radius) {
            inode->Real(phi_tag) = 0.0;
            inode->Integer(active_tag) = 1;
        } else {
            inode->Real(phi_tag) = 1.0e12;
            inode->Integer(active_tag) = 0;
        }
    }
    
    std::cout << "Eikonal data initialized with source at (" 
              << source_point[0] << ", " << source_point[1] << ", " << source_point[2] 
              << "), radius = " << source_radius << std::endl;
}

void Mesh::SetIsotropicSpeed(double speed) {
    speed_tag = mesh->CreateTag("Speed", INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
    
    for (INMOST::Mesh::iteratorCell inode = mesh->BeginCell(); inode != mesh->EndCell(); ++inode) {
        inode->Real(speed_tag) = speed;
    }
    
    std::cout << "Isotropic speed set to: " << speed << std::endl;
}