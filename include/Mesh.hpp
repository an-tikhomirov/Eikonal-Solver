// include/Mesh.hpp
#pragma once
#include <inmost.h>
#include <string>
#include <array>

class Mesh {
public:
    Mesh();
    ~Mesh();

    bool Load(const std::string& filename);
    bool Save(const std::string& filename);
    

    void InitializeEikonalData(const std::array<double, 3>& source_point, double source_radius = 0.1);
    

    void SetIsotropicSpeed(double speed);

    INMOST::Mesh* GetMesh() { return mesh; }
    INMOST::Tag GetPhiTag() { return phi_tag; }
    INMOST::Tag GetActiveTag() { return active_tag; }
    INMOST::Tag GetSpeedTag() { return speed_tag; }

private:
    INMOST::Mesh* mesh; 

    INMOST::Tag phi_tag;
    INMOST::Tag active_tag; 
    INMOST::Tag speed_tag;
};