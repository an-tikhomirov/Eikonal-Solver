#pragma once

#include <string>
#include <array>
#include "inmost.h"

// Mesh.hpp (фрагмент)

class Mesh
{
private:
    INMOST::Mesh* mesh;

    INMOST::Tag at_tag;        // activation time
    INMOST::Tag rt_tag;        // repolarization time    
    INMOST::Tag di_tag;        // diagnostic interval
    INMOST::Tag active_tag;    // 0 или 1
    INMOST::Tag cycle_tag;     // номер цикла
    INMOST::Tag vm_tag;        // transmembrane voltage
    INMOST::Tag w_tag;         // 
    INMOST::Tag idiff_tag;     // 
    INMOST::Tag stim_tag;      // 
 
    INMOST::Tag excited_tag;   // 
    INMOST::Tag src_tag;       // 

    INMOST::Tag speed_tag;     // 
    INMOST::Tag sigma_tag;     // 

public:
    Mesh();
    ~Mesh();

    INMOST::Mesh* GetMesh() { return mesh; }

    INMOST::Tag GetATTag() const { return at_tag; }
    INMOST::Tag GetRTTag() const { return rt_tag; }
    INMOST::Tag GetDITag() const { return di_tag; }
    INMOST::Tag GetActiveTag() const { return active_tag; }
    INMOST::Tag GetCycleTag() const { return cycle_tag; }

    INMOST::Tag GetVmTag() const { return vm_tag; }
    INMOST::Tag GetWTag() const { return w_tag; }
    INMOST::Tag GetIdiffTag() const { return idiff_tag; }
    INMOST::Tag GetStimTag() const { return stim_tag; }

    INMOST::Tag GetExcitedTag() const { return excited_tag; }
    INMOST::Tag GetSrcTag() const { return src_tag; }

    INMOST::Tag GetSpeedTag() const { return speed_tag; }
    INMOST::Tag GetSigmaTag() const { return sigma_tag; }

    void InitializeFields();

    bool Load(const std::string& filename);
    bool Save(const std::string& filename);

    void SetIsotropicSpeed(double speed);

    void BuildConductivityTensor(double lambda_longitudinal,
                                 double lambda_transversal,
                                 bool force_fiber_x);
                                 
    void ApplyStimulusSphere(const std::array<double,3>& center,
                             double radius,
                             double t_stim);

    void UpdateRefractory(double ERP);
};