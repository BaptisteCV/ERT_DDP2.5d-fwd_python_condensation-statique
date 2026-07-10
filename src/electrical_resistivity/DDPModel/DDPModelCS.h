#ifndef DDPMODELCS_H_
#define DDPMODELCS_H_

#include <memory>
#include "../EPMModel/EPMDefinitionCSEigen.h"
#include "../Utilitaries/Storage/EigenCSStructures.h"
#include "../DFNModel/NetworkMeshes.h"
#include "../InputOutput/Parameters.h"
#include "../InputOutput/BoundaryConditions.h"
#include <vector>
#include <map>

namespace eigen_cs {


class DDPModelCS : public EPMSystemCSEigen {
public:
    
    int nb_nodes_fract;
    SubNetworkMap subnetwork_map;
    BoundaryConditionsDFN bc_map_dfn;

    Domain domain_stored;
    BoundaryConditionsDef bc_map_classic;
    BordersMap global_border_map;
    NodesMap global_nodes_map;
    
    
    int Nfi;
    int Nfb;
    std::vector<int> interior_fract_global_indices;
    std::vector<int> border_fract_global_indices;
    std::map<int, int> global_to_interior_fract;
    std::map<int, int> global_to_border_fract;

    
    
    SparseMatrix block_Afi_i2i2;
    SparseMatrix block_Afi_b2b2;
    SparseMatrix block_Afi_i2b2;
    SparseMatrix block_Afi_b2i2;

    SparseMatrix block_Afm_i2i1;
    SparseMatrix block_Afm_i2b1;
    SparseMatrix block_Afm_b2i1;
    SparseMatrix block_Afm_b2b1;
    
    SparseMatrix block_Amf_i1i2;
    SparseMatrix block_Amf_i1b2;
    SparseMatrix block_Amf_b1i2;
    SparseMatrix block_Amf_b1b2;

    SparseMatrix Aii_g;
    SparseMatrix Aib_g;
    SparseMatrix Abi_g;
    SparseMatrix Abb_g;

    DenseMatrix matrix_Y_DDP;
    DenseMatrix matrix_S_DDP_precomputed;
    DenseMatrix matrix_S_DDP;
    
    Vector vector_delta_diag_DDP;

    Vector vector_zi_DDP;
    Vector vector_zb_DDP;
    Vector vector_rhs_DDP;
    Vector vector_phi_i_DDP;
    Vector vector_phi_b_DDP;
    

    Vector vector_b_fi;
    Vector vector_b_fb;

    DenseMatrix potential_fracture;

    double current_w;
    bool singularity;


    std::shared_ptr<LinearSystemSparseLLTFactorization> factor_Ai1plusi2;
    bool factor_Afi_i2i2_ready;

public:

    DDPModelCS();
    
    DDPModelCS(const Parameters& param, 
               const Domain& domain, 
               const NetworkMeshes& net_mesh,
               bool fourier = false,
               bool singularity_treatment = false);
    
    DDPModelCS(const DDPModelCS& other);

    DDPModelCS& operator=(const DDPModelCS& other);
    virtual ~DDPModelCS() = default;

    void InitializeDFNSCIndexing();

    void DetermineBasicLinearSystemDDPCS();
    void UpdateDDPSCBlocksFourier();
    void UpdateDDPSCBoundaryConditionsAndSources(bool fourier);
    void PrecomputeDDPSCForCurrentFrequency();
    void PrecomputeDDPSCForCurrentFrequency2();
    void SolveCondensedDDPSystem();
    void ReconstructDDPInteriorFromBorder();

    Vector BuildGlobalDDPSolutionVector() const;
    DenseMatrix PotentialComputationDDPCS();

    bool IsInitializedDFNSC() const;
    void ExportDDPPotential(const std::string& output_dir, 
                              int elec_nb, 
                              double current_w);
};

}  // namespace eigen_cs

#endif /* DDPMODELCS_H_ */
