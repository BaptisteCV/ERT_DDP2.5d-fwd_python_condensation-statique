/*
 * DDPModelCS.cpp
 *
 *  Skeleton implementation for static condensation (DDPModelCS)
 */

#include "DDPModelCS.h"
#include "DDPDFNFunctions.h"
#include "ExchangeFunctions.h"
#include "SystLinearFunctions.h"
#include "../DFNModel/NetworkMeshes.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <ctime>
#include <omp.h>

namespace eigen_cs {

// Default constructor
DDPModelCS::DDPModelCS()
: EPMSystemCSEigen()
{
    current_w = 0.0;
    singularity = false;
    factor_Afi_i2i2_ready = false;
    factor_Ai1plusi2 = std::make_shared<LinearSystemSparseLLTFactorization>();
}

// Parametrized constructor
DDPModelCS::DDPModelCS(const Parameters& param,
                       const Domain& domain,
                       const NetworkMeshes& net_mesh,
                       bool fourier,
                       bool singularity_treatment)
: EPMSystemCSEigen()
{
    Nx = param.simu_param.DDP_param.Nx;
    Ny = param.simu_param.DDP_param.Ny;
    delta_x = param.domain_param.Lx / Nx;
    delta_y = param.domain_param.Ly / Ny;
    potential = DenseMatrix::Zero(Nx, Ny);

    if (param.domain_param.mat_cond != -1) {
        porous_cond = DenseMatrix::Constant(Nx, Ny, param.domain_param.mat_cond);
    } else {
        porous_cond = DenseMatrix::Zero(Nx, Ny);
        for (int i = 0; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j)
                porous_cond(i, j) = param.domain_param.mat_cond_dist(i, j);
    }

    current_w = 0.0;
    singularity = singularity_treatment;
    factor_Afi_i2i2_ready = false;
    factor_Ai1plusi2 = std::make_shared<LinearSystemSparseLLTFactorization>();

    //double t_ctor = omp_get_wtime();
    domain_stored = domain;
    bc_map_classic = param.simu_param.bc_map;
    SetTopBorderCondensation(bc_map_classic);
    global_border_map = net_mesh.border_map;
    global_nodes_map = net_mesh.nodes_map;
    //std::cout << "[PROFILE][DDPModelCS ctor] copies domain/bc/maps: "
    //          << omp_get_wtime() - t_ctor << " sec" << std::endl;

    //double t_netcopy = omp_get_wtime();
    NetworkMeshes net_mesh_copy = net_mesh;
    //std::cout << "[PROFILE][DDPModelCS ctor] NetworkMeshes copy: "
    //          << omp_get_wtime() - t_netcopy << " sec" << std::endl;

    //double t_ddpnet = omp_get_wtime();
    subnetwork_map = DDPNetworkMeshes(net_mesh_copy, Nx, Ny);
    nb_nodes_fract = ReturnNbNodes(subnetwork_map);
    //std::cout << "[PROFILE][DDPModelCS ctor] DDPNetworkMeshes + ReturnNbNodes: "
    //          << omp_get_wtime() - t_ddpnet << " sec" << std::endl;
}

DDPModelCS::DDPModelCS(const DDPModelCS& other)
    : EPMSystemCSEigen(other),
      nb_nodes_fract(other.nb_nodes_fract),
      subnetwork_map(other.subnetwork_map),
      bc_map_dfn(other.bc_map_dfn),
      domain_stored(other.domain_stored),
      bc_map_classic(other.bc_map_classic),
      global_border_map(other.global_border_map),
      global_nodes_map(other.global_nodes_map),
      Nfi(other.Nfi), Nfb(other.Nfb),
      interior_fract_global_indices(other.interior_fract_global_indices),
      border_fract_global_indices(other.border_fract_global_indices),
      global_to_interior_fract(other.global_to_interior_fract),
      global_to_border_fract(other.global_to_border_fract),
      block_Afi_i2i2(other.block_Afi_i2i2), block_Afi_b2b2(other.block_Afi_b2b2),
      block_Afi_i2b2(other.block_Afi_i2b2), block_Afi_b2i2(other.block_Afi_b2i2),
      block_Afm_i2i1(other.block_Afm_i2i1), block_Afm_i2b1(other.block_Afm_i2b1),
      block_Afm_b2i1(other.block_Afm_b2i1), block_Afm_b2b1(other.block_Afm_b2b1),
      block_Amf_i1i2(other.block_Amf_i1i2), block_Amf_i1b2(other.block_Amf_i1b2),
      block_Amf_b1i2(other.block_Amf_b1i2), block_Amf_b1b2(other.block_Amf_b1b2),
      Aii_g(other.Aii_g), Aib_g(other.Aib_g), Abi_g(other.Abi_g), Abb_g(other.Abb_g),
      matrix_Y_DDP(other.matrix_Y_DDP),
      matrix_S_DDP_precomputed(other.matrix_S_DDP_precomputed),
      matrix_S_DDP(other.matrix_S_DDP),
      vector_delta_diag_DDP(other.vector_delta_diag_DDP),
      vector_zi_DDP(other.vector_zi_DDP), vector_zb_DDP(other.vector_zb_DDP),
      vector_rhs_DDP(other.vector_rhs_DDP),
      vector_phi_i_DDP(other.vector_phi_i_DDP), vector_phi_b_DDP(other.vector_phi_b_DDP),
      vector_b_fi(other.vector_b_fi), vector_b_fb(other.vector_b_fb),
      potential_fracture(other.potential_fracture),
      current_w(other.current_w), singularity(other.singularity),
      factor_Ai1plusi2(other.factor_Ai1plusi2),   // shared_ptr : partage la factorisation (lecture seule)
      factor_Afi_i2i2_ready(other.factor_Afi_i2i2_ready)
{}

DDPModelCS& DDPModelCS::operator=(const DDPModelCS& other) {
    if (this != &other) {
        EPMSystemCSEigen::operator=(other);
        nb_nodes_fract = other.nb_nodes_fract;
        subnetwork_map = other.subnetwork_map;
        bc_map_dfn = other.bc_map_dfn;
        domain_stored = other.domain_stored;
        bc_map_classic = other.bc_map_classic;
        global_border_map = other.global_border_map;
        global_nodes_map = other.global_nodes_map;
        Nfi = other.Nfi; Nfb = other.Nfb;
        interior_fract_global_indices = other.interior_fract_global_indices;
        border_fract_global_indices = other.border_fract_global_indices;
        global_to_interior_fract = other.global_to_interior_fract;
        global_to_border_fract = other.global_to_border_fract;
        block_Afi_i2i2 = other.block_Afi_i2i2; block_Afi_b2b2 = other.block_Afi_b2b2;
        block_Afi_i2b2 = other.block_Afi_i2b2; block_Afi_b2i2 = other.block_Afi_b2i2;
        block_Afm_i2i1 = other.block_Afm_i2i1; block_Afm_i2b1 = other.block_Afm_i2b1;
        block_Afm_b2i1 = other.block_Afm_b2i1; block_Afm_b2b1 = other.block_Afm_b2b1;
        block_Amf_i1i2 = other.block_Amf_i1i2; block_Amf_i1b2 = other.block_Amf_i1b2;
        block_Amf_b1i2 = other.block_Amf_b1i2; block_Amf_b1b2 = other.block_Amf_b1b2;
        Aii_g = other.Aii_g; Aib_g = other.Aib_g;
        Abi_g = other.Abi_g; Abb_g = other.Abb_g;
        matrix_Y_DDP = other.matrix_Y_DDP;
        matrix_S_DDP_precomputed = other.matrix_S_DDP_precomputed;
        matrix_S_DDP = other.matrix_S_DDP;
        vector_delta_diag_DDP = other.vector_delta_diag_DDP;
        vector_zi_DDP = other.vector_zi_DDP; vector_zb_DDP = other.vector_zb_DDP;
        vector_rhs_DDP = other.vector_rhs_DDP;
        vector_phi_i_DDP = other.vector_phi_i_DDP; vector_phi_b_DDP = other.vector_phi_b_DDP;
        vector_b_fi = other.vector_b_fi; vector_b_fb = other.vector_b_fb;
        potential_fracture = other.potential_fracture;
        current_w = other.current_w; singularity = other.singularity;
        factor_Ai1plusi2 = other.factor_Ai1plusi2;
        factor_Afi_i2i2_ready = other.factor_Afi_i2i2_ready;
    }
    return *this;
}

// Destructor is defaulted in the header

// Initialization of DFN indexing for the condensed system
void DDPModelCS::InitializeDFNSCIndexing() {
    //double t_index_total = omp_get_wtime();

    // 1. Start from the EPM condensed-grid indexing.
    //    This tells us which porous blocks are interior vs exterior.
    //    We keep the inherited EPM indexing untouched.
    //double t_sc = omp_get_wtime();
    InitializeSCIndexing();
    //std::cout << "[PROFILE][InitializeDFNSCIndexing] InitializeSCIndexing: "
    //          << omp_get_wtime() - t_sc << " sec" << std::endl;

    // 2. Rebuild the DFN numbering from scratch.
    //    Every fracture node gets a global index, but only nodes located
    //    on the outer boundary of the global domain are classified as border.
    std::set<int> all_nodes;
    std::set<int> border_nodes;

    interior_fract_global_indices.clear();
    border_fract_global_indices.clear();
    global_to_interior_fract.clear();
    global_to_border_fract.clear();

    // 3. Traverse all porous blocks.
    //    - interior porous blocks: all associated fracture nodes are interior.
    //    - exterior porous blocks: fracture nodes are border only if they lie on
    //      the global domain boundary; otherwise they remain interior.
    for (SubNetworkMap::iterator it = subnetwork_map.begin(); it != subnetwork_map.end(); ++it) {
        NetworkMeshes &net = it->second;
        const bool porous_block_is_interior = (global_to_interior.find(it->first) != global_to_interior.end());

        for (std::vector<FractureMesh>::iterator im = net.meshes.begin(); im != net.meshes.end(); ++im) {
            int node_ids[2] = { im->p_ori.index, im->p_tar.index };

            for (int k = 0; k < 2; ++k) {
                const int node_id = node_ids[k];
                all_nodes.insert(node_id);

                // Avoid changing an already assigned node: a fracture node can belong
                // to several porous blocks, but its class and global index must stay stable.
                if (global_to_interior_fract.find(node_id) != global_to_interior_fract.end() ||
                    global_to_border_fract.find(node_id) != global_to_border_fract.end()) {
                    continue;
                }

                bool is_global_border_node = (net.border_map.find(node_id) != net.border_map.end());

                if (porous_block_is_interior) {
                    // Interior porous block: all its fracture nodes are interior.
                    global_to_interior_fract[node_id] = (int)interior_fract_global_indices.size();
                    interior_fract_global_indices.push_back(node_id);
                } else {
                    // Exterior porous block: only nodes on the global boundary are border.
                    if (is_global_border_node) {
                        global_to_border_fract[node_id] = (int)border_fract_global_indices.size();
                        border_fract_global_indices.push_back(node_id);
                        border_nodes.insert(node_id);
                    } else {
                        global_to_interior_fract[node_id] = (int)interior_fract_global_indices.size();
                        interior_fract_global_indices.push_back(node_id);
                    }
                }
            }
        }
    }

    // 4. Ensure every fracture node has a global index, without modifying the
    //    interior/border classification already assigned above.
    for (int node_id : all_nodes) {
        if (global_to_interior_fract.find(node_id) == global_to_interior_fract.end() &&
            global_to_border_fract.find(node_id) == global_to_border_fract.end()) {
            // If a node was not encountered in the first pass because it was shared
            // in an unexpected way, default it to interior unless it is a true global border.
            bool is_global_border_node = (border_nodes.find(node_id) != border_nodes.end());
            if (is_global_border_node) {
                global_to_border_fract[node_id] = (int)border_fract_global_indices.size();
                border_fract_global_indices.push_back(node_id);
            } else {
                global_to_interior_fract[node_id] = (int)interior_fract_global_indices.size();
                interior_fract_global_indices.push_back(node_id);
            }
        }
    }

    Nfi = (int)interior_fract_global_indices.size();
    Nfb = (int)border_fract_global_indices.size();

    // Compute bc_map_dfn: classify each DFN border node by its BC type.
    // Frequency-independent: for MIXED_BC only block_index is stored here;
    // alpha is recalculated per electrode in UpdateDDPSCBoundaryConditionsAndSources.
    bc_map_dfn.clear();
    if (!global_border_map.empty()) {
        SourceTermsDef dummy_source;
        SourceTermsCSEigen dummy_epm_source;
        BoundaryConditionsCSEigen bc_dummy = BoundaryConditionsDefinitionCSEigen(
            domain_stored, Nx, Ny, bc_map_classic, dummy_source, dummy_epm_source, true, 0.0);

        for (BordersMap::const_iterator it = global_border_map.begin(); it != global_border_map.end(); ++it) {
            const int index = it->first;
            BoundaryConditionsDef::const_iterator it_classic = bc_map_classic.find(it->second);
            if (it_classic == bc_map_classic.end()) continue;
            const std::string& bc_type = it_classic->second.bc_type;

            if (bc_type == DIRICHLET || bc_type == NEUMANN) {
                bc_map_dfn[index] = it_classic->second;
            } else if (bc_type == DIR_GRADIENT) {
                double border_dist = 0.0, position = 0.0;
                std::string border1, border2;
                position = domain_stored.ReturnBorderCoordinate(global_nodes_map.at(index), border_dist, border1, border2);
                BoundaryConditionsDef::const_iterator it_b1 = bc_map_classic.find(border1);
                BoundaryConditionsDef::const_iterator it_b2 = bc_map_classic.find(border2);
                if (it_b1 != bc_map_classic.end() && it_b2 != bc_map_classic.end() &&
                    it_b1->second.bc_type == DIRICHLET && it_b2->second.bc_type == DIRICHLET) {
                    double bc_value = (it_b2->second.bc_value - it_b1->second.bc_value)
                                      / border_dist * position + it_b1->second.bc_value;
                    bc_map_dfn[index] = BoundaryConditions(DIRICHLET, bc_value);
                }
            } else if (bc_type == MIXED_BC) {
                CgalPoint2D node_pos = global_nodes_map.at(index);
                double x = node_pos.x(), y = node_pos.y();
                if (x == domain_stored.max_pt.i) x -= EPSILON;
                if (y == domain_stored.max_pt.j) y -= EPSILON;
                std::pair<int,int> ij = ::return_indices(x, y, delta_x, delta_y, Nx, Ny);
                BoundaryConditionsCSEigen::const_iterator it_bc = bc_dummy.find(ij);
                if (it_bc != bc_dummy.end()) {
                    const std::map<std::string, BoundaryConditions>& cell_bc = it_bc->second;
                    std::map<std::string, BoundaryConditions>::const_iterator it_cell = cell_bc.find(it->second);
                    if (it_cell != cell_bc.end()) {
                        bc_map_dfn[index] = BoundaryConditions(MIXED_BC, it_cell->second.bc_value);
                        bc_map_dfn[index].block_index = return_index(ij.first, ij.second, Ny);
                    }
                }
            }
        }
    }
    //std::cout << "[PROFILE][InitializeDFNSCIndexing] total: "
    //          << omp_get_wtime() - t_index_total << " sec"
    //          << "  (Nfi=" << Nfi << " Nfb=" << Nfb << ")" << std::endl;
}

// Determine basic blocks / linear system independent parts for DDP condensed system
void DDPModelCS::DetermineBasicLinearSystemDDPCS() {
    //double t_basic_total = omp_get_wtime();

    // 1. Assemble the EPM condensed system blocks (Aii, Aib, Abi, Abb)
    //double t_epm = omp_get_wtime();
    DetermineBasicLinearSystemCS();
    //std::cout << "[PROFILE][DetermineBasicLinearSystemDDPCS] DetermineBasicLinearSystemCS (EPM): "
    //          << omp_get_wtime() - t_epm << " sec" << std::endl;

    // 2. Initialize vectors to collect triplets for the 12 DDP blocks
    std::vector<Triplet> trips_Afi_i2i2, trips_Afi_b2b2, trips_Afi_i2b2, trips_Afi_b2i2;
    std::vector<Triplet> trips_Afm_i2i1, trips_Afm_i2b1, trips_Afm_b2i1, trips_Afm_b2b1;
    std::vector<Triplet> trips_Amf_i1i2, trips_Amf_i1b2, trips_Amf_b1i2, trips_Amf_b1b2;

    trips_Afi_i2i2.reserve(Nfi * 2);
    trips_Afi_b2b2.reserve(Nfb * 2);
    trips_Afi_i2b2.reserve(Nfi * Nfb);
    trips_Afi_b2i2.reserve(Nfb * Nfi);
    trips_Afm_i2i1.reserve(Nfi * 2);
    trips_Afm_i2b1.reserve(Nfi * Nb);
    trips_Afm_b2i1.reserve(Nfb * 2);
    trips_Afm_b2b1.reserve(Nfb * Nb);
    trips_Amf_i1i2.reserve(Ni * 2);
    trips_Amf_i1b2.reserve(Ni * Nfb);
    trips_Amf_b1i2.reserve(Nb * 2);
    trips_Amf_b1b2.reserve(Nb * Nfb);

    // 3. Loop over all subnetworks and assemble fracture-related blocks
    for (SubNetworkMap::iterator it_subnet = subnetwork_map.begin(); 
         it_subnet != subnetwork_map.end(); ++it_subnet) {
        int block_index = it_subnet->first;
        NetworkMeshes &net = it_subnet->second;

        // Get the porous block's local index (interior or border)
        bool block_is_interior = (global_to_interior.count(block_index) > 0);
        int local_block_idx = block_is_interior ? global_to_interior[block_index] : global_to_border[block_index];

        // 3.1. Process each fracture in this subnetwork
           for (std::vector<FractureMesh>::iterator it_fract = net.meshes.begin(); 
               it_fract != net.meshes.end(); ++it_fract) {
              FractureMesh &fract = *it_fract;
            int node1_global = fract.p_ori.index;
            int node2_global = fract.p_tar.index;

            // Determine local indices for the fracture nodes
            bool node1_interior = (global_to_interior_fract.count(node1_global) > 0);
            bool node2_interior = (global_to_interior_fract.count(node2_global) > 0);

            int node1_local = node1_interior ? global_to_interior_fract[node1_global] 
                                              : global_to_border_fract[node1_global];
            int node2_local = node2_interior ? global_to_interior_fract[node2_global] 
                                              : global_to_border_fract[node2_global];

            // Fracture properties
            double coeff = fract.aperture * fract.conductivity;
            double L = fract.ReturnLength();

            // Coefficients for fracture-fracture relationships (basic case, k=0 for now)
            double k_zero = 0.0;
            double a1_ori, a2_ori, a3_ori, a1_tar, a2_tar, a3_tar;
            ReturnDerivativeCoeffFull(0.0, L, k_zero, a1_ori, a2_ori, a3_ori);
            ReturnDerivativeCoeffFull(L, L, k_zero, a1_tar, a2_tar, a3_tar);

            // 3.2. Assemble fracture-fracture blocks (Afi_*)
            // Node 1 - Node 1
            if (node1_interior && node1_interior) {
                trips_Afi_i2i2.emplace_back(node1_local, node1_local, coeff * a1_ori);
            } else if (!node1_interior && !node1_interior) {
                trips_Afi_b2b2.emplace_back(node1_local, node1_local, coeff * a1_ori);
            }

            // Node 1 - Node 2
            if (node1_interior && node2_interior) {
                trips_Afi_i2i2.emplace_back(node1_local, node2_local, coeff * a2_ori);
            } else if (node1_interior && !node2_interior) {
                trips_Afi_i2b2.emplace_back(node1_local, node2_local, coeff * a2_ori);
            } else if (!node1_interior && node2_interior) {
                trips_Afi_b2i2.emplace_back(node1_local, node2_local, coeff * a2_ori);
            } else if (!node1_interior && !node2_interior) {
                trips_Afi_b2b2.emplace_back(node1_local, node2_local, coeff * a2_ori);
            }

            // Node 2 - Node 1
            if (node2_interior && node1_interior) {
                trips_Afi_i2i2.emplace_back(node2_local, node1_local, -coeff * a1_tar);
            } else if (node2_interior && !node1_interior) {
                trips_Afi_i2b2.emplace_back(node2_local, node1_local, -coeff * a1_tar);
            } else if (!node2_interior && node1_interior) {
                trips_Afi_b2i2.emplace_back(node2_local, node1_local, -coeff * a1_tar);
            } else if (!node2_interior && !node1_interior) {
                trips_Afi_b2b2.emplace_back(node2_local, node1_local, -coeff * a1_tar);
            }

            // Node 2 - Node 2
            if (node2_interior && node2_interior) {
                trips_Afi_i2i2.emplace_back(node2_local, node2_local, -coeff * a2_tar);
            } else if (!node2_interior && !node2_interior) {
                trips_Afi_b2b2.emplace_back(node2_local, node2_local, -coeff * a2_tar);
            }

            // 3.3. Assemble fracture-matrix blocks (Afm_* and Amf_*)
            // Node 1 -> Block
            if (node1_interior && block_is_interior) {
                trips_Afm_i2i1.emplace_back(node1_local, local_block_idx, coeff * a3_ori);
            } else if (node1_interior && !block_is_interior) {
                trips_Afm_i2b1.emplace_back(node1_local, local_block_idx, coeff * a3_ori);
            } else if (!node1_interior && block_is_interior) {
                trips_Afm_b2i1.emplace_back(node1_local, local_block_idx, coeff * a3_ori);
            } else if (!node1_interior && !block_is_interior) {
                trips_Afm_b2b1.emplace_back(node1_local, local_block_idx, coeff * a3_ori);
            }

            // Node 2 -> Block
            if (node2_interior && block_is_interior) {
                trips_Afm_i2i1.emplace_back(node2_local, local_block_idx, -coeff * a3_tar);
            } else if (node2_interior && !block_is_interior) {
                trips_Afm_i2b1.emplace_back(node2_local, local_block_idx, -coeff * a3_tar);
            } else if (!node2_interior && block_is_interior) {
                trips_Afm_b2i1.emplace_back(node2_local, local_block_idx, -coeff * a3_tar);
            } else if (!node2_interior && !block_is_interior) {
                trips_Afm_b2b1.emplace_back(node2_local, local_block_idx, -coeff * a3_tar);
            }

            // Block -> Node 1
            if (block_is_interior && node1_interior) {
                trips_Amf_i1i2.emplace_back(local_block_idx, node1_local, coeff * a3_ori);
            } else if (block_is_interior && !node1_interior) {
                trips_Amf_i1b2.emplace_back(local_block_idx, node1_local, coeff * a3_ori);
            } else if (!block_is_interior && node1_interior) {
                trips_Amf_b1i2.emplace_back(local_block_idx, node1_local, coeff * a3_ori);
            } else if (!block_is_interior && !node1_interior) {
                trips_Amf_b1b2.emplace_back(local_block_idx, node1_local, coeff * a3_ori);
            }

            // Block -> Node 2
            if (block_is_interior && node2_interior) {
                trips_Amf_i1i2.emplace_back(local_block_idx, node2_local, -coeff * a3_tar);
            } else if (block_is_interior && !node2_interior) {
                trips_Amf_i1b2.emplace_back(local_block_idx, node2_local, -coeff * a3_tar);
            } else if (!block_is_interior && node2_interior) {
                trips_Amf_b1i2.emplace_back(local_block_idx, node2_local, -coeff * a3_tar);
            } else if (!block_is_interior && !node2_interior) {
                trips_Amf_b1b2.emplace_back(local_block_idx, node2_local, -coeff * a3_tar);
            }
        }
    }

    //std::cout << "[PROFILE][DetermineBasicLinearSystemDDPCS] fracture loop (triplet assembly): "
    //          << omp_get_wtime() - t_basic_total << " sec" << std::endl;

    // 4. Resize and assemble all 12 blocks
    //double t_assemble = omp_get_wtime();
    block_Afi_i2i2.resize(Nfi, Nfi);
    block_Afi_b2b2.resize(Nfb, Nfb);
    block_Afi_i2b2.resize(Nfi, Nfb);
    block_Afi_b2i2.resize(Nfb, Nfi);

    block_Afm_i2i1.resize(Nfi, Ni);
    block_Afm_i2b1.resize(Nfi, Nb);
    block_Afm_b2i1.resize(Nfb, Ni);
    block_Afm_b2b1.resize(Nfb, Nb);

    block_Amf_i1i2.resize(Ni, Nfi);
    block_Amf_i1b2.resize(Ni, Nfb);
    block_Amf_b1i2.resize(Nb, Nfi);
    block_Amf_b1b2.resize(Nb, Nfb);

    block_Afi_i2i2.setFromTriplets(trips_Afi_i2i2.begin(), trips_Afi_i2i2.end());
    block_Afi_b2b2.setFromTriplets(trips_Afi_b2b2.begin(), trips_Afi_b2b2.end());
    block_Afi_i2b2.setFromTriplets(trips_Afi_i2b2.begin(), trips_Afi_i2b2.end());
    block_Afi_b2i2.setFromTriplets(trips_Afi_b2i2.begin(), trips_Afi_b2i2.end());

    block_Afm_i2i1.setFromTriplets(trips_Afm_i2i1.begin(), trips_Afm_i2i1.end());
    block_Afm_i2b1.setFromTriplets(trips_Afm_i2b1.begin(), trips_Afm_i2b1.end());
    block_Afm_b2i1.setFromTriplets(trips_Afm_b2i1.begin(), trips_Afm_b2i1.end());
    block_Afm_b2b1.setFromTriplets(trips_Afm_b2b1.begin(), trips_Afm_b2b1.end());

    block_Amf_i1i2.setFromTriplets(trips_Amf_i1i2.begin(), trips_Amf_i1i2.end());
    block_Amf_i1b2.setFromTriplets(trips_Amf_i1b2.begin(), trips_Amf_i1b2.end());
    block_Amf_b1i2.setFromTriplets(trips_Amf_b1i2.begin(), trips_Amf_b1i2.end());
    block_Amf_b1b2.setFromTriplets(trips_Amf_b1b2.begin(), trips_Amf_b1b2.end());

    // 5. Compress matrices for efficient storage
    block_Afi_i2i2.makeCompressed();
    block_Afi_b2b2.makeCompressed();
    block_Afi_i2b2.makeCompressed();
    block_Afi_b2i2.makeCompressed();
    block_Afm_i2i1.makeCompressed();
    block_Afm_i2b1.makeCompressed();
    block_Afm_b2i1.makeCompressed();
    block_Afm_b2b1.makeCompressed();
    block_Amf_i1i2.makeCompressed();
    block_Amf_i1b2.makeCompressed();
    block_Amf_b1i2.makeCompressed();
    block_Amf_b1b2.makeCompressed();

    //std::cout << "[PROFILE][DetermineBasicLinearSystemDDPCS] setFromTriplets + makeCompressed (12 blocs): "
    //          << omp_get_wtime() - t_assemble << " sec" << std::endl;
    //std::cout << "[PROFILE][DetermineBasicLinearSystemDDPCS] total: "
    //          << omp_get_wtime() - t_basic_total << " sec" << std::endl;
}

void DDPModelCS::UpdateDDPSCBlocksFourier() {
    // 1. Update EPM matrix blocks that depend on Fourier variable
    UpdateSCBlocksFourier();

    // 2. Prepare triplet lists for fracture-related blocks
    std::vector<Triplet> trips_Afi_i2i2, trips_Afi_b2b2, trips_Afi_i2b2, trips_Afi_b2i2;
    std::vector<Triplet> trips_Afm_i2i1, trips_Afm_i2b1, trips_Afm_b2i1, trips_Afm_b2b1;
    std::vector<Triplet> trips_Amf_i1i2, trips_Amf_i1b2, trips_Amf_b1i2, trips_Amf_b1b2;

    trips_Afi_i2i2.reserve(Nfi * 4);
    trips_Afi_b2b2.reserve(Nfb * 4);
    trips_Afi_i2b2.reserve(std::max(1, Nfi) * std::max(1, Nfb));
    trips_Afi_b2i2.reserve(std::max(1, Nfb) * std::max(1, Nfi));

    // 3. Loop on subnetworks and assemble Fourier-dependent contributions
    for (SubNetworkMap::iterator it_sub = subnetwork_map.begin(); it_sub != subnetwork_map.end(); ++it_sub) {
        int block_index = it_sub->first;
        NetworkMeshes &net = it_sub->second;

        // determine porous block index (interior or border in EPM indexing)
        bool block_is_interior = (global_to_interior.count(block_index) > 0);
        int local_block_idx = block_is_interior ? global_to_interior[block_index] : global_to_border[block_index];

        // get porous conductivity for this block
        std::pair<int,int> ij = return_indices(block_index, Ny);
        double porous_cond_block = porous_cond(ij.first, ij.second);

        // compute alpha for this block using ExchangeFunctions
        double alpha_block = ReturnAlphaBlock(porous_cond_block, delta_x, delta_y, net);

        // process fractures in subnetwork
        for (std::vector<FractureMesh>::iterator itf = net.meshes.begin(); itf != net.meshes.end(); ++itf) {
            FractureMesh &fract = *itf;
            int g1 = fract.p_ori.index;
            int g2 = fract.p_tar.index;

            bool n1_in = (global_to_interior_fract.count(g1) > 0);
            bool n2_in = (global_to_interior_fract.count(g2) > 0);
            int l1 = n1_in ? global_to_interior_fract[g1] : global_to_border_fract[g1];
            int l2 = n2_in ? global_to_interior_fract[g2] : global_to_border_fract[g2];

            double coeff = fract.aperture * fract.conductivity;
            double L = fract.ReturnLength();

            // k_fract for Fourier case (matches existing DFNDDPLinearSystemFourier)
            double k_fract = 0.0;
            if (fract.conductivity * fract.aperture != 0.0) {
                k_fract = alpha_block / (fract.conductivity * fract.aperture);
            }

            double a1_0, a2_0, a3_0;
            double a1_L, a2_L, a3_L;
            ReturnDerivativeCoeffFullFourier(0.0, L, k_fract, current_w, a1_0, a2_0, a3_0);
            ReturnDerivativeCoeffFullFourier(L, L, k_fract, current_w, a1_L, a2_L, a3_L);

            // node1 diagonal: mass-like term minus derivative part
            double diag1 = coeff * fract.aperture * current_w * current_w - coeff * a1_0;
            if (n1_in) trips_Afi_i2i2.emplace_back(l1, l1, diag1);
            else trips_Afi_b2b2.emplace_back(l1, l1, diag1);

            // node1 -> node2 coupling
            double c12 = -coeff * a2_0;
            if (n1_in && n2_in) trips_Afi_i2i2.emplace_back(l1, l2, c12);
            else if (n1_in && !n2_in) trips_Afi_i2b2.emplace_back(l1, l2, c12);
            else if (!n1_in && n2_in) trips_Afi_b2i2.emplace_back(l1, l2, c12);
            else trips_Afi_b2b2.emplace_back(l1, l2, c12);

            // node2 diagonal and couplings
            double diag2 = coeff * a2_L; // note: DFNDDPLinearSystemFourier mixes signs; we follow same pattern
            if (n2_in) trips_Afi_i2i2.emplace_back(l2, l2, diag2);
            else trips_Afi_b2b2.emplace_back(l2, l2, diag2);

            double c21 = coeff * a1_L;
            if (n2_in && n1_in) trips_Afi_i2i2.emplace_back(l2, l1, c21);
            else if (n2_in && !n1_in) trips_Afi_i2b2.emplace_back(l2, l1, c21);
            else if (!n2_in && n1_in) trips_Afi_b2i2.emplace_back(l2, l1, c21);
            else trips_Afi_b2b2.emplace_back(l2, l1, c21);

            // Fracture-matrix coupling via a3 coefficients
            double fm1 = -coeff * a3_0;
            double fm2 = coeff * a3_L;
            if (n1_in) {
                if (block_is_interior) trips_Afm_i2i1.emplace_back(l1, local_block_idx, fm1);
                else trips_Afm_i2b1.emplace_back(l1, local_block_idx, fm1);
            } else {
                if (block_is_interior) trips_Afm_b2i1.emplace_back(l1, local_block_idx, fm1);
                else trips_Afm_b2b1.emplace_back(l1, local_block_idx, fm1);
            }

            if (n2_in) {
                if (block_is_interior) trips_Afm_i2i1.emplace_back(l2, local_block_idx, fm2);
                else trips_Afm_i2b1.emplace_back(l2, local_block_idx, fm2);
            } else {
                if (block_is_interior) trips_Afm_b2i1.emplace_back(l2, local_block_idx, fm2);
                else trips_Afm_b2b1.emplace_back(l2, local_block_idx, fm2);
            }

            // Matrix -> fracture (transpose-like)
            if (block_is_interior) {
                if (n1_in) trips_Amf_i1i2.emplace_back(local_block_idx, l1, fm1);
                else trips_Amf_i1b2.emplace_back(local_block_idx, l1, fm1);
                if (n2_in) trips_Amf_i1i2.emplace_back(local_block_idx, l2, fm2);
                else trips_Amf_i1b2.emplace_back(local_block_idx, l2, fm2);
            } else {
                if (n1_in) trips_Amf_b1i2.emplace_back(local_block_idx, l1, fm1);
                else trips_Amf_b1b2.emplace_back(local_block_idx, l1, fm1);
                if (n2_in) trips_Amf_b1i2.emplace_back(local_block_idx, l2, fm2);
                else trips_Amf_b1b2.emplace_back(local_block_idx, l2, fm2);
            }

            // Diagonal exchange: alpha_block * (L - coeff3) on the matrix block diagonal.
            // In DDP (PorousDDPLinearSystemUpdateFourier), the matrix equation contains the
            // integral of phi_f along the fracture, which contributes a self-coupling term
            // alpha*(L - coeff3)*phi_m. This term is absent from the Amf transpose and must
            // be added explicitly to block_Aii / block_Abb.
            if (alpha_block != 0.0 && (k_fract > 0.0 || current_w > 0.0)) {
                double c1, c2, c3;
                CumulatePotFractCoeffFourier(L, k_fract, current_w, c1, c2, c3);
                const double diag_exchange = alpha_block * (L - c3);
                if (block_is_interior) {
                    block_Aii.coeffRef(local_block_idx, local_block_idx) += diag_exchange;
                } else {
                    block_Abb.coeffRef(local_block_idx, local_block_idx) += diag_exchange;
                }
            }
        }
    }

    // 4. Assign to sparse matrices and compress
    block_Afi_i2i2.resize(Nfi, Nfi);
    block_Afi_b2b2.resize(Nfb, Nfb);
    block_Afi_i2b2.resize(Nfi, Nfb);
    block_Afi_b2i2.resize(Nfb, Nfi);

    block_Afm_i2i1.resize(Nfi, Ni);
    block_Afm_i2b1.resize(Nfi, Nb);
    block_Afm_b2i1.resize(Nfb, Ni);
    block_Afm_b2b1.resize(Nfb, Nb);

    block_Amf_i1i2.resize(Ni, Nfi);
    block_Amf_i1b2.resize(Ni, Nfb);
    block_Amf_b1i2.resize(Nb, Nfi);
    block_Amf_b1b2.resize(Nb, Nfb);

    block_Afi_i2i2.setFromTriplets(trips_Afi_i2i2.begin(), trips_Afi_i2i2.end());
    block_Afi_b2b2.setFromTriplets(trips_Afi_b2b2.begin(), trips_Afi_b2b2.end());
    block_Afi_i2b2.setFromTriplets(trips_Afi_i2b2.begin(), trips_Afi_i2b2.end());
    block_Afi_b2i2.setFromTriplets(trips_Afi_b2i2.begin(), trips_Afi_b2i2.end());

    block_Afm_i2i1.setFromTriplets(trips_Afm_i2i1.begin(), trips_Afm_i2i1.end());
    block_Afm_i2b1.setFromTriplets(trips_Afm_i2b1.begin(), trips_Afm_i2b1.end());
    block_Afm_b2i1.setFromTriplets(trips_Afm_b2i1.begin(), trips_Afm_b2i1.end());
    block_Afm_b2b1.setFromTriplets(trips_Afm_b2b1.begin(), trips_Afm_b2b1.end());

    block_Amf_i1i2.setFromTriplets(trips_Amf_i1i2.begin(), trips_Amf_i1i2.end());
    block_Amf_i1b2.setFromTriplets(trips_Amf_i1b2.begin(), trips_Amf_i1b2.end());
    block_Amf_b1i2.setFromTriplets(trips_Amf_b1i2.begin(), trips_Amf_b1i2.end());
    block_Amf_b1b2.setFromTriplets(trips_Amf_b1b2.begin(), trips_Amf_b1b2.end());

    block_Afi_i2i2.makeCompressed();
    block_Afi_b2b2.makeCompressed();
    block_Afi_i2b2.makeCompressed();
    block_Afi_b2i2.makeCompressed();
    block_Afm_i2i1.makeCompressed();
    block_Afm_i2b1.makeCompressed();
    block_Afm_b2i1.makeCompressed();
    block_Afm_b2b1.makeCompressed();
    block_Amf_i1i2.makeCompressed();
    block_Amf_i1b2.makeCompressed();
    block_Amf_b1i2.makeCompressed();
    block_Amf_b1b2.makeCompressed();
}

void DDPModelCS::UpdateDDPSCBoundaryConditionsAndSources(bool fourier) {
    clock_t t_total = clock();

    // Boundary-condition update for the condensed DDP border system.
    //
    // Border ordering:
    //   [phi_fb ; phi_mb]
    // where phi_fb are fracture border unknowns and phi_mb are matrix border unknowns.

    if (matrix_S_DDP_precomputed.rows() == 0 || matrix_S_DDP_precomputed.cols() == 0) {
        matrix_S_DDP.resize(0, 0);
        vector_rhs_DDP = Vector::Zero(0);
        return;
    }

    clock_t t_copy = clock();
    matrix_S_DDP = matrix_S_DDP_precomputed;
    // Start RHS from the precomputed condensed contribution (b_b - Abi_g * z_i),
    // then BC corrections are added on top — same pattern as EPM UpdateSCBoundaryConditionsAndSources.
    vector_rhs_DDP = vector_zb_DDP;
    vector_b_fi = Vector::Zero(Nfi);
    vector_b_fb = Vector::Zero(Nfb);
    vector_delta_diag_DDP = Vector::Zero(Nb);
    // std::cout << "[PROFILE][DDPCS]     Step 1 - Copy precomputed matrix_S_DDP: "
    //           << (double)(clock() - t_copy) / CLOCKS_PER_SEC << " sec" << std::endl;

    // ------------------------------------------------------------------
    // 1) Matrix border unknowns: same spirit as EPMSystemCSEigen.
    //    We add the diagonal correction and the source/BC contribution
    //    on the matrix-border part of the condensed DDP system.
    // ------------------------------------------------------------------
    clock_t t_matrix_bc = clock();
    for (int lb = 0; lb < Nb; ++lb) {
        const int g = border_global_indices[lb];
        const std::pair<int, int> ij = return_indices(g, Ny);
        const int i = ij.first;
        const int j = ij.second;

        std::map<std::string, BoundaryConditions> bc_cell = bc_def_epm[std::make_pair(i, j)];

        if (i == 0) {
            double coeff1 = 0.0, coeff3 = 0.0;
            BoundCondCoeffsCSEigen(LEFT_BORDER, bc_cell[LEFT_BORDER], porous_cond(i, j), delta_x, delta_y, coeff1, coeff3);
            vector_delta_diag_DDP(lb) += coeff1;
            vector_rhs_DDP(Nfb + lb) += coeff3;
        }
        if (i == Nx - 1) {
            double coeff1 = 0.0, coeff3 = 0.0;
            BoundCondCoeffsCSEigen(RIGHT_BORDER, bc_cell[RIGHT_BORDER], porous_cond(i, j), delta_x, delta_y, coeff1, coeff3);
            vector_delta_diag_DDP(lb) += coeff1;
            vector_rhs_DDP(Nfb + lb) += coeff3;
        }
        if (j == 0) {
            double coeff1 = 0.0, coeff3 = 0.0;
            BoundCondCoeffsCSEigen(BOTTOM_BORDER, bc_cell[BOTTOM_BORDER], porous_cond(i, j), delta_x, delta_y, coeff1, coeff3);
            vector_delta_diag_DDP(lb) += coeff1;
            vector_rhs_DDP(Nfb + lb) += coeff3;
        }
        if (j == Ny - 1) {
            double coeff1 = 0.0, coeff3 = 0.0;
            BoundCondCoeffsCSEigen(TOP_BORDER, bc_cell[TOP_BORDER], porous_cond(i, j), delta_x, delta_y, coeff1, coeff3);
            vector_delta_diag_DDP(lb) += coeff1;
            vector_rhs_DDP(Nfb + lb) += coeff3;
        }

        // Source terms on matrix border cells (same philosophy as EPM update).
        for (SourceTermsCSEigen::iterator it = source_terms_epm.begin(); it != source_terms_epm.end(); ++it) {
            if (return_index(it->first.first, it->first.second, Ny) == g) {
                double val = it->second;
                if (fourier) {
                    val *= 0.5;
                }
                vector_rhs_DDP(Nfb + lb) += val;
            }
        }
    }

    // Inject the matrix-border diagonal correction into the condensed matrix.
    for (int lb = 0; lb < Nb; ++lb) {
        matrix_S_DDP.coeffRef(Nfb + lb, Nfb + lb) += vector_delta_diag_DDP(lb);
    }

    // ------------------------------------------------------------------
    // Interior matrix source terms: these are NOT in the border loop above.
    // They must be propagated to the condensed RHS via the Schur formula:
    //   delta_zi  = Aii_g^{-1} * delta_rhs_i
    //   rhs_b    -= Abi_g * delta_zi
    // This matches the EPM UpdateSCBoundaryConditionsAndSources interior logic.
    // ------------------------------------------------------------------
    {
        const int N_int_global = Nfi + Ni;
        const int N_bor_global = Nfb + Nb;
        Vector delta_rhs_i = Vector::Zero(N_int_global);
        bool has_interior_source = false;
        for (SourceTermsCSEigen::const_iterator it = source_terms_epm.begin();
             it != source_terms_epm.end(); ++it) {
            const int g = return_index(it->first.first, it->first.second, Ny);
            if (global_to_interior.count(g) > 0) {
                double val = it->second;
                if (fourier) val *= 0.5;
                delta_rhs_i(Nfi + global_to_interior.at(g)) += val;
                has_interior_source = true;
            }
        }
        if (has_interior_source && N_int_global > 0) {
            if (!factor_Afi_i2i2_ready) {
                factor_Ai1plusi2 = std::make_shared<LinearSystemSparseLLTFactorization>();
                LinearSystemFactorizeLLT(Aii_g, *factor_Ai1plusi2);
                factor_Afi_i2i2_ready = true;
            }
            Vector delta_zi = LinearSystemSolveFactorizedLLT(*factor_Ai1plusi2, delta_rhs_i);
            vector_zi_DDP += delta_zi;
            if (N_bor_global > 0) {
                vector_rhs_DDP.noalias() -= (Abi_g * delta_zi);
            }
        }
    }

    // std::cout << "[PROFILE][DDPCS]     Step 2 - Matrix border BC (delta_diag + sources): "
    //           << (double)(clock() - t_matrix_bc) / CLOCKS_PER_SEC << " sec" << std::endl;

    // ------------------------------------------------------------------
    // 2) Fracture border unknowns.
    //    Apply the DFN boundary conditions directly on the first block.
    // ------------------------------------------------------------------
    clock_t t_frac_bc = clock();
    for (BoundaryConditionsDFN::iterator it = bc_map_dfn.begin(); it != bc_map_dfn.end(); ++it) {
        const int node = it->first;
        const BoundaryConditions& bc = it->second;

        if (global_to_border_fract.count(node) == 0) {
            continue;
        }

        const int row = global_to_border_fract[node];

        if (bc.bc_type == DIRICHLET) {
            // Symmetric elimination: move column contribution to RHS before zeroing.
            for (int r = 0; r < matrix_S_DDP.rows(); ++r) {
                if (r != row) {
                    vector_rhs_DDP(r) -= matrix_S_DDP.coeffRef(r, row) * bc.bc_value;
                    matrix_S_DDP.coeffRef(r, row) = 0.0;
                }
            }
            matrix_S_DDP.row(row).setZero();
            matrix_S_DDP.coeffRef(row, row) = 1.0;
            vector_rhs_DDP(row) = bc.bc_value;
            vector_b_fb(row) = bc.bc_value;
        } else if (bc.bc_type == NEUMANN) {
            vector_rhs_DDP(row) += bc.bc_value;
            vector_b_fb(row) += bc.bc_value;
        } else if (bc.bc_type == MIXED_BC) {
            if (fourier) {
                int block_index = bc.block_index;
                if (block_index < 0) {
                    block_index = ReturnBlockIndex(subnetwork_map, node);
                }
                if (block_index >= 0 && subnetwork_map.count(block_index) > 0) {
                    NetworkMeshes& net = subnetwork_map[block_index];
                    double alpha = ::ReturnAlphaBCMixed(node,
                                                        block_index,
                                                        delta_x,
                                                        delta_y,
                                                        current_w,
                                                        source_terms_epm,
                                                        subnetwork_map);
                    double aper = net.return_aperture_from_node(node);
                    double cond = net.return_conductivity_from_node(node);
                    double corr = 2.0 * aper * cond * alpha / (2.0 + alpha * delta_x);
                    matrix_S_DDP.coeffRef(row, row) += corr;
                }
            }
        } else {
            std::cout << "WARNING in DDPModelCS::UpdateDDPSCBoundaryConditionsAndSources: fracture boundary condition not implemented" << std::endl;
        }
    }
    // std::cout << "[PROFILE][DDPCS]     Step 3 - Fracture border BC: "
    //           << (double)(clock() - t_frac_bc) / CLOCKS_PER_SEC << " sec" << std::endl;

    // std::cout << "[PROFILE][DDPCS]   Total UpdateDDPSCBoundaryConditionsAndSources: "
    //           << (double)(clock() - t_total) / CLOCKS_PER_SEC << " sec" << std::endl;
}

void DDPModelCS::PrecomputeDDPSCForCurrentFrequency() {
    clock_t t_total = clock();

    // 1. Ensure indexing and basic blocks are ready.
    clock_t t_init = clock();
    if (!IsInitializedDFNSC()) {
        InitializeDFNSCIndexing();
    }
    if (block_Afi_i2i2.rows() == 0 || block_Aii.rows() == 0) {
        DetermineBasicLinearSystemDDPCS();
    }
    // std::cout << "[PROFILE][DDPCS]     Step 1 - Init/indexing check: "
    //           << (double)(clock() - t_init) / CLOCKS_PER_SEC << " sec" << std::endl;

    // 2. Build the 4 global blocks used by the static condensation stage.
    //    Global ordering:
    //      x_i = [phi_fi ; phi_mi]
    //      x_b = [phi_fb ; phi_mb]
    const int N_int_global = Nfi + Ni;
    const int N_bor_global = Nfb + Nb;

    clock_t t_blocks = clock();

    std::vector<Triplet> trips_Aii_g;
    std::vector<Triplet> trips_Aib_g;
    std::vector<Triplet> trips_Abi_g;
    std::vector<Triplet> trips_Abb_g;

    trips_Aii_g.reserve(block_Afi_i2i2.nonZeros() + block_Afm_i2i1.nonZeros() +
                        block_Amf_i1i2.nonZeros() + block_Aii.nonZeros());
    trips_Aib_g.reserve(block_Afi_i2b2.nonZeros() + block_Afm_i2b1.nonZeros() +
                        block_Amf_i1b2.nonZeros() + block_Aib.nonZeros());
    trips_Abi_g.reserve(block_Afi_b2i2.nonZeros() + block_Afm_b2i1.nonZeros() +
                        block_Amf_b1i2.nonZeros() + block_Abi.nonZeros());
    trips_Abb_g.reserve(block_Afi_b2b2.nonZeros() + block_Afm_b2b1.nonZeros() +
                        block_Amf_b1b2.nonZeros() + block_Abb.nonZeros());

    const auto append_block = [](const SparseMatrix& src,
                                 int row_offset,
                                 int col_offset,
                                 std::vector<Triplet>& trips) {
        for (int k = 0; k < src.outerSize(); ++k) {
            for (SparseMatrix::InnerIterator it(src, k); it; ++it) {
                trips.emplace_back(row_offset + it.row(), col_offset + it.col(), it.value());
            }
        }
    };

    // Aii_g = [A_fi_fi  A_fi_mi; A_mi_fi  A_mi_mi]
    append_block(block_Afi_i2i2, 0, 0, trips_Aii_g);
    append_block(block_Afm_i2i1, 0, Nfi, trips_Aii_g);
    append_block(block_Amf_i1i2, Nfi, 0, trips_Aii_g);
    append_block(block_Aii, Nfi, Nfi, trips_Aii_g);

    // Aib_g = [A_fi_fb  A_fi_mb; A_mi_fb  A_mi_mb]
    append_block(block_Afi_i2b2, 0, 0, trips_Aib_g);
    append_block(block_Afm_i2b1, 0, Nfb, trips_Aib_g);
    append_block(block_Amf_i1b2, Nfi, 0, trips_Aib_g);
    append_block(block_Aib, Nfi, Nfb, trips_Aib_g);

    // Abi_g = [A_fb_fi  A_fb_mi; A_mb_fi  A_mb_mi]
    append_block(block_Afi_b2i2, 0, 0, trips_Abi_g);
    append_block(block_Afm_b2i1, 0, Nfi, trips_Abi_g);
    append_block(block_Amf_b1i2, Nfb, 0, trips_Abi_g);
    append_block(block_Abi, Nfb, Nfi, trips_Abi_g);

    // Abb_g = [A_fb_fb  A_fb_mb; A_mb_fb  A_mb_mb]
    append_block(block_Afi_b2b2, 0, 0, trips_Abb_g);
    append_block(block_Afm_b2b1, 0, Nfb, trips_Abb_g);
    append_block(block_Amf_b1b2, Nfb, 0, trips_Abb_g);
    append_block(block_Abb, Nfb, Nfb, trips_Abb_g);

    Aii_g.resize(N_int_global, N_int_global);
    Aib_g.resize(N_int_global, N_bor_global);
    Abi_g.resize(N_bor_global, N_int_global);
    Abb_g.resize(N_bor_global, N_bor_global);

    Aii_g.setFromTriplets(trips_Aii_g.begin(), trips_Aii_g.end());
    Aib_g.setFromTriplets(trips_Aib_g.begin(), trips_Aib_g.end());
    Abi_g.setFromTriplets(trips_Abi_g.begin(), trips_Abi_g.end());
    Abb_g.setFromTriplets(trips_Abb_g.begin(), trips_Abb_g.end());

    Aii_g.makeCompressed();
    Aib_g.makeCompressed();
    Abi_g.makeCompressed();
    Abb_g.makeCompressed();

    // std::cout << "[PROFILE][DDPCS]     Step 2 - Build global blocks (Aii_g, Aib_g, Abi_g, Abb_g): "
    //           << (double)(clock() - t_blocks) / CLOCKS_PER_SEC << " sec" << std::endl;

    // 3. Static condensation precomputation.
    //    Compute Y = Aii_g^{-1} Aib_g and the condensed matrix
    //    S = Abb_g - Abi_g * Y.
    clock_t t_factorize = clock();
    if (N_int_global > 0) {
        factor_Ai1plusi2 = std::make_shared<LinearSystemSparseLLTFactorization>();
        LinearSystemFactorizeLLT(Aii_g, *factor_Ai1plusi2);
        factor_Afi_i2i2_ready = true;
    }
    // std::cout << "[PROFILE][DDPCS]     Step 3 - Factorize Aii_g: "
            //   << (double)(clock() - t_factorize) / CLOCKS_PER_SEC << " sec" << std::endl;

    clock_t t_matrix_y = clock();
    matrix_Y_DDP = DenseMatrix::Zero(N_int_global, N_bor_global);
    if (N_int_global > 0 && N_bor_global > 0) {
        // SimplicialLLT::solve() est const → thread-safe pour des RHS différents.
        // schedule(static) + vecteurs thread-locaux pré-alloués pour éviter la contention malloc.
        // Note : si le facteur LLT est très grand (>L3 cache), réduire nb_threads peut être plus rapide.
        #pragma omp parallel num_threads(14)
        {
            Vector rhs_col(N_int_global);
            Vector x_col(N_int_global);
            #pragma omp for schedule(static)
            for (int c = 0; c < N_bor_global; ++c) {
                rhs_col.setZero();
                for (SparseMatrix::InnerIterator it(Aib_g, c); it; ++it) {
                    rhs_col(it.row()) = it.value();
                }
                x_col = LinearSystemSolveFactorizedLLT(*factor_Ai1plusi2, rhs_col);
                matrix_Y_DDP.col(c) = x_col;
            }
        }
    }
    // std::cout << "[PROFILE][DDPCS]     Step 4 - Compute matrix_Y_DDP: "
    //           << (double)(clock() - t_matrix_y) / CLOCKS_PER_SEC << " sec" << std::endl;

    clock_t t_matrix_s = clock();
    matrix_S_DDP = DenseMatrix::Zero(N_bor_global, N_bor_global);
    for (int j = 0; j < N_bor_global; ++j) {
        for (SparseMatrix::InnerIterator it(Abb_g, j); it; ++it) {
            matrix_S_DDP(it.row(), it.col()) += it.value();
        }
    }
    if (N_bor_global > 0 && N_int_global > 0) {
        matrix_S_DDP.noalias() -= (Abi_g * matrix_Y_DDP);
    }
    matrix_S_DDP_precomputed = matrix_S_DDP;
    matrix_Y_DDP = DenseMatrix();  // ~266 MB freed here — not needed after this point
    // std::cout << "[PROFILE][DDPCS]     Step 5 - Compute matrix_S_DDP: "
    //           << (double)(clock() - t_matrix_s) / CLOCKS_PER_SEC << " sec" << std::endl;

    // 4. Build split RHS vectors and store the intermediate condensed vectors.
    clock_t t_rhs = clock();
    Vector rhs_i = Vector::Zero(N_int_global);

    if (Nfi > 0 && vector_b_fi.size() == Nfi) {
        rhs_i.head(Nfi) = vector_b_fi;
    }
    if (Ni > 0 && vector_bi.size() == Ni) {
        rhs_i.tail(Ni) = vector_bi;
    }

    vector_zi_DDP = (N_int_global > 0) ? LinearSystemSolveFactorizedLLT(*factor_Ai1plusi2, rhs_i) : Vector::Zero(0);

    // rhs_b = [b_fb ; b_mb] — initialize before subtracting the interior contribution
    vector_zb_DDP = Vector::Zero(N_bor_global);
    if (Nfb > 0 && vector_b_fb.size() == Nfb) {
        vector_zb_DDP.head(Nfb) = vector_b_fb;
    }
    if (Nb > 0 && vector_bb.size() == Nb) {
        vector_zb_DDP.tail(Nb) = vector_bb;
    }
    if (N_bor_global > 0 && N_int_global > 0) {
        vector_zb_DDP.noalias() -= (Abi_g * vector_zi_DDP);
    }
    vector_rhs_DDP = vector_zb_DDP;
    // std::cout << "[PROFILE][DDPCS]     Step 6 - Compute vector_zi_DDP / vector_zb_DDP: "
    //           << (double)(clock() - t_rhs) / CLOCKS_PER_SEC << " sec" << std::endl;

    // std::cout << "[PROFILE][DDPCS]   Total PrecomputeDDPSCForCurrentFrequency: "
    //           << (double)(clock() - t_total) / CLOCKS_PER_SEC << " sec" << std::endl;
}

void DDPModelCS::PrecomputeDDPSCForCurrentFrequency2() {
    // Nested Schur condensation on the 16 native blocks.
    //
    // Ordering:
    //   x_i = [phi_fi ; phi_mi]
    //   x_b = [phi_fb ; phi_mb]
    //
    // We eliminate phi_fi first, then phi_mi, without building Aii_g/Aib_g/Abi_g/Abb_g.

    if (!IsInitializedDFNSC()) {
        InitializeDFNSCIndexing();
    }

    if (block_Afi_i2i2.rows() == 0 || block_Aii.rows() == 0) {
        DetermineBasicLinearSystemDDPCS();
    }

    // Safety fallback for degenerate cases.
    if (Nfi == 0 || Ni == 0) {
        PrecomputeDDPSCForCurrentFrequency();
        return;
    }

    const int N_int_global = Nfi + Ni;
    const int N_bor_global = Nfb + Nb;

    // Local aliases for readability.
    const SparseMatrix& A = block_Afi_i2i2; // A_fi_fi
    const SparseMatrix& B = block_Afm_i2i1;  // A_fi_mi
    const SparseMatrix& C = block_Afi_i2b2;  // A_fi_fb
    const SparseMatrix& D = block_Afm_i2b1;  // A_fi_mb

    const SparseMatrix& E = block_Amf_i1i2;  // A_mi_fi
    const SparseMatrix& F = block_Aii;       // A_mi_mi
    const SparseMatrix& Gblk = block_Amf_i1b2;  // A_mi_fb
    const SparseMatrix& H = block_Aib;       // A_mi_mb

    const SparseMatrix& K = block_Afi_b2i2;  // A_fb_fi
    const SparseMatrix& L = block_Afm_b2i1;   // A_fb_mi
    const SparseMatrix& M = block_Afi_b2b2;   // A_fb_fb
    const SparseMatrix& N = block_Afm_b2b1;   // A_fb_mb

    const SparseMatrix& P = block_Amf_b1i2;   // A_mb_fi
    const SparseMatrix& Q = block_Abi;        // A_mb_mi
    const SparseMatrix& R = block_Amf_b1b2;   // A_mb_fb
    const SparseMatrix& T = block_Abb;        // A_mb_mb

    // ---------------------------------------------------------------------
    // 1) Factorize A = A_fi_fi
    // ---------------------------------------------------------------------
    factor_Afi_i2i2_ready = false;
    if (A.rows() > 0) {
        factor_Ai1plusi2 = std::make_shared<LinearSystemSparseLLTFactorization>();
        LinearSystemFactorizeLLT(A, *factor_Ai1plusi2);
        factor_Afi_i2i2_ready = true;
    }

    const int nA = Nfi;
    const int nB = Ni;
    const int nC = Nfb;
    const int nD = Nb;

    // Helper lambdas for sparse -> dense extraction.
    const auto sparse_to_dense = [](const SparseMatrix& src, int rows, int cols) -> DenseMatrix {
        DenseMatrix dense = DenseMatrix::Zero(rows, cols);
        for (int c = 0; c < src.outerSize(); ++c) {
            for (SparseMatrix::InnerIterator it(src, c); it; ++it) {
                dense(it.row(), it.col()) += it.value();
            }
        }
        return dense;
    };

    const auto solve_A = [&](const DenseMatrix& rhs) -> DenseMatrix {
        DenseMatrix out = DenseMatrix::Zero(rhs.rows(), rhs.cols());
        if (nA == 0 || rhs.cols() == 0) {
            return out;
        }
        Vector rhs_col = Vector::Zero(nA);
        for (int c = 0; c < rhs.cols(); ++c) {
            for (int r = 0; r < nA; ++r) {
                rhs_col(r) = rhs(r, c);
            }
            Vector x = LinearSystemSolveFactorizedLLT(*factor_Ai1plusi2, rhs_col);
            for (int r = 0; r < nA; ++r) {
                out(r, c) = x(r);
            }
        }
        return out;
    };

    // ---------------------------------------------------------------------
    // 2) Compute Ainv_bfi, Ainv_B, Ainv_C, Ainv_D
    // ---------------------------------------------------------------------
    DenseMatrix B_dense = sparse_to_dense(B, nA, nB);
    DenseMatrix C_dense = sparse_to_dense(C, nA, nC);
    DenseMatrix D_dense = sparse_to_dense(D, nA, nD);

    Vector bfi = Vector::Zero(nA);
    if (vector_b_fi.size() == nA) {
        bfi = vector_b_fi;
    }

    DenseMatrix Ainv_B = solve_A(B_dense);
    DenseMatrix Ainv_C = solve_A(C_dense);
    DenseMatrix Ainv_D = solve_A(D_dense);
    Vector Ainv_bfi = (nA > 0) ? LinearSystemSolveFactorizedLLT(*factor_Ai1plusi2, bfi) : Vector::Zero(0);

    // ---------------------------------------------------------------------
    // 3) Build Fhat = F - E * Ainv_B
    // ---------------------------------------------------------------------
    DenseMatrix F_dense = sparse_to_dense(F, nB, nB);
    DenseMatrix E_dense = sparse_to_dense(E, nB, nA);
    DenseMatrix G_dense = sparse_to_dense(Gblk, nB, nC);
    DenseMatrix H_dense = sparse_to_dense(H, nB, nD);

    DenseMatrix Fhat = F_dense;
    if (nB > 0 && nA > 0) {
        Fhat.noalias() -= (E_dense * Ainv_B);
    }

    // ---------------------------------------------------------------------
    // 4) Build reduced RHS and reduced couplings for phi_mi
    // ---------------------------------------------------------------------
    Vector bmi = Vector::Zero(nB);
    if (vector_bi.size() == nB) {
        bmi = vector_bi;
    }

    Vector rhs_mi_hat = bmi;
    if (nB > 0 && nA > 0) {
        rhs_mi_hat.noalias() -= (E_dense * Ainv_bfi);
    }

    DenseMatrix Ghat = G_dense;
    DenseMatrix Hhat = H_dense;
    if (nB > 0 && nA > 0) {
        Ghat.noalias() -= (E_dense * Ainv_C);
        Hhat.noalias() -= (E_dense * Ainv_D);
    }

    // ---------------------------------------------------------------------
    // 5) Factorize Fhat and solve for z_mi, Y_mi_fb, Y_mi_mb
    // ---------------------------------------------------------------------
    DenseMatrix Y_mi_fb = DenseMatrix::Zero(nB, nC);
    DenseMatrix Y_mi_mb = DenseMatrix::Zero(nB, nD);
    Vector z_mi = Vector::Zero(nB);

    if (nB > 0) {
        SparseMatrix Fhat_sparse(nB, nB);
        std::vector<Triplet> trips_Fhat;
        trips_Fhat.reserve(static_cast<size_t>(nB) * static_cast<size_t>(nB));
        for (int i = 0; i < nB; ++i) {
            for (int j = 0; j < nB; ++j) {
                const double v = Fhat(i, j);
                if (v != 0.0) {
                    trips_Fhat.emplace_back(i, j, v);
                }
            }
        }
        Fhat_sparse.setFromTriplets(trips_Fhat.begin(), trips_Fhat.end());
        Fhat_sparse.makeCompressed();

        LinearSystemSparseLLTFactorization factor_Fhat;
        LinearSystemFactorizeLLT(Fhat_sparse, factor_Fhat);

        z_mi = LinearSystemSolveFactorizedLLT(factor_Fhat, rhs_mi_hat);
        Y_mi_fb = solve_A(DenseMatrix::Zero(0, 0)); // keep initialization style explicit
        Y_mi_mb = solve_A(DenseMatrix::Zero(0, 0));

        if (nC > 0) {
            Y_mi_fb = DenseMatrix::Zero(nB, nC);
            for (int c = 0; c < nC; ++c) {
                Vector rhs_col = Vector::Zero(nB);
                for (int r = 0; r < nB; ++r) {
                    rhs_col(r) = Ghat(r, c);
                }
                Vector x_col = LinearSystemSolveFactorizedLLT(factor_Fhat, rhs_col);
                for (int r = 0; r < nB; ++r) {
                    Y_mi_fb(r, c) = x_col(r);
                }
            }
        }
        if (nD > 0) {
            Y_mi_mb = DenseMatrix::Zero(nB, nD);
            for (int c = 0; c < nD; ++c) {
                Vector rhs_col = Vector::Zero(nB);
                for (int r = 0; r < nB; ++r) {
                    rhs_col(r) = Hhat(r, c);
                }
                Vector x_col = LinearSystemSolveFactorizedLLT(factor_Fhat, rhs_col);
                for (int r = 0; r < nB; ++r) {
                    Y_mi_mb(r, c) = x_col(r);
                }
            }
        }
    }

    // ---------------------------------------------------------------------
    // 6) Reconstruct phi_fi = z_fi - Y_fi_fb * phi_fb - Y_fi_mb * phi_mb
    // ---------------------------------------------------------------------
    Vector z_fi = Ainv_bfi;
    DenseMatrix Y_fi_fb = Ainv_C;
    DenseMatrix Y_fi_mb = Ainv_D;

    if (nA > 0 && nB > 0) {
        z_fi.noalias() -= (Ainv_B * z_mi);
    }
    if (nA > 0 && nC > 0) {
        Y_fi_fb.noalias() -= (Ainv_B * Y_mi_fb);
    }
    if (nA > 0 && nD > 0) {
        Y_fi_mb.noalias() -= (Ainv_B * Y_mi_mb);
    }

    // Full interior reconstruction data.
    vector_zi_DDP = Vector::Zero(nA + nB);
    if (nA > 0) {
        vector_zi_DDP.head(nA) = z_fi;
    }
    if (nB > 0) {
        vector_zi_DDP.tail(nB) = z_mi;
    }

    matrix_Y_DDP = DenseMatrix::Zero(nA + nB, nC + nD);
    if (nA > 0 && nC > 0) {
        matrix_Y_DDP.block(0, 0, nA, nC) = Y_fi_fb;
    }
    if (nA > 0 && nD > 0) {
        matrix_Y_DDP.block(0, nC, nA, nD) = Y_fi_mb;
    }
    if (nB > 0 && nC > 0) {
        matrix_Y_DDP.block(nA, 0, nB, nC) = Y_mi_fb;
    }
    if (nB > 0 && nD > 0) {
        matrix_Y_DDP.block(nA, nC, nB, nD) = Y_mi_mb;
    }

    // ---------------------------------------------------------------------
    // 7) Build the final condensed border system.
    //     [phi_fb ; phi_mb]
    // ---------------------------------------------------------------------
    DenseMatrix M_dense = sparse_to_dense(M, nC, nC);
    DenseMatrix N_dense = sparse_to_dense(N, nC, nD);
    DenseMatrix P_dense = sparse_to_dense(P, nD, nA);
    DenseMatrix Q_dense = sparse_to_dense(Q, nD, nB);
    DenseMatrix R_dense = sparse_to_dense(R, nD, nC);
    DenseMatrix T_dense = sparse_to_dense(T, nD, nD);

    DenseMatrix S_fb_fb = M_dense;
    DenseMatrix S_fb_mb = N_dense;
    DenseMatrix S_mb_fb = R_dense;
    DenseMatrix S_mb_mb = T_dense;

    if (nC > 0 && nA > 0) {
        S_fb_fb.noalias() -= (K * Y_fi_fb);
        S_fb_mb.noalias() -= (K * Y_fi_mb);
    }
    if (nC > 0 && nB > 0) {
        S_fb_fb.noalias() -= (L * Y_mi_fb);
        S_fb_mb.noalias() -= (L * Y_mi_mb);
    }
    if (nD > 0 && nA > 0) {
        S_mb_fb.noalias() -= (P_dense * Y_fi_fb);
        S_mb_mb.noalias() -= (P_dense * Y_fi_mb);
    }
    if (nD > 0 && nB > 0) {
        S_mb_fb.noalias() -= (Q_dense * Y_mi_fb);
        S_mb_mb.noalias() -= (Q_dense * Y_mi_mb);
    }

    Vector g_fb = Vector::Zero(nC);
    Vector g_mb = Vector::Zero(nD);
    if (nC > 0) {
        g_fb = Vector::Zero(nC);
        if (vector_b_fb.size() == nC) {
            g_fb = vector_b_fb;
        }
        if (nA > 0) {
            g_fb.noalias() -= (K * z_fi);
        }
        if (nB > 0) {
            g_fb.noalias() -= (L * z_mi);
        }
    }
    if (nD > 0) {
        g_mb = Vector::Zero(nD);
        if (vector_bb.size() == nD) {
            g_mb = vector_bb;
        }
        if (nA > 0) {
            g_mb.noalias() -= (P_dense * z_fi);
        }
        if (nB > 0) {
            g_mb.noalias() -= (Q_dense * z_mi);
        }
    }

    // ---------------------------------------------------------------------
    // 8) Store final condensed system and RHS.
    // ---------------------------------------------------------------------
    DenseMatrix dense_S = DenseMatrix::Zero(nC + nD, nC + nD);
    if (nC > 0) {
        dense_S.block(0, 0, nC, nC) = S_fb_fb;
        dense_S.block(0, nC, nC, nD) = S_fb_mb;
    }
    if (nD > 0) {
        dense_S.block(nC, 0, nD, nC) = S_mb_fb;
        dense_S.block(nC, nC, nD, nD) = S_mb_mb;
    }

    std::vector<Triplet> trips_S;
    trips_S.reserve(static_cast<size_t>(nC + nD) * static_cast<size_t>(nC + nD));
    for (int i = 0; i < nC + nD; ++i) {
        for (int j = 0; j < nC + nD; ++j) {
            const double v = dense_S(i, j);
            if (v != 0.0) {
                trips_S.emplace_back(i, j, v);
            }
        }
    }

    // matrix_S_DDP_precomputed.resize(nC + nD, nC + nD);
    // matrix_S_DDP_precomputed.setFromTriplets(trips_S.begin(), trips_S.end());
    // matrix_S_DDP_precomputed.makeCompressed();
    // matrix_S_DDP = matrix_S_DDP_precomputed;

    vector_zb_DDP = Vector::Zero(nC + nD);
    if (nC > 0) {
        vector_zb_DDP.head(nC) = g_fb;
    }
    if (nD > 0) {
        vector_zb_DDP.tail(nD) = g_mb;
    }
    vector_rhs_DDP = vector_zb_DDP;
}

void DDPModelCS::SolveCondensedDDPSystem() {
    clock_t t_total = clock();
    const int N_bor_global = Nfb + Nb;

    if (N_bor_global == 0) {
        vector_phi_b_DDP = Vector::Zero(0);
        return;
    }

    // Solve S * [phi_fb ; phi_mb] = rhs_DDP
    vector_phi_b_DDP = LinearSystemSolvingDense(matrix_S_DDP, vector_rhs_DDP);

    // Propagate matrix border part to the inherited EPM member
    if (Nb > 0) {
        vector_phi_b = vector_phi_b_DDP.tail(Nb);
    }

    // std::cout << "[PROFILE][DDPCS]   SolveCondensedDDPSystem: "
    //           << (double)(clock() - t_total) / CLOCKS_PER_SEC << " sec" << std::endl;
}

void DDPModelCS::ReconstructDDPInteriorFromBorder() {
    clock_t t_total = clock();
    const int N_int_global = Nfi + Ni;
    const int N_bor_global = Nfb + Nb;

    vector_phi_i_DDP = Vector::Zero(N_int_global);
    if (N_int_global == 0) {
        return;
    }

    // phi_i = z_i - Aii_g^{-1} * (Aib_g * phi_b)
    // Uses the sparse factorization precomputed in PrecomputeDDPSCForCurrentFrequency(),
    // avoiding the dense matrix_Y_DDP multiply (O(N_int*N_bor) -> O(nnz + N_int)).
    // std::cout << "[PROFILE][DDPCS]     Reconstruct sizes: N_int=" << N_int_global << " N_bor=" << N_bor_global
    //           << "  (Nfi=" << Nfi << " Ni=" << Ni << " | Nfb=" << Nfb << " Nb=" << Nb << ")" << std::endl;

    clock_t t_zero = clock();
    vector_phi_i_DDP = vector_zi_DDP;
    // std::cout << "[PROFILE][DDPCS]     Step 0 - Copy zi: "
    //           << (double)(clock() - t_zero) / CLOCKS_PER_SEC << " sec" << std::endl;

    clock_t t_reconstruct = clock();
    if (N_bor_global > 0 && vector_phi_b_DDP.size() == N_bor_global && factor_Ai1plusi2) {
        Vector tmp = Aib_g * vector_phi_b_DDP;
        Vector x   = LinearSystemSolveFactorizedLLT(*factor_Ai1plusi2, tmp);
        vector_phi_i_DDP.noalias() -= x;
    }
    // std::cout << "[PROFILE][DDPCS]     Step 1 - Sparse Aib*phi_b + solve Aii: "
    //           << (double)(clock() - t_reconstruct) / CLOCKS_PER_SEC << " sec" << std::endl;

    // Propagate each part to the inherited EPM member
    clock_t t_tail = clock();
    if (Ni > 0) {
        vector_phi_i = vector_phi_i_DDP.tail(Ni);
    }
    // std::cout << "[PROFILE][DDPCS]     Step 2 - Tail copy to phi_i: "
    //           << (double)(clock() - t_tail) / CLOCKS_PER_SEC << " sec" << std::endl;

//     std::cout << "[PROFILE][DDPCS]   Total ReconstructDDPInteriorFromBorder: "
//               << (double)(clock() - t_total) / CLOCKS_PER_SEC << " sec" << std::endl;
}

Vector DDPModelCS::BuildGlobalDDPSolutionVector() const {
    // Build a global DFN potential vector of size nb_nodes_fract
    if (nb_nodes_fract <= 0) { return Vector(); }
    Vector global = Vector::Zero(nb_nodes_fract);
    // fill interior fracture unknowns: phi_fi = vector_phi_i_DDP.head(Nfi)
    if ((unsigned)Nfi == interior_fract_global_indices.size() && vector_phi_i_DDP.size() == Nfi + Ni) {
        for (int k = 0; k < Nfi; ++k) {
            int gidx = interior_fract_global_indices[k];
            if (gidx >= 0 && gidx < global.size()) global(gidx) = vector_phi_i_DDP(k);
        }
    }
    // fill border fracture unknowns: phi_fb = vector_phi_b_DDP.head(Nfb)
    if ((unsigned)Nfb == border_fract_global_indices.size() && vector_phi_b_DDP.size() == Nfb + Nb) {
        for (int k = 0; k < Nfb; ++k) {
            int gidx = border_fract_global_indices[k];
            if (gidx >= 0 && gidx < global.size()) global(gidx) = vector_phi_b_DDP(k);
        }
    }
    return global;
}

DenseMatrix DDPModelCS::PotentialComputationDDPCS() {
    // Return the computed fracture potential matrix if available
    return potential_fracture;
}

bool DDPModelCS::IsInitializedDFNSC() const {
    return (Nfi > 0 || Nfb > 0) && (interior_fract_global_indices.size() + border_fract_global_indices.size() == (size_t)nb_nodes_fract);
}

void DDPModelCS::ExportDDPPotential(const std::string& output_dir,
                                    int elec_nb,
                                    double current_w)
{
    // Build the global matrix potential vector using the inherited EPM method.
    // vector_phi_i and vector_phi_b must have been populated by
    // ReconstructDDPInteriorFromBorder() and SolveCondensedDDPSystem() beforehand.
    Vector phi_global = BuildGlobalSolutionVectorFromSC();

    std::stringstream ss;
    ss << elec_nb;
    const std::string output_file = output_dir + std::string("ResultsDDP") + ss.str() + ".txt";
    std::ofstream output(output_file.c_str(), std::ofstream::app);
    if (!output.is_open()) {
        std::cerr << "[ERROR][DDPModelCS] Cannot open output file: " << output_file << std::endl;
        return;
    }

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            output << phi_global(return_index(i, (Ny - 1) - j, Ny)) << " ";
        }
        output << "\n";
    }
    output.close();
}

} // namespace eigen_cs
