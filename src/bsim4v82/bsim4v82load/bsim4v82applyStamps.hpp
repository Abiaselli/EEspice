#pragma once
#include "bsim4v82/bsim4v82.hpp"
#include <armadillo>

namespace bsim4 {
// For safe modification of MNA matrix
/*  
    These checks are critical to prevent out-of-bounds access, 
    which could lead to runtime errors or undefined behavior, especially
    since negative indices are invalid for Armadillo matrices.
    Each stamp follows this template:
    Stamp(LHS, {rowNode}, {colNode},
        BSIM4{Row}{Col}, nodeValid, 
        BSIM4V82::{ROW_NODE_TYPE}, BSIM4V82::{COL_NODE_TYPE});
*/
inline void Stamp(arma::mat &mat, int row, int col, double val,
                    const std::array<bool, 12> &BSIM4nodeValid, BSIM4V82::NodeType row_node, BSIM4V82::NodeType col_node)
{
    if (BSIM4nodeValid[row_node] && BSIM4nodeValid[col_node]){
        mat(row, col) += val;
    }
}
inline void StampRHS(arma::vec &vec, int index, double val,
                     const std::array<bool, 12> &BSIM4nodeValid, BSIM4V82::NodeType node)
{
    if (BSIM4nodeValid[node]) {
        vec[index] += val;
    }
}
inline void StampOffset(arma::mat &mat, int row, int col, double val,
                    const std::array<bool, 12> &BSIM4nodeValid, BSIM4V82::NodeType row_node, BSIM4V82::NodeType col_node)
{
    if (BSIM4nodeValid[row_node] && BSIM4nodeValid[col_node]){
        double* p = mat.memptr();
        p[row + col * mat.n_rows] += val; // column-major storage
    }
}

/**
 * @brief Applies the calculated conductances and currents to the system matrices.
 *
 * This function takes a BSIM4stamp struct containing pre-calculated values and "stamps" them
 * into the appropriate locations in the Left-Hand-Side (LHS) and Right-Hand-Side (RHS) matrices
 * of the Modified Nodal Analysis (MNA) system.
 * @param instance The BSIM4 device instance, used to get node indices.
 * @param stamps A const reference to the struct holding all the values to be stamped.
 * @param LHS The MNA conductance matrix (to be modified).
 * @param RHS The MNA current/voltage source vector (to be modified).
 */
void bsim4applyStamps(const BSIM4V82 &instance, const BSIM4stamp &s, 
                      arma::mat &LHS, arma::vec &RHS)
{
    if (!s.load) return; // do not load stamps into the matrix

    // Node indexes in mna matrix
    auto get_NodeIndex = [](int node) -> int {
        return   node - 1;
    };
    const int dNode = get_NodeIndex(instance.BSIM4dNode);
    const int gNodeExt = get_NodeIndex(instance.BSIM4gNodeExt);
    const int sNode = get_NodeIndex(instance.BSIM4sNode);
    const int bNode = get_NodeIndex(instance.BSIM4bNode);
    const int dNodePrime = get_NodeIndex(instance.BSIM4dNodePrime);
    const int gNodePrime = get_NodeIndex(instance.BSIM4gNodePrime);
    const int gNodeMid = get_NodeIndex(instance.BSIM4gNodeMid);
    const int sNodePrime = get_NodeIndex(instance.BSIM4sNodePrime);
    const int bNodePrime = get_NodeIndex(instance.BSIM4bNodePrime);
    const int dbNode = get_NodeIndex(instance.BSIM4dbNode);
    const int sbNode = get_NodeIndex(instance.BSIM4sbNode);
    const int qNode = get_NodeIndex(instance.BSIM4qNode);

    // Now stamping to the MNA matrixes
    const std::array<bool,12> &nodeValid = instance.BSIM4nodeValid;
    // RHS Stamps
    // RHSdNodePrime 
    StampRHS(RHS, dNodePrime, s.RHSdNodePrime, nodeValid, BSIM4V82::D_NODE_PRIME);
    
    // RHSgNodePrime 
    StampRHS(RHS, gNodePrime, s.RHSgNodePrime, nodeValid, BSIM4V82::G_NODE_PRIME);
    
    // RHSgNodeExt 
    StampRHS(RHS, gNodeExt, s.RHSgNodeExt, nodeValid, BSIM4V82::G_NODE_EXT);
    
    // RHSgNodeMid 
    StampRHS(RHS, gNodeMid, s.RHSgNodeMid, nodeValid, BSIM4V82::G_NODE_MID);
    
    // RHSbNodePrime 
    StampRHS(RHS, bNodePrime, s.RHSbNodePrime, nodeValid, BSIM4V82::B_NODE_PRIME);
    
    // RHSsNodePrime 
    StampRHS(RHS, sNodePrime, s.RHSsNodePrime, nodeValid, BSIM4V82::S_NODE_PRIME);
    
    // RHSdbNode 
    StampRHS(RHS, dbNode, s.RHSdbNode, nodeValid, BSIM4V82::DB_NODE);
    
    // RHSsbNode
    StampRHS(RHS, sbNode, s.RHSsbNode, nodeValid, BSIM4V82::SB_NODE);
        
    // RHSdNode 
    StampRHS(RHS, dNode, s.RHSdNode, nodeValid, BSIM4V82::D_NODE);
    
    // RHSsNode 
    StampRHS(RHS, sNode, s.RHSsNode, nodeValid, BSIM4V82::S_NODE);
    
    // RHSqNode
    StampRHS(RHS, qNode, s.RHSqNode, nodeValid, BSIM4V82::Q_NODE);

    // LHS Stamps
    // BSIM4DPdPtr
    Stamp(LHS, dNodePrime, dNode,
        s.BSIM4DPd, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::D_NODE);
    // BSIM4DPdpPtr
    Stamp(LHS, dNodePrime, dNodePrime,
        s.BSIM4DPdp, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::D_NODE_PRIME);
    
    // BSIM4DPgpPtr
    Stamp(LHS, dNodePrime, gNodePrime,
        s.BSIM4DPgp, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::G_NODE_PRIME);
    // BSIM4DPgmPtr
    Stamp(LHS, dNodePrime, gNodeMid,
        s.BSIM4DPgm, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::G_NODE_MID);
    // BSIM4DPspPtr
    Stamp(LHS, dNodePrime, sNodePrime,
        s.BSIM4DPsp, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::S_NODE_PRIME);
    // BSIM4DPbpPtr
    Stamp(LHS, dNodePrime, bNodePrime,
        s.BSIM4DPbp, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::B_NODE_PRIME);
    // BSIM4DPdbPtr
    Stamp(LHS, dNodePrime, dbNode,
        s.BSIM4DPdb, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::DB_NODE);

    // BSIM4DdPtr
    Stamp(LHS, dNode, dNode,
        s.BSIM4Dd, nodeValid, 
        BSIM4V82::D_NODE, BSIM4V82::D_NODE);
    // BSIM4DdpPtr
    Stamp(LHS, dNode, dNodePrime,
        s.BSIM4Ddp, nodeValid, 
        BSIM4V82::D_NODE, BSIM4V82::D_NODE_PRIME);

    // BSIM4GPdpPtr
    Stamp(LHS, gNodePrime, dNodePrime,
        s.BSIM4GPdp, nodeValid, 
        BSIM4V82::G_NODE_PRIME, BSIM4V82::D_NODE_PRIME);
    // BSIM4GPgpPtr
    Stamp(LHS, gNodePrime, gNodePrime,
        s.BSIM4GPgp, nodeValid, 
        BSIM4V82::G_NODE_PRIME, BSIM4V82::G_NODE_PRIME);
    // BSIM4GPgmPtr
    Stamp(LHS, gNodePrime, gNodeMid,
        s.BSIM4GPgm, nodeValid, 
        BSIM4V82::G_NODE_PRIME, BSIM4V82::G_NODE_MID);
    // BSIM4GPgePtr
    Stamp(LHS, gNodePrime, gNodeExt,
        s.BSIM4GPge, nodeValid, 
        BSIM4V82::G_NODE_PRIME, BSIM4V82::G_NODE_EXT);
    // BSIM4GPspPtr
    Stamp(LHS, gNodePrime, sNodePrime,
        s.BSIM4GPsp, nodeValid, 
        BSIM4V82::G_NODE_PRIME, BSIM4V82::S_NODE_PRIME);
    // BSIM4GPbpPtr
    Stamp(LHS, gNodePrime, bNodePrime,
        s.BSIM4GPbp, nodeValid, 
        BSIM4V82::G_NODE_PRIME, BSIM4V82::B_NODE_PRIME);

    // BSIM4GMdpPtr
    Stamp(LHS, gNodeMid, dNodePrime,
        s.BSIM4GMdp, nodeValid, 
        BSIM4V82::G_NODE_MID, BSIM4V82::D_NODE_PRIME);
    // BSIM4GMgpPtr
    Stamp(LHS, gNodeMid, gNodePrime,
        s.BSIM4GMgp, nodeValid, 
        BSIM4V82::G_NODE_MID, BSIM4V82::G_NODE_PRIME);
    // BSIM4GMgmPtr
    Stamp(LHS, gNodeMid, gNodeMid,
        s.BSIM4GMgm, nodeValid, 
        BSIM4V82::G_NODE_MID, BSIM4V82::G_NODE_MID);
    // BSIM4GMgePtr
    Stamp(LHS, gNodeMid, gNodeExt,
        s.BSIM4GMge, nodeValid, 
        BSIM4V82::G_NODE_MID, BSIM4V82::G_NODE_EXT);
    // BSIM4GMspPtr
    Stamp(LHS, gNodeMid, sNodePrime,
        s.BSIM4GMsp, nodeValid, 
        BSIM4V82::G_NODE_MID, BSIM4V82::S_NODE_PRIME);
    // BSIM4GMbpPtr
    Stamp(LHS, gNodeMid, bNodePrime,
        s.BSIM4GMbp, nodeValid, 
        BSIM4V82::G_NODE_MID, BSIM4V82::B_NODE_PRIME);

    // BSIM4GEdpPtr
    Stamp(LHS, gNodeExt, dNodePrime,
        s.BSIM4GEdp, nodeValid, 
        BSIM4V82::G_NODE_EXT, BSIM4V82::D_NODE_PRIME);
    // BSIM4GEgpPtr
    Stamp(LHS, gNodeExt, gNodePrime,
        s.BSIM4GEgp, nodeValid, 
        BSIM4V82::G_NODE_EXT, BSIM4V82::G_NODE_PRIME);
    // BSIM4GEgmPtr
    Stamp(LHS, gNodeExt, gNodeMid,
        s.BSIM4GEgm, nodeValid, 
        BSIM4V82::G_NODE_EXT, BSIM4V82::G_NODE_MID);
    // BSIM4GEgePtr
    Stamp(LHS, gNodeExt, gNodeExt,
        s.BSIM4GEge, nodeValid, 
        BSIM4V82::G_NODE_EXT, BSIM4V82::G_NODE_EXT);
    // BSIM4GEspPtr
    Stamp(LHS, gNodeExt, sNodePrime,
        s.BSIM4GEsp, nodeValid, 
        BSIM4V82::G_NODE_EXT, BSIM4V82::S_NODE_PRIME);
    // BSIM4GEbpPtr
    Stamp(LHS, gNodeExt, bNodePrime,
        s.BSIM4GEbp, nodeValid, 
        BSIM4V82::G_NODE_EXT, BSIM4V82::B_NODE_PRIME);

    // BSIM4SPdpPtr
    Stamp(LHS, sNodePrime, dNodePrime,
        s.BSIM4SPdp, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::D_NODE_PRIME);
    // BSIM4SPgpPtr
    Stamp(LHS, sNodePrime, gNodePrime,
        s.BSIM4SPgp, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::G_NODE_PRIME);
    // BSIM4SPgmPtr
    Stamp(LHS, sNodePrime, gNodeMid,
        s.BSIM4SPgm, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::G_NODE_MID);
    // BSIM4SPsPtr
    Stamp(LHS, sNodePrime, sNode,
        s.BSIM4SPs, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::S_NODE);
    // BSIM4SPspPtr
    Stamp(LHS, sNodePrime, sNodePrime,
        s.BSIM4SPsp, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::S_NODE_PRIME);
    // BSIM4SPbpPtr
    Stamp(LHS, sNodePrime, bNodePrime,
        s.BSIM4SPbp, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::B_NODE_PRIME);
    // BSIM4SPsbPtr
    Stamp(LHS, sNodePrime, sbNode,
        s.BSIM4SPsb, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::SB_NODE);

    // BSIM4SspPtr
    Stamp(LHS, sNode, sNodePrime,
        s.BSIM4Ssp, nodeValid, 
        BSIM4V82::S_NODE, BSIM4V82::S_NODE_PRIME);
    // BSIM4SsPtr
    Stamp(LHS, sNode, sNode,
        s.BSIM4Ss, nodeValid, 
        BSIM4V82::S_NODE, BSIM4V82::S_NODE);

    // BSIM4BPdpPtr
    Stamp(LHS, bNodePrime, dNodePrime,
        s.BSIM4BPdp, nodeValid, 
        BSIM4V82::B_NODE_PRIME, BSIM4V82::D_NODE_PRIME);
    // BSIM4BPgpPtr
    Stamp(LHS, bNodePrime, gNodePrime,
        s.BSIM4BPgp, nodeValid, 
        BSIM4V82::B_NODE_PRIME, BSIM4V82::G_NODE_PRIME);
    // BSIM4BPgmPtr
    Stamp(LHS, bNodePrime, gNodeMid,
        s.BSIM4BPgm, nodeValid, 
        BSIM4V82::B_NODE_PRIME, BSIM4V82::G_NODE_MID);
    // BSIM4BPspPtr
    Stamp(LHS, bNodePrime, sNodePrime,
        s.BSIM4BPsp, nodeValid, 
        BSIM4V82::B_NODE_PRIME, BSIM4V82::S_NODE_PRIME);
    // BSIM4BPdbPtr
    Stamp(LHS, bNodePrime, dbNode,
        s.BSIM4BPdb, nodeValid, 
        BSIM4V82::B_NODE_PRIME, BSIM4V82::DB_NODE);
    // BSIM4BPbPtr
    Stamp(LHS, bNodePrime, bNode,
        s.BSIM4BPb, nodeValid, 
        BSIM4V82::B_NODE_PRIME, BSIM4V82::B_NODE);
    // BSIM4BPsbPtr
    Stamp(LHS, bNodePrime, sbNode,
        s.BSIM4BPsb, nodeValid, 
        BSIM4V82::B_NODE_PRIME, BSIM4V82::SB_NODE);
    // BSIM4BPbpPtr
    Stamp(LHS, bNodePrime, bNodePrime,
        s.BSIM4BPbp, nodeValid, 
        BSIM4V82::B_NODE_PRIME, BSIM4V82::B_NODE_PRIME);

    // BSIM4DBdpPtr
    Stamp(LHS, dbNode, dNodePrime,
        s.BSIM4DBdp, nodeValid, 
        BSIM4V82::DB_NODE, BSIM4V82::D_NODE_PRIME);
    // BSIM4DBdbPtr
    Stamp(LHS, dbNode, dbNode,
        s.BSIM4DBdb, nodeValid, 
        BSIM4V82::DB_NODE, BSIM4V82::DB_NODE);
    // BSIM4DBbpPtr
    Stamp(LHS, dbNode, bNodePrime,
        s.BSIM4DBbp, nodeValid, 
        BSIM4V82::DB_NODE, BSIM4V82::B_NODE_PRIME);
    // BSIM4DBbPtr
    Stamp(LHS, dbNode, bNode,
        s.BSIM4DBb, nodeValid, 
        BSIM4V82::DB_NODE, BSIM4V82::B_NODE);

    // BSIM4SBspPtr
    Stamp(LHS, sbNode, sNodePrime,
        s.BSIM4SBsp, nodeValid, 
        BSIM4V82::SB_NODE, BSIM4V82::S_NODE_PRIME);
    // BSIM4SBbpPtr
    Stamp(LHS, sbNode, bNodePrime,
        s.BSIM4SBbp, nodeValid, 
        BSIM4V82::SB_NODE, BSIM4V82::B_NODE_PRIME);
    // BSIM4SBbPtr
    Stamp(LHS, sbNode, bNode,
        s.BSIM4SBb, nodeValid, 
        BSIM4V82::SB_NODE, BSIM4V82::B_NODE);
    // BSIM4SBsbPtr
    Stamp(LHS, sbNode, sbNode,
        s.BSIM4SBsb, nodeValid, 
        BSIM4V82::SB_NODE, BSIM4V82::SB_NODE);

    // BSIM4BdbPtr
    Stamp(LHS, bNode, dbNode,
        s.BSIM4Bdb, nodeValid, 
        BSIM4V82::B_NODE, BSIM4V82::DB_NODE);
    // BSIM4BbpPtr
    Stamp(LHS, bNode, bNodePrime,
        s.BSIM4Bbp, nodeValid, 
        BSIM4V82::B_NODE, BSIM4V82::B_NODE_PRIME);
    // BSIM4BsbPtr
    Stamp(LHS, bNode, sbNode,
        s.BSIM4Bsb, nodeValid, 
        BSIM4V82::B_NODE, BSIM4V82::SB_NODE);
    // BSIM4BbPtr
    Stamp(LHS, bNode, bNode,
        s.BSIM4Bb, nodeValid, 
        BSIM4V82::B_NODE, BSIM4V82::B_NODE);

    // BSIM4DgpPtr
    Stamp(LHS, dNode, gNodePrime,
        s.BSIM4Dgp, nodeValid, 
        BSIM4V82::D_NODE, BSIM4V82::G_NODE_PRIME);
    // BSIM4DspPtr
    Stamp(LHS, dNode, sNodePrime,
        s.BSIM4Dsp, nodeValid, 
        BSIM4V82::D_NODE, BSIM4V82::S_NODE_PRIME);
    // BSIM4DbpPtr
    Stamp(LHS, dNode, bNodePrime,
        s.BSIM4Dbp, nodeValid, 
        BSIM4V82::D_NODE, BSIM4V82::B_NODE_PRIME);
    // BSIM4SdpPtr
    Stamp(LHS, sNode, dNodePrime,
        s.BSIM4Sdp, nodeValid, 
        BSIM4V82::S_NODE, BSIM4V82::D_NODE_PRIME);
    // BSIM4SgpPtr
    Stamp(LHS, sNode, gNodePrime,
        s.BSIM4Sgp, nodeValid, 
        BSIM4V82::S_NODE, BSIM4V82::G_NODE_PRIME);
    // BSIM4SbpPtr
    Stamp(LHS, sNode, bNodePrime,
        s.BSIM4Sbp, nodeValid, 
        BSIM4V82::S_NODE, BSIM4V82::B_NODE_PRIME);

    // BSIM4QdpPtr
    Stamp(LHS, qNode, dNodePrime,
        s.BSIM4Qdp, nodeValid, 
        BSIM4V82::Q_NODE, BSIM4V82::D_NODE_PRIME);
    // BSIM4QgpPtr
    Stamp(LHS, qNode, gNodePrime,
        s.BSIM4Qgp, nodeValid, 
        BSIM4V82::Q_NODE, BSIM4V82::G_NODE_PRIME);
    // BSIM4QspPtr
    Stamp(LHS, qNode, sNodePrime,
        s.BSIM4Qsp, nodeValid, 
        BSIM4V82::Q_NODE, BSIM4V82::S_NODE_PRIME);
    // BSIM4QbpPtr
    Stamp(LHS, qNode, bNodePrime,
        s.BSIM4Qbp, nodeValid, 
        BSIM4V82::Q_NODE, BSIM4V82::B_NODE_PRIME);
    // BSIM4QqPtr
    Stamp(LHS, qNode, qNode,
        s.BSIM4Qq, nodeValid, 
        BSIM4V82::Q_NODE, BSIM4V82::Q_NODE);
    // BSIM4DPqPtr
    Stamp(LHS, dNodePrime, qNode,
        s.BSIM4DPq, nodeValid, 
        BSIM4V82::D_NODE_PRIME, BSIM4V82::Q_NODE);
    // BSIM4GPqPtr
    Stamp(LHS, gNodePrime, qNode,
        s.BSIM4GPq, nodeValid, 
        BSIM4V82::G_NODE_PRIME, BSIM4V82::Q_NODE);
    // BSIM4SPqPtr
    Stamp(LHS, sNodePrime, qNode,
        s.BSIM4SPq, nodeValid, 
        BSIM4V82::S_NODE_PRIME, BSIM4V82::Q_NODE);

}

} // namespace bsim4