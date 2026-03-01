#pragma once
#include "bsim4v82/bsim4v82.hpp"
#include <armadillo>
#include "hybrid_matrix.hpp"

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

// Record position for pattern discovery (no value stamped)
inline void RecordPosition(HybridMatrix &mat, int row, int col,
                    const std::array<bool, 12> &BSIM4nodeValid,
                    BSIM4V82::NodeType row_node, BSIM4V82::NodeType col_node)
{
    if (BSIM4nodeValid[row_node] && BSIM4nodeValid[col_node]){
        mat.record_position(row, col);
    }
}

inline void Stamp(HybridMatrix &mat, int row, int col, double val,
                    const std::array<bool, 12> &BSIM4nodeValid, BSIM4V82::NodeType row_node, BSIM4V82::NodeType col_node)
{
    if (BSIM4nodeValid[row_node] && BSIM4nodeValid[col_node]){
        mat.add_stamp_indexed(row, col, val);
    }
}
inline void StampRHS(arma::vec &vec, int index, double val,
                     const std::array<bool, 12> &BSIM4nodeValid, BSIM4V82::NodeType node)
{
    if (BSIM4nodeValid[node]) {
        vec[index] += val;
    }
}
inline void StampOffset(HybridMatrix &mat, int row, int col, double val,
                    const std::array<bool, 12> &BSIM4nodeValid, BSIM4V82::NodeType row_node, BSIM4V82::NodeType col_node)
{
    if (BSIM4nodeValid[row_node] && BSIM4nodeValid[col_node]){
        // Use HybridMatrix's stamping method for portability
        mat.add_stamp(row, col, val);
    }
}

// Single source of truth for BSIM4 LHS stamps (field, rowIndexVar, colIndexVar, rowNodeType, colNodeType)
#define BSIM4_LHS_STAMP_LIST(X) \
    X(BSIM4DPd, dNodePrime, dNode, BSIM4V82::D_NODE_PRIME, BSIM4V82::D_NODE) \
    X(BSIM4DPdp, dNodePrime, dNodePrime, BSIM4V82::D_NODE_PRIME, BSIM4V82::D_NODE_PRIME) \
    X(BSIM4DPgp, dNodePrime, gNodePrime, BSIM4V82::D_NODE_PRIME, BSIM4V82::G_NODE_PRIME) \
    X(BSIM4DPgm, dNodePrime, gNodeMid, BSIM4V82::D_NODE_PRIME, BSIM4V82::G_NODE_MID) \
    X(BSIM4DPsp, dNodePrime, sNodePrime, BSIM4V82::D_NODE_PRIME, BSIM4V82::S_NODE_PRIME) \
    X(BSIM4DPbp, dNodePrime, bNodePrime, BSIM4V82::D_NODE_PRIME, BSIM4V82::B_NODE_PRIME) \
    X(BSIM4DPdb, dNodePrime, dbNode, BSIM4V82::D_NODE_PRIME, BSIM4V82::DB_NODE) \
    X(BSIM4Dd, dNode, dNode, BSIM4V82::D_NODE, BSIM4V82::D_NODE) \
    X(BSIM4Ddp, dNode, dNodePrime, BSIM4V82::D_NODE, BSIM4V82::D_NODE_PRIME) \
    X(BSIM4GPdp, gNodePrime, dNodePrime, BSIM4V82::G_NODE_PRIME, BSIM4V82::D_NODE_PRIME) \
    X(BSIM4GPgp, gNodePrime, gNodePrime, BSIM4V82::G_NODE_PRIME, BSIM4V82::G_NODE_PRIME) \
    X(BSIM4GPgm, gNodePrime, gNodeMid, BSIM4V82::G_NODE_PRIME, BSIM4V82::G_NODE_MID) \
    X(BSIM4GPge, gNodePrime, gNodeExt, BSIM4V82::G_NODE_PRIME, BSIM4V82::G_NODE_EXT) \
    X(BSIM4GPsp, gNodePrime, sNodePrime, BSIM4V82::G_NODE_PRIME, BSIM4V82::S_NODE_PRIME) \
    X(BSIM4GPbp, gNodePrime, bNodePrime, BSIM4V82::G_NODE_PRIME, BSIM4V82::B_NODE_PRIME) \
    X(BSIM4GMdp, gNodeMid, dNodePrime, BSIM4V82::G_NODE_MID, BSIM4V82::D_NODE_PRIME) \
    X(BSIM4GMgp, gNodeMid, gNodePrime, BSIM4V82::G_NODE_MID, BSIM4V82::G_NODE_PRIME) \
    X(BSIM4GMgm, gNodeMid, gNodeMid, BSIM4V82::G_NODE_MID, BSIM4V82::G_NODE_MID) \
    X(BSIM4GMge, gNodeMid, gNodeExt, BSIM4V82::G_NODE_MID, BSIM4V82::G_NODE_EXT) \
    X(BSIM4GMsp, gNodeMid, sNodePrime, BSIM4V82::G_NODE_MID, BSIM4V82::S_NODE_PRIME) \
    X(BSIM4GMbp, gNodeMid, bNodePrime, BSIM4V82::G_NODE_MID, BSIM4V82::B_NODE_PRIME) \
    X(BSIM4GEdp, gNodeExt, dNodePrime, BSIM4V82::G_NODE_EXT, BSIM4V82::D_NODE_PRIME) \
    X(BSIM4GEgp, gNodeExt, gNodePrime, BSIM4V82::G_NODE_EXT, BSIM4V82::G_NODE_PRIME) \
    X(BSIM4GEgm, gNodeExt, gNodeMid, BSIM4V82::G_NODE_EXT, BSIM4V82::G_NODE_MID) \
    X(BSIM4GEge, gNodeExt, gNodeExt, BSIM4V82::G_NODE_EXT, BSIM4V82::G_NODE_EXT) \
    X(BSIM4GEsp, gNodeExt, sNodePrime, BSIM4V82::G_NODE_EXT, BSIM4V82::S_NODE_PRIME) \
    X(BSIM4GEbp, gNodeExt, bNodePrime, BSIM4V82::G_NODE_EXT, BSIM4V82::B_NODE_PRIME) \
    X(BSIM4SPdp, sNodePrime, dNodePrime, BSIM4V82::S_NODE_PRIME, BSIM4V82::D_NODE_PRIME) \
    X(BSIM4SPgp, sNodePrime, gNodePrime, BSIM4V82::S_NODE_PRIME, BSIM4V82::G_NODE_PRIME) \
    X(BSIM4SPgm, sNodePrime, gNodeMid, BSIM4V82::S_NODE_PRIME, BSIM4V82::G_NODE_MID) \
    X(BSIM4SPs, sNodePrime, sNode, BSIM4V82::S_NODE_PRIME, BSIM4V82::S_NODE) \
    X(BSIM4SPsp, sNodePrime, sNodePrime, BSIM4V82::S_NODE_PRIME, BSIM4V82::S_NODE_PRIME) \
    X(BSIM4SPbp, sNodePrime, bNodePrime, BSIM4V82::S_NODE_PRIME, BSIM4V82::B_NODE_PRIME) \
    X(BSIM4SPsb, sNodePrime, sbNode, BSIM4V82::S_NODE_PRIME, BSIM4V82::SB_NODE) \
    X(BSIM4Ssp, sNode, sNodePrime, BSIM4V82::S_NODE, BSIM4V82::S_NODE_PRIME) \
    X(BSIM4Ss, sNode, sNode, BSIM4V82::S_NODE, BSIM4V82::S_NODE) \
    X(BSIM4BPdp, bNodePrime, dNodePrime, BSIM4V82::B_NODE_PRIME, BSIM4V82::D_NODE_PRIME) \
    X(BSIM4BPgp, bNodePrime, gNodePrime, BSIM4V82::B_NODE_PRIME, BSIM4V82::G_NODE_PRIME) \
    X(BSIM4BPgm, bNodePrime, gNodeMid, BSIM4V82::B_NODE_PRIME, BSIM4V82::G_NODE_MID) \
    X(BSIM4BPsp, bNodePrime, sNodePrime, BSIM4V82::B_NODE_PRIME, BSIM4V82::S_NODE_PRIME) \
    X(BSIM4BPdb, bNodePrime, dbNode, BSIM4V82::B_NODE_PRIME, BSIM4V82::DB_NODE) \
    X(BSIM4BPb, bNodePrime, bNode, BSIM4V82::B_NODE_PRIME, BSIM4V82::B_NODE) \
    X(BSIM4BPsb, bNodePrime, sbNode, BSIM4V82::B_NODE_PRIME, BSIM4V82::SB_NODE) \
    X(BSIM4BPbp, bNodePrime, bNodePrime, BSIM4V82::B_NODE_PRIME, BSIM4V82::B_NODE_PRIME) \
    X(BSIM4DBdp, dbNode, dNodePrime, BSIM4V82::DB_NODE, BSIM4V82::D_NODE_PRIME) \
    X(BSIM4DBdb, dbNode, dbNode, BSIM4V82::DB_NODE, BSIM4V82::DB_NODE) \
    X(BSIM4DBbp, dbNode, bNodePrime, BSIM4V82::DB_NODE, BSIM4V82::B_NODE_PRIME) \
    X(BSIM4DBb, dbNode, bNode, BSIM4V82::DB_NODE, BSIM4V82::B_NODE) \
    X(BSIM4SBsp, sbNode, sNodePrime, BSIM4V82::SB_NODE, BSIM4V82::S_NODE_PRIME) \
    X(BSIM4SBbp, sbNode, bNodePrime, BSIM4V82::SB_NODE, BSIM4V82::B_NODE_PRIME) \
    X(BSIM4SBb, sbNode, bNode, BSIM4V82::SB_NODE, BSIM4V82::B_NODE) \
    X(BSIM4SBsb, sbNode, sbNode, BSIM4V82::SB_NODE, BSIM4V82::SB_NODE) \
    X(BSIM4Bdb, bNode, dbNode, BSIM4V82::B_NODE, BSIM4V82::DB_NODE) \
    X(BSIM4Bbp, bNode, bNodePrime, BSIM4V82::B_NODE, BSIM4V82::B_NODE_PRIME) \
    X(BSIM4Bsb, bNode, sbNode, BSIM4V82::B_NODE, BSIM4V82::SB_NODE) \
    X(BSIM4Bb, bNode, bNode, BSIM4V82::B_NODE, BSIM4V82::B_NODE) \
    X(BSIM4Dgp, dNode, gNodePrime, BSIM4V82::D_NODE, BSIM4V82::G_NODE_PRIME) \
    X(BSIM4Dsp, dNode, sNodePrime, BSIM4V82::D_NODE, BSIM4V82::S_NODE_PRIME) \
    X(BSIM4Dbp, dNode, bNodePrime, BSIM4V82::D_NODE, BSIM4V82::B_NODE_PRIME) \
    X(BSIM4Sdp, sNode, dNodePrime, BSIM4V82::S_NODE, BSIM4V82::D_NODE_PRIME) \
    X(BSIM4Sgp, sNode, gNodePrime, BSIM4V82::S_NODE, BSIM4V82::G_NODE_PRIME) \
    X(BSIM4Sbp, sNode, bNodePrime, BSIM4V82::S_NODE, BSIM4V82::B_NODE_PRIME) \
    X(BSIM4Qdp, qNode, dNodePrime, BSIM4V82::Q_NODE, BSIM4V82::D_NODE_PRIME) \
    X(BSIM4Qgp, qNode, gNodePrime, BSIM4V82::Q_NODE, BSIM4V82::G_NODE_PRIME) \
    X(BSIM4Qsp, qNode, sNodePrime, BSIM4V82::Q_NODE, BSIM4V82::S_NODE_PRIME) \
    X(BSIM4Qbp, qNode, bNodePrime, BSIM4V82::Q_NODE, BSIM4V82::B_NODE_PRIME) \
    X(BSIM4Qq, qNode, qNode, BSIM4V82::Q_NODE, BSIM4V82::Q_NODE) \
    X(BSIM4DPq, dNodePrime, qNode, BSIM4V82::D_NODE_PRIME, BSIM4V82::Q_NODE) \
    X(BSIM4GPq, gNodePrime, qNode, BSIM4V82::G_NODE_PRIME, BSIM4V82::Q_NODE) \
    X(BSIM4SPq, sNodePrime, qNode, BSIM4V82::S_NODE_PRIME, BSIM4V82::Q_NODE)

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
                      HybridMatrix &LHS, arma::vec &RHS)
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
#define BSIM4_APPLY_LHS_STAMP(field, rowNode, colNode, rowType, colType) \
    Stamp(LHS, rowNode, colNode, s.field, nodeValid, rowType, colType);
    BSIM4_LHS_STAMP_LIST(BSIM4_APPLY_LHS_STAMP)
#undef BSIM4_APPLY_LHS_STAMP

}

/**
 * @brief Records all stamp positions for a BSIM4 instance (pattern discovery phase)
 *
 * This function records all matrix positions that will be stamped by bsim4applyStamps,
 * without actually writing any values. Used to build the pre-allocated sparse structure.
 */
void bsim4RecordPattern(const BSIM4V82 &instance, HybridMatrix &LHS)
{
    auto get_NodeIndex = [](int node) -> int {
        return node - 1;
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

    const std::array<bool,12> &nodeValid = instance.BSIM4nodeValid;

    // Record all LHS stamp positions (same pattern as bsim4applyStamps)
#define BSIM4_RECORD_LHS_STAMP(field, rowNode, colNode, rowType, colType) \
    RecordPosition(LHS, rowNode, colNode, nodeValid, rowType, colType);
    BSIM4_LHS_STAMP_LIST(BSIM4_RECORD_LHS_STAMP)
#undef BSIM4_RECORD_LHS_STAMP
}

/**
 * @brief Builds per-instance cache of sparse CSC values[] indices for BSIM4 LHS stamps.
 *
 * Must be called after pattern discovery and LHS.lock_pattern().
 * If any required stamped position is missing from the locked pattern, cache.built=false.
 */
void bsim4BuildStampIndexCache(const BSIM4V82 &instance, const HybridMatrix &LHS, BSIM4StampIndexCache &cache)
{
    cache = BSIM4StampIndexCache{};
    if (!LHS.is_sparse() || !LHS.is_pattern_locked()) {
        return;
    }

    auto get_NodeIndex = [](int node) -> int {
        return node - 1;
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

    const std::array<bool, 12> &nodeValid = instance.BSIM4nodeValid;

    bool ok = true;

#define BSIM4_CACHE_LHS_STAMP(field, rowNode, colNode, rowType, colType)                                     \
    do {                                                                                                    \
        if (nodeValid[rowType] && nodeValid[colType]) {                                                     \
            size_t idx = 0;                                                                                 \
            if (LHS.try_get_position_index(static_cast<size_t>(rowNode), static_cast<size_t>(colNode), idx)) { \
                cache.field = idx;                                                                          \
            } else {                                                                                        \
                cache.field = BSIM4StampIndexCache::kInvalid;                                               \
                ok = false;                                                                                 \
            }                                                                                               \
        } else {                                                                                            \
            cache.field = BSIM4StampIndexCache::kInvalid;                                                   \
        }                                                                                                   \
    } while (0);
    BSIM4_LHS_STAMP_LIST(BSIM4_CACHE_LHS_STAMP)
#undef BSIM4_CACHE_LHS_STAMP

    cache.built = ok;
}

/**
 * @brief Fast stamping path using pre-cached CSC values[] indices (sparse only).
 *
 * Falls back to bsim4applyStamps if cache is not usable.
 */
void bsim4applyStampsCached(const BSIM4V82 &instance, const BSIM4stamp &s,
                      const BSIM4StampIndexCache &cache,
                      HybridMatrix &LHS, arma::vec &RHS)
{
    if (!s.load) return; // do not load stamps into the matrix
    if (!cache.built || !LHS.is_sparse() || !LHS.is_pattern_locked()) {
        bsim4applyStamps(instance, s, LHS, RHS);
        return;
    }

    // Node indexes in mna matrix
    auto get_NodeIndex = [](int node) -> int {
        return node - 1;
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

    const std::array<bool,12> &nodeValid = instance.BSIM4nodeValid;

    // RHS Stamps
    StampRHS(RHS, dNodePrime, s.RHSdNodePrime, nodeValid, BSIM4V82::D_NODE_PRIME);
    StampRHS(RHS, gNodePrime, s.RHSgNodePrime, nodeValid, BSIM4V82::G_NODE_PRIME);
    StampRHS(RHS, gNodeExt, s.RHSgNodeExt, nodeValid, BSIM4V82::G_NODE_EXT);
    StampRHS(RHS, gNodeMid, s.RHSgNodeMid, nodeValid, BSIM4V82::G_NODE_MID);
    StampRHS(RHS, bNodePrime, s.RHSbNodePrime, nodeValid, BSIM4V82::B_NODE_PRIME);
    StampRHS(RHS, sNodePrime, s.RHSsNodePrime, nodeValid, BSIM4V82::S_NODE_PRIME);
    StampRHS(RHS, dbNode, s.RHSdbNode, nodeValid, BSIM4V82::DB_NODE);
    StampRHS(RHS, sbNode, s.RHSsbNode, nodeValid, BSIM4V82::SB_NODE);
    StampRHS(RHS, dNode, s.RHSdNode, nodeValid, BSIM4V82::D_NODE);
    StampRHS(RHS, sNode, s.RHSsNode, nodeValid, BSIM4V82::S_NODE);
    StampRHS(RHS, qNode, s.RHSqNode, nodeValid, BSIM4V82::Q_NODE);

    arma::sp_mat &sp = LHS.get_sparse();
    double *vals = arma::access::rwp(sp.values);

#define BSIM4_APPLY_LHS_CACHED(field, rowNode, colNode, rowType, colType) \
    do {                                                                 \
        if (cache.field != BSIM4StampIndexCache::kInvalid) {             \
            vals[cache.field] += s.field;                                \
        }                                                                \
    } while (0);
    BSIM4_LHS_STAMP_LIST(BSIM4_APPLY_LHS_CACHED)
#undef BSIM4_APPLY_LHS_CACHED
}

} // namespace bsim4

#undef BSIM4_LHS_STAMP_LIST
