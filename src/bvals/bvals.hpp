#ifndef BOUNDARY_VALUES_HPP
#define BOUNDARY_VALUES_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file bvals.hpp
//  \brief defines BoundaryValues class used for setting BCs on all data types
//======================================================================================

// C++ headers
#include <string>   // string

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// forward declarations
class MeshBlock;
class Hydro;
class Field;
class ParameterInput;
class Coordinates;
class Radiation;
struct FaceField;
struct NeighborBlock;
struct PolarNeighborBlock;

// identifiers for all 6 faces of a MeshBlock on which boundary conditions are applied
enum BoundaryFace {FACE_UNDEF=-1, INNER_X1=0, OUTER_X1=1, INNER_X2=2, OUTER_X2=3, 
  INNER_X3=4, OUTER_X3=5};

// identifiers for boundary conditions
enum BoundaryFlag {BLOCK_BNDRY=-1, BNDRY_UNDEF=0, REFLECTING_BNDRY=1, OUTFLOW_BNDRY=2,
  USER_BNDRY=3, PERIODIC_BNDRY=4, POLAR_BNDRY=5, VACUUM_BNDRY=6};

// identifiers for different Cell Center Variables
enum CellCenterPhysics {HYDRO=1, RAD=2};

//-------------------- prototypes for all BC functions ---------------------------------
void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);

void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);

// For radiation boundary condition
void RadReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);

void RadOutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadOutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadOutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadOutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadOutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadOutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);


void RadVacuumInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadVacuumInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadVacuumInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadVacuumOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadVacuumOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);
void RadVacuumOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    int is, int ie, int js, int je, int ks, int ke);


void RotateHPi_InnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke);

void RotateHPi_OuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke);

void RotateHPi_InnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                     int is, int ie, int js, int je, int ks, int ke);

void RotateHPi_OuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                     int is, int ie, int js, int je, int ks, int ke);


void RotatePi_InnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                     int is, int ie, int js, int je, int ks, int ke);

void RotatePi_OuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                     int is, int ie, int js, int je, int ks, int ke);


// function to return boundary flag given input string
enum BoundaryFlag GetBoundaryFlag(std::string input_string);


//! \class BoundaryValues
//  \brief BVals data and functions

class BoundaryValues {
public:
  BoundaryValues(MeshBlock *pmb, ParameterInput *pin);
  ~BoundaryValues();

  void Initialize(void);
  void StartReceivingForInit(void);
  void StartReceivingAll(void);

  void CheckBoundary(void);

  
  void SendCenterBoundaryBuffers(AthenaArray<Real> &src,
        AthenaArray<Real> &src_rad, int step, bool conserved_values);
  
  void LoadCenterBoundaryBufferSameLevel(AthenaArray<Real> &src, Real *buf,
                            const NeighborBlock& nb, int phys, int &p);
  void LoadCenterBoundaryBufferToCoarser(AthenaArray<Real> &src, Real *buf,
               const NeighborBlock& nb, int phys, bool conserved_values, int &p);
  void LoadCenterBoundaryBufferToFiner(AthenaArray<Real> &src, Real *buf,
                            const NeighborBlock& nb, int phys, int &p);
  
  bool ReceiveCenterBoundaryBuffers(AthenaArray<Real> &dst,
                           AthenaArray<Real> &dst_rad, int step);
  void ReceiveCenterBoundaryBuffersWithWait(AthenaArray<Real> &dst,
                           AthenaArray<Real> &dst_rad, int step);
  
  void SetCenterBoundarySameLevel(AthenaArray<Real> &dst,
               AthenaArray<Real> &dst_rad, Real *buf, const NeighborBlock& nb);
  
  void SetCenterBoundaryFromCoarser(Real *buf,
                        const NeighborBlock& nb, bool conserved_values);
  
  void SetCenterBoundaryFromFiner(AthenaArray<Real> &dst,
               AthenaArray<Real> &dst_rad, Real *buf, const NeighborBlock& nb);


  void SendFluxCorrection(int step);
  bool ReceiveFluxCorrection(int step);
  

  int LoadFieldBoundaryBufferSameLevel(FaceField &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadFieldBoundaryBufferToCoarser(FaceField &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadFieldBoundaryBufferToFiner(FaceField &src, Real *buf,
                                     const NeighborBlock& nb);
  void SendFieldBoundaryBuffers(FaceField &src, int step);
  void SetFieldBoundarySameLevel(FaceField &dst, Real *buf, const NeighborBlock& nb);
  void SetFieldBoundaryFromCoarser(Real *buf, const NeighborBlock& nb);
  void SetFieldBoundaryFromFiner(FaceField &dst, Real *buf, const NeighborBlock& nb);
  bool ReceiveFieldBoundaryBuffers(FaceField &dst, int step);
  void ReceiveFieldBoundaryBuffersWithWait(FaceField &dst, int step);

  int LoadEMFBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryPolarBuffer(Real *buf, const PolarNeighborBlock &nb);
  void SendEMFCorrection(int step);
  void SetEMFBoundarySameLevel(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryFromFiner(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryPolar(Real **buf_list, int num_bufs, bool north);
  void ClearCoarseEMFBoundary(void);
  void AverageEMFBoundary(void);
  void PolarSingleEMF(void);
  void PolarSingleHydro(AthenaArray<Real> &dst);
  void PolarSingleRad(AthenaArray<Real> &dst);
  void PolarSingleField(FaceField &dst);
  bool ReceiveEMFCorrection(int step);

  void ProlongateBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst, 
                            FaceField &bfdst, AthenaArray<Real> &bcdst);
  void ProlongateRadBoundaries(AthenaArray<Real> &dst);


  void ApplyPhysicalBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst,
                               FaceField &bfdst, AthenaArray<Real> &bcdst);

  void ApplyRadPhysicalBoundaries(AthenaArray<Real> &pdst);

  void ClearBoundaryForInit(void);
  void ClearBoundaryAll(void);

private:
  MeshBlock *pmy_mblock_;  // ptr to MeshBlock containing this BVals

  BValFunc_t BoundaryFunction_[6];
  
  RadBValFunc_t RadBoundaryFunction_[6]; // Function Pointer for radiation

  int nface_, nedge_;
  bool edge_flag_[12];
  int nedge_fine_[12];
  bool firsttime_[NSTEP];

  // cc for cell_center,
  // MPI communication once for hydro and radiation boundary
  enum boundary_status cc_flag_[NSTEP][56], field_flag_[NSTEP][56];
  enum boundary_status cc_flcor_flag_[NSTEP][6][2][2];
  enum boundary_status emfcor_flag_[NSTEP][48];
  enum boundary_status *emf_north_flag_[NSTEP];
  enum boundary_status *emf_south_flag_[NSTEP];
  Real *cc_send_[NSTEP][56],  *cc_recv_[NSTEP][56];
  Real *field_send_[NSTEP][56],  *field_recv_[NSTEP][56];
  Real *cc_flcor_send_[NSTEP][6],   *cc_flcor_recv_[NSTEP][6][2][2];
  Real *emfcor_send_[NSTEP][48], *emfcor_recv_[NSTEP][48];
  Real **emf_north_send_[NSTEP], **emf_north_recv_[NSTEP];
  Real **emf_south_send_[NSTEP], **emf_south_recv_[NSTEP];
  AthenaArray<Real> sarea_[2];
  AthenaArray<Real> exc_, rad_exc_;
  int num_north_polar_blocks_, num_south_polar_blocks_;

#ifdef MPI_PARALLEL
  MPI_Request req_cc_send_[NSTEP][56],  req_cc_recv_[NSTEP][56];
  MPI_Request req_field_send_[NSTEP][56],  req_field_recv_[NSTEP][56];
  MPI_Request req_cc_flcor_send_[NSTEP][6],   req_cc_flcor_recv_[NSTEP][6][2][2];
  MPI_Request req_emfcor_send_[NSTEP][48], req_emfcor_recv_[NSTEP][48];
  MPI_Request *req_emf_north_send_[NSTEP], *req_emf_north_recv_[NSTEP];
  MPI_Request *req_emf_south_send_[NSTEP], *req_emf_south_recv_[NSTEP];
#endif
  // temporary
  friend class Mesh;
};

unsigned int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2);
int BufferID(int dim, bool multilevel, bool face_only);
int FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2, int bmax);

typedef struct NeighborIndexes
{
  int ox1, ox2, ox3, fi1, fi2;
  enum neighbor_type type;
} NeighborIndexes;


#endif // BOUNDARY_VALUES_HPP
