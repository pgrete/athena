
#ifndef TASK_LIST_RADIATION_TASK_LIST_HPP_
#define TASK_LIST_RADIATION_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file radiation_task_list.hpp
//    \brief provides functionality to control dynamic execution using tasks

// C headers

// C++ headers
#include <cstdint>      // std::uint64_t

// Athena++ headers
#include "../athena.hpp"
#include "task_list.hpp"

// forward declarations
class Mesh;
class MeshBlock;

//--------------------------------------------------------------------------------------
//! \class RadiationIntegratorTaskList
//  \brief data and function definitions for RadiationIntegratorTaskList derived class
//
class RadiationIntegratorTaskList : public TaskList {
 public:
  RadiationIntegratorTaskList(ParameterInput *pin, Mesh *pm);
  ~RadiationIntegratorTaskList() {};
  std::string integrator;
  void AddRadiationIntegratorTask(uint64_t id, uint64_t dep);
  enum TaskStatus LocalIntegratorJeans(MeshBlock *pmb, int step);
  enum TaskStatus ConstRadiation(MeshBlock *pmb, int step);
  enum TaskStatus StartSixrayReceive(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB0(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB1(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB2(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB3(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB4(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB5(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend0(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend1(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend2(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend3(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend4(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend5(MeshBlock *pmb, int step);
  enum TaskStatus UpdateRadiation(MeshBlock *pmb, int step);
  enum TaskStatus ClearSixrayReceive(MeshBlock *pmb, int step);
 private:
  enum TaskStatus RecvAndSend_direction(MeshBlock *pmb, int step, int direction);
  void StartupTaskList(MeshBlock *pmb, int stage) override;
};

namespace RadiationIntegratorTaskNames {
const std::uint64_t NONE            = 0ULL;
const std::uint64_t INT_CONST       = 1ULL<<1; //constant radiation, do nothing
const std::uint64_t INT_LOC_JEANS   = 1ULL<<2; //local jeans shielding
const std::uint64_t CLEAR_SIXRAY_RECV=1ULL<<3; 
const std::uint64_t GET_COL_MB0     = 1ULL<<4; 
const std::uint64_t GET_COL_MB1     = 1ULL<<5; 
const std::uint64_t GET_COL_MB2     = 1ULL<<6;
const std::uint64_t GET_COL_MB3     = 1ULL<<7;
const std::uint64_t GET_COL_MB4     = 1ULL<<8;
const std::uint64_t GET_COL_MB5     = 1ULL<<9;
const std::uint64_t RECV_SEND_COL0  = 1ULL<<10; 
const std::uint64_t RECV_SEND_COL1  = 1ULL<<11;
const std::uint64_t RECV_SEND_COL2  = 1ULL<<12;
const std::uint64_t RECV_SEND_COL3  = 1ULL<<13;
const std::uint64_t RECV_SEND_COL4  = 1ULL<<14;
const std::uint64_t RECV_SEND_COL5  = 1ULL<<15;
const std::uint64_t UPDATE_RAD      = 1ULL<<16;
};

#endif // TASK_LIST_RADIATION_TASK_LIST_HPP_
