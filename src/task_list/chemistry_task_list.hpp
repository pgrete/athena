#ifndef TASK_LIST_CHEMISTRY_TASK_LIST_HPP_
#define TASK_LIST_CHEMISTRY_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file chemistry_task_list.hpp
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
//! \class ChemistryIntegratorTaskList
//  \brief data and function definitions for ChemistryIntegratorTaskList derived class
//
class ChemistryIntegratorTaskList : public TaskList {
 public:
  ChemistryIntegratorTaskList(ParameterInput *pin, Mesh *pm);
  ~ChemistryIntegratorTaskList() {};
  void AddChemistryIntegratorTask(uint64_t id, uint64_t dep);
  enum TaskStatus IntegrateSourceTerm(MeshBlock *pmb, int step);
  enum TaskStatus StartSpeciesReceive(MeshBlock *pmb, int step);
  enum TaskStatus ClearSpeciesReceive(MeshBlock *pmb, int step);
  enum TaskStatus SpeciesSend(MeshBlock *pmb, int step);
  enum TaskStatus SpeciesReceive(MeshBlock *pmb, int step);
  //add advection term here. 
 private:
  void StartupTaskList(MeshBlock *pmb, int stage) override;
};


namespace ChemistryIntegratorTaskNames {
const std::uint64_t NONE            = 0ULL;
const std::uint64_t INT_CHEM_SRC    = 1ULL<<1; //chemistry source term
const std::uint64_t CLEAR_SPEC_RECV = 1ULL<<2; 
const std::uint64_t SEND_SPEC       = 1ULL<<3; 
const std::uint64_t RECV_SPEC       = 1ULL<<4; 
  //add advection term here
};

#endif // TASK_LIST_CHEMISTRY_TASK_LIST_HPP_
