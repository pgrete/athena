//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file radiation_integrator.cpp
//  \brief derived class for radiation integrator task list.
//======================================================================================

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#include "../defs.hpp"
#include "../radiation/integrators/rad_integrators.hpp"
#include "../radiation/radiation.hpp"
#include "../utils/cgk_utils.hpp"


// this class header
#include "task_list.hpp"

//--------------------------------------------------------------------------------------
//  RadiationIntegratorTaskList constructor
RadiationIntegratorTaskList::RadiationIntegratorTaskList(ParameterInput *pin, Mesh *pm)
  : TaskList(pm)
{
  nstages = 1;
  integrator = RADIATION_INTEGRATOR;
  // Now assemble list of tasks for each step of chemistry integrator
  {using namespace RadiationIntegratorTaskNames;
    if (integrator == "loc_jeans") {
      AddRadiationIntegratorTask(INT_LOC_JEANS,NONE);
    } else if (integrator == "six_ray") {
      //add six ray
      AddRadiationIntegratorTask(START_SIXRAY_RECV,NONE);
      AddRadiationIntegratorTask(GET_COL_MB0,NONE);
      AddRadiationIntegratorTask(RECV_SEND_COL0,GET_COL_MB0);
      AddRadiationIntegratorTask(GET_COL_MB1,NONE);
      AddRadiationIntegratorTask(RECV_SEND_COL1,GET_COL_MB1);
      AddRadiationIntegratorTask(GET_COL_MB2,NONE);
      AddRadiationIntegratorTask(RECV_SEND_COL2,GET_COL_MB2);
      AddRadiationIntegratorTask(GET_COL_MB3,NONE);
      AddRadiationIntegratorTask(RECV_SEND_COL3,GET_COL_MB3);
      AddRadiationIntegratorTask(GET_COL_MB4,NONE);
      AddRadiationIntegratorTask(RECV_SEND_COL4,GET_COL_MB4);
      AddRadiationIntegratorTask(GET_COL_MB5,NONE);
      AddRadiationIntegratorTask(RECV_SEND_COL5,GET_COL_MB5);
      AddRadiationIntegratorTask(UPDATE_RAD,
          RECV_SEND_COL0|RECV_SEND_COL1|RECV_SEND_COL2|
          RECV_SEND_COL3|RECV_SEND_COL4|RECV_SEND_COL5);
      AddRadiationIntegratorTask(CLEAR_SIXRAY_RECV,UPDATE_RAD);
    } else if (integrator == "const") {
      //do nothing, radiation field constant, remain initial value
      AddRadiationIntegratorTask(INT_CONST,NONE);
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in Radiation task list" << std::endl
        << "integrator=" << integrator << " not valid radiation integrator, " << std::endl
        << "choose from {jeans, six_ray, const}" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  } // end of using namespace block
}

//--------------------------------------------------------------------------------------
//! \fn
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  
void RadiationIntegratorTaskList::AddRadiationIntegratorTask(uint64_t id, uint64_t dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;
  using namespace RadiationIntegratorTaskNames;
  switch((id)) {
    case (INT_LOC_JEANS):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::LocalIntegratorJeans);
      break;
    case (INT_CONST):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::ConstRadiation);
      break;
    //add six ray here
    case (START_SIXRAY_RECV):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::StartSixrayReceive);
      break;
    case (GET_COL_MB0):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB0);
      break;
    case (GET_COL_MB1):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB1);
      break;
    case (GET_COL_MB2):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB2);
      break;
    case (GET_COL_MB3):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB3);
      break;
    case (GET_COL_MB4):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB4);
      break;
    case (GET_COL_MB5):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB5);
      break;
    case (RECV_SEND_COL0):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend0);
      break;
    case (RECV_SEND_COL1):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend1);
      break;
    case (RECV_SEND_COL2):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend2);
      break;
    case (RECV_SEND_COL3):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend3);
      break;
    case (RECV_SEND_COL4):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend4);
      break;
    case (RECV_SEND_COL5):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend5);
      break;
    case (UPDATE_RAD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::UpdateRadiation);
      break;
    case (CLEAR_SIXRAY_RECV):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::ClearSixrayReceive);
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in Add Radiation Task" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}

//--------------------------------------------------------------------------------------
// Functions to integrate Radiation
enum TaskStatus RadiationIntegratorTaskList::LocalIntegratorJeans(MeshBlock *pmb,
                                                                 int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->UpdateRadiation(0);
  pmb->prad->pradintegrator->CopyToOutput();
#endif
  return TASK_SUCCESS;
}

enum TaskStatus RadiationIntegratorTaskList::ConstRadiation(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->CopyToOutput();
#endif
  return TASK_SUCCESS;
}

//add six ray here

enum TaskStatus RadiationIntegratorTaskList::StartSixrayReceive(MeshBlock *pmb, 
                                                                int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->pbval->StartReceivingSixray();
#endif
  return TASK_SUCCESS;
}

enum TaskStatus RadiationIntegratorTaskList::GetColMB0(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(0);
#endif
  return TASK_SUCCESS;
}

enum TaskStatus RadiationIntegratorTaskList::GetColMB1(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(1);
#endif
  return TASK_SUCCESS;
}

enum TaskStatus RadiationIntegratorTaskList::GetColMB2(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(2);
#endif
  return TASK_SUCCESS;
}

enum TaskStatus RadiationIntegratorTaskList::GetColMB3(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(3);
#endif
  return TASK_SUCCESS;
}

enum TaskStatus RadiationIntegratorTaskList::GetColMB4(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(4);
#endif
  return TASK_SUCCESS;
}

enum TaskStatus RadiationIntegratorTaskList::GetColMB5(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(5);
#endif
  return TASK_SUCCESS;
}

enum TaskStatus RadiationIntegratorTaskList::RecvAndSend_direction(
    MeshBlock *pmb, int step, int direction)
{
#ifdef INCLUDE_CHEMISTRY
  int opp_direction = CGKUtility::GetOppositeDirection(direction);
  bool ret = true;
  NeighborBlock *nb = pmb->prad->pradintegrator->pfacenb_[direction];
  NeighborBlock *nb_opp = pmb->prad->pradintegrator->pfacenb_[opp_direction];
  if (nb == NULL) {
    if (nb_opp != NULL) {
      pmb->pbval->SendSixrayBoundaryBuffers(pmb->prad->pradintegrator->col, opp_direction);
    }
  } else {
    ret = pmb->pbval->ReceiveSixrayBoundaryBuffers(pmb->prad->pradintegrator->col,
                                                  direction);
    if (ret == true) {
      pmb->prad->pradintegrator->UpdateCol(direction);
      if (nb_opp != NULL) {
        pmb->pbval->SendSixrayBoundaryBuffers(pmb->prad->pradintegrator->col,
                                              opp_direction);
      }
    } else {
      return TASK_FAIL;
    }
  }
#endif
  return TASK_SUCCESS;
}

enum TaskStatus RadiationIntegratorTaskList::RecvAndSend0(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, 0);
}

enum TaskStatus RadiationIntegratorTaskList::RecvAndSend1(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, 1);
}

enum TaskStatus RadiationIntegratorTaskList::RecvAndSend2(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, 2);
}

enum TaskStatus RadiationIntegratorTaskList::RecvAndSend3(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, 3);
}

enum TaskStatus RadiationIntegratorTaskList::RecvAndSend4(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, 4);
}

enum TaskStatus RadiationIntegratorTaskList::RecvAndSend5(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, 5);
}

enum TaskStatus RadiationIntegratorTaskList::UpdateRadiation(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  for (int i=0; i<6; i++) {
    pmb->prad->pradintegrator->UpdateRadiation(i);
  }
  pmb->prad->pradintegrator->CopyToOutput();
#endif
  return TASK_SUCCESS;
}

enum TaskStatus RadiationIntegratorTaskList::ClearSixrayReceive(MeshBlock *pmb,
                                                                 int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->pbval->ClearBoundarySixray();
#endif
  return TASK_SUCCESS;
}
