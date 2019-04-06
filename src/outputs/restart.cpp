//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file restart.cpp
//  \brief writes restart files

// C headers

// C++ headers
#include <cstdio>    // snprintf()
#include <cstring>   // memcpy()
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../particles/particles.hpp"
#include "outputs.hpp"

//----------------------------------------------------------------------------------------
// RestartOutput constructor
// destructor - not needed for this derived class

RestartOutput::RestartOutput(OutputParameters oparams)
    : OutputType(oparams) {
}

//----------------------------------------------------------------------------------------
//! \fn void RestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag)
//  \brief Cycles over all MeshBlocks and writes data to a single restart file.

void RestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool force_write) {
  MeshBlock *pmb;
  IOWrapper resfile;

  // create single output filename:"file_basename"+"."+XXXXX+".rst",
  // where XXXXX = 5-digit file_number
  std::string fname;
  char number[6];
  std::snprintf(number, sizeof(number), "%05d", output_params.file_number);

  fname.assign(output_params.file_basename);
  fname.append(".");
  // add file number to name, unless write is forced by terminate signal, in which case
  // replace number in the name by the string "final".  This keeps the restart file
  // numbers consistent with output.dt when a job is restarted many times.
  if (force_write==false)
    fname.append(number);
  else
    fname.append("final");
  fname.append(".rst");

  // increment counters now so values for *next* dump are stored in restart file
  if (force_write==false) {
    output_params.file_number++;
    output_params.next_time += output_params.dt;
    pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
    pin->SetReal(output_params.block_name, "next_time", output_params.next_time);
  }
  resfile.Open(fname.c_str(), IOWrapper::FileMode::write);

  // prepare the input parameters
  std::stringstream ost;
  pin->ParameterDump(ost);
  std::string sbuf=ost.str();

  // calculate the header size
  IOWrapperSizeT udsize=0;
  for (int n=0; n<pm->nint_user_mesh_data_; n++)
    udsize+=pm->iuser_mesh_data[n].GetSizeInBytes();
  for (int n=0; n<pm->nreal_user_mesh_data_; n++)
    udsize+=pm->ruser_mesh_data[n].GetSizeInBytes();

  IOWrapperSizeT headeroffset=sbuf.size()*sizeof(char)+3*sizeof(int)
                 +sizeof(RegionSize)+2*sizeof(Real)+udsize;
  // the size of each MeshBlock
  int nbtotal=pm->nbtotal;
  int myns=pm->nslist[Globals::my_rank];
  int mynb=pm->nblist[Globals::my_rank];

  // construct the size and offset lists
  IOWrapperSizeT *sizelist = new IOWrapperSizeT[nbtotal];
  pmb=pm->pblock;
  while (pmb!=nullptr) {
    sizelist[pmb->gid]=pmb->GetBlockSizeInBytes();
    pmb=pmb->next;
  }
#ifdef MPI_PARALLEL
  // collect the size list - assuming IOWrapperSizeT is 64bit integer
  MPI_Allgatherv(MPI_IN_PLACE, mynb, MPI_UINT64_T, sizelist, pm->nblist, pm->nslist,
                 MPI_UINT64_T, MPI_COMM_WORLD);
#endif

  // write the header; this part is serial
  if (Globals::my_rank==0) {
    // output the input parameters
    resfile.Write(sbuf.c_str(),sizeof(char),sbuf.size());

    // output Mesh information
    resfile.Write(&(pm->nbtotal), sizeof(int), 1);
    resfile.Write(&(pm->root_level), sizeof(int), 1);
    resfile.Write(&(pm->mesh_size), sizeof(RegionSize), 1);
    resfile.Write(&(pm->time), sizeof(Real), 1);
    resfile.Write(&(pm->dt), sizeof(Real), 1);
    resfile.Write(&(pm->ncycle), sizeof(int), 1);

    // collect and write user Mesh data
    if (udsize!=0) {
      char *ud = new char[udsize];
      IOWrapperSizeT udoffset = 0;
      for (int n=0; n<pm->nint_user_mesh_data_; n++) {
        std::memcpy(&(ud[udoffset]), pm->iuser_mesh_data[n].data(),
                    pm->iuser_mesh_data[n].GetSizeInBytes());
        udoffset+=pm->iuser_mesh_data[n].GetSizeInBytes();
      }
      for (int n=0; n<pm->nreal_user_mesh_data_; n++) {
        std::memcpy(&(ud[udoffset]), pm->ruser_mesh_data[n].data(),
                    pm->ruser_mesh_data[n].GetSizeInBytes());
        udoffset+=pm->ruser_mesh_data[n].GetSizeInBytes();
      }
      resfile.Write(ud, 1, udsize);
      delete [] ud;
    }
  }

  // allocate memory for the metadata of each MeshBlock
  IOWrapperSizeT mdsize=sizeof(LogicalLocation)+sizeof(Real)+sizeof(IOWrapperSizeT);
  char *mbmetadata = new char[mdsize*mynb];

  // Loop over MeshBlocks and pack the meta data
  pmb=pm->pblock;
  int os=0;
  while (pmb!=nullptr) {
    std::memcpy(&(mbmetadata[os]), &(pmb->loc), sizeof(LogicalLocation));
    os+=sizeof(LogicalLocation);
    std::memcpy(&(mbmetadata[os]), &(pmb->cost), sizeof(Real));
    os+=sizeof(Real);
    std::memcpy(&(mbmetadata[os]), &(sizelist[pmb->gid]), sizeof(IOWrapperSizeT));
    os+=sizeof(IOWrapperSizeT);
    pmb=pmb->next;
  }

  // write the ID list collectively
  IOWrapperSizeT myoffset=headeroffset+mdsize*myns;
  resfile.Write_at_all(mbmetadata,mdsize,mynb,myoffset);

  delete [] mbmetadata;

  // calculate the size and offset for this rank
  IOWrapperSizeT mysize=0;
  for (int n=myns; n<myns+mynb; n++)
    mysize+=sizelist[n];
  myoffset=headeroffset+mdsize*nbtotal;
  for (int n=0; n<myns; n++)
    myoffset+=sizelist[n];

  delete [] sizelist;

  // allocate memory for the output of this rank
  char *data = new char[mysize];
  char *pdata=data;

  // loop over MeshBlocks and pack the data
  pmb=pm->pblock;
  while (pmb != nullptr) {
    std::memcpy(pdata,pmb->phydro->u.data(), pmb->phydro->u.GetSizeInBytes());
    pdata+=pmb->phydro->u.GetSizeInBytes();
    if (GENERAL_RELATIVITY) {
      std::memcpy(pdata,pmb->phydro->w.data(), pmb->phydro->w.GetSizeInBytes());
      pdata+=pmb->phydro->w.GetSizeInBytes();
      std::memcpy(pdata,pmb->phydro->w1.data(), pmb->phydro->w1.GetSizeInBytes());
      pdata+=pmb->phydro->w1.GetSizeInBytes();
    }
    if (MAGNETIC_FIELDS_ENABLED) {
      std::memcpy(pdata,pmb->pfield->b.x1f.data(),pmb->pfield->b.x1f.GetSizeInBytes());
      pdata+=pmb->pfield->b.x1f.GetSizeInBytes();
      std::memcpy(pdata,pmb->pfield->b.x2f.data(),pmb->pfield->b.x2f.GetSizeInBytes());
      pdata+=pmb->pfield->b.x2f.GetSizeInBytes();
      std::memcpy(pdata,pmb->pfield->b.x3f.data(),pmb->pfield->b.x3f.GetSizeInBytes());
      pdata+=pmb->pfield->b.x3f.GetSizeInBytes();
    }

    // NEW_PHYSICS: add output of additional physics to restarts here
    // also update MeshBlock::GetBlockSizeInBytes accordingly and MeshBlock constructor
    // for restarts.

    // Self-gravity outout is removed, as the implemented solvers do not need to save phi

    if (PARTICLES) pmb->ppar->PackParticlesForRestart(pdata);

    // pack the user MeshBlock data
    for (int n=0; n<pmb->nint_user_meshblock_data_; n++) {
      std::memcpy(pdata, pmb->iuser_meshblock_data[n].data(),
                  pmb->iuser_meshblock_data[n].GetSizeInBytes());
      pdata+=pmb->iuser_meshblock_data[n].GetSizeInBytes();
    }
    for (int n=0; n<pmb->nreal_user_meshblock_data_; n++) {
      std::memcpy(pdata, pmb->ruser_meshblock_data[n].data(),
                  pmb->ruser_meshblock_data[n].GetSizeInBytes());
      pdata+=pmb->ruser_meshblock_data[n].GetSizeInBytes();
    }

    pmb=pmb->next;
  }

  // now write restart data in parallel
  resfile.Write_at_all(data,mysize,1,myoffset);
  resfile.Close();

  delete [] data;
}
