//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particles.cpp
//  \brief implements functions in particle classes

#include <sstream>
#include <string>
#include "../athena_arrays.hpp"
#include "particles.hpp"

bool Particles::initialized = false;
int Particles::nint = 0;
int Particles::nreal = 0;
int Particles::ipid = -1;
int Particles::ixp1 = -1, Particles::ixp2 = -1, Particles::ixp3 = -1;
int Particles::ivp1 = -1, Particles::ivp2 = -1, Particles::ivp3 = -1;

void _ErrorIfInitialized(const std::string& calling_function, bool initialized);

//--------------------------------------------------------------------------------------
//! \fn Particles::Initialize()
//  \brief initializes the class.

void Particles::Initialize()
{
  // Add particle ID.
  ipid = AddIntProperty();

  // Add particle position.
  ixp1 = AddRealProperty();
  ixp2 = AddRealProperty();
  ixp3 = AddRealProperty();

  // Add particle velocity.
  ivp1 = AddRealProperty();
  ivp2 = AddRealProperty();
  ivp3 = AddRealProperty();

  initialized = true;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
//  \brief constructs a Particles instance.

Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
{
  // Initialize the class if first called.
  if (!initialized) Particles::Initialize();

  // Point to the calling MeshBlock.
  pmy_block = pmb;
  npar = 1;
  nparmax = 1;  // TODO: dynamically adjust npar and nparmax.

  // Allocate integer properties.
  intprop.NewAthenaArray(nint,nparmax);
  pid.InitWithShallowSlice(intprop, 1, ipid, 1);

  // Allocate integer properties.
  realprop.NewAthenaArray(nreal,nparmax);
  xp1.InitWithShallowSlice(realprop, 1, ixp1, 1);
  xp2.InitWithShallowSlice(realprop, 1, ixp2, 1);
  xp3.InitWithShallowSlice(realprop, 1, ixp3, 1);
  vp1.InitWithShallowSlice(realprop, 1, ivp1, 1);
  vp2.InitWithShallowSlice(realprop, 1, ivp2, 1);
  vp3.InitWithShallowSlice(realprop, 1, ivp3, 1);
}

//--------------------------------------------------------------------------------------
//! \fn Particles::~Particles()
//  \brief destroys a Particles instance.

Particles::~Particles()
{
  // Delete integer properties.
  intprop.DeleteAthenaArray();

  // Delete real properties.
  realprop.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Drift(Real t, Real dt)
//  \brief updates the particle positions from t to t + dt given velocities at t.

void Particles::Drift(Real t, Real dt)
{
  for (long k = 0; k < npar; ++k) {
    xp1(k) += dt * vp1(k);
    xp2(k) += dt * vp2(k);
    xp3(k) += dt * vp3(k);
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Kick(Real t, Real dt)
//  \brief updates the particle velocities from t to t + dt given accelerations at t.

#include <cmath>

void Particles::Kick(Real t, Real dt)
{
  Real a1 = 0.0, a2 = 0.0, a3 = 0.0;  // TODO: might need to be put inside the loop for
                                      // vectorization.

  for (long k = 0; k < npar; ++k) {
    vp1(k) += dt * a1;
    vp2(k) += dt * a2;
    vp3(k) += dt * a3;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Update(Mesh *pm)
//  \brief updates all particle positions and velocities from t to t + dt.
//======================================================================================
#include <iostream>

void Particles::Update(Mesh *pm)
{
  MeshBlock *pmb = pm->pblock;
  Real dth = 0.5 * pm->dt;

  // Loop over MeshBlocks
  while (pmb != NULL) {
    pmb->ppar->Drift(pm->time, dth);
    pmb->ppar->Kick(pm->time, pm->dt);
    pmb->ppar->Drift(pm->time + dth, dth);
    pmb = pmb->next;
  }
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddIntProperty()
//  \brief adds one integer property to the particles and returns the index.

int Particles::AddIntProperty()
{
  _ErrorIfInitialized("Particles::AddIntProperty", initialized);
  return nint++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddRealProperty()
//  \brief adds one real property to the particles and returns the index.

int Particles::AddRealProperty()
{
  _ErrorIfInitialized("Particles::AddRealProperty", initialized);
  return nreal++;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::GetSizeInBytes()
//  \brief returns the data size in bytes in the meshblock.

size_t Particles::GetSizeInBytes()
{
  size_t size = sizeof(npar);
  if (npar > 0) size += npar * (nint * sizeof(long) + nreal * sizeof(Real));
  return size;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::ReadRestart()
//  \brief reads the particle data from the restart file.

#include <cstring>

void Particles::ReadRestart(char *mbdata, int &os)
{
  // Read number of particles.
  std::memcpy(&npar, &(mbdata[os]), sizeof(npar));
  os += sizeof(npar);
  if (npar > nparmax)
  {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::ReadRestart]" << std::endl
        << "npar = " << npar << " > nparmax = " << nparmax << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  if (npar > 0) {
    // Read integer properties.
    size_t size = npar * sizeof(long);
    for (int k = 0; k < nint; ++k) {
      std::memcpy(&(intprop(k,0)), &(mbdata[os]), size);
      os += size;
    }

    // Read real properties.
    size = npar * sizeof(Real);
    for (int k = 0; k < nreal; ++k) {
      std::memcpy(&(realprop(k,0)), &(mbdata[os]), size);
      os += size;
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::WriteRestart()
//  \brief writes the particle data to the restart file.

void Particles::WriteRestart(char *&pdata)
{
  // Write number of particles.
  memcpy(pdata, &npar, sizeof(npar));
  pdata += sizeof(npar);

  if (npar > 0) {
    // Write integer properties.
    size_t size = npar * sizeof(long);
    for (int k = 0; k < nint; ++k) {
      std::memcpy(pdata, &(intprop(k,0)), size);
      pdata += size;
    }

    // Write real properties.
    size = npar * sizeof(Real);
    for (int k = 0; k < nreal; ++k) {
      std::memcpy(pdata, &(realprop(k,0)), size);
      pdata += size;
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::FormattedTableOutput()
//  \brief outputs the particle data in tabulated format.

#include <fstream>
#include <iomanip>

void Particles::FormattedTableOutput(Mesh *pm, OutputParameters op)
{
  MeshBlock *pmb = pm->pblock;
  Particles *ppar;
  std::stringstream fname, msg;
  std::ofstream os;

  // Loop over MeshBlocks
  while (pmb != NULL) {
    ppar = pmb->ppar;

    // Create the filename.
    fname << op.file_basename
          << ".block" << pmb->gid << '.' << op.file_id
          << '.' << std::setw(5) << std::right << std::setfill('0') << op.file_number
          << '.' << "par.tab";

    // Open the file for write.
    os.open(fname.str());
    if (!os.is_open()) {
      msg << "### FATAL ERROR in function [Particles::FormattedTableOutput]"
          << std::endl << "Output file '" << fname.str() << "' could not be opened"
          << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    // Write the particle data in the meshblock.
    for (long k = 0; k < ppar->npar; ++k)
      os << ppar->pid(k) << "  "
         << ppar->xp1(k) << "  " << ppar->xp2(k) << "  " << ppar->xp3(k) << "  "
         << ppar->vp1(k) << "  " << ppar->vp2(k) << "  " << ppar->vp3(k) << std::endl;

    // Close the file and get the next meshblock.
    os.close();
    pmb = pmb->next;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void _ErrorIfInitialized(const string& calling_function)
//  \brief throws an error when the class has already been initialized.

void _ErrorIfInitialized(const std::string& calling_function, bool initialized)
{
  if (initialized)
  {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [" << calling_function << "]" << std::endl
        << "The Particles class has already been initialized. " << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
}
