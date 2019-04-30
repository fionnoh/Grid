/*************************************************************************************
 
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: Hadrons/Modules/MDistil/Perambulator.hpp
 
 Copyright (C) 2019
 
 Author: Felix Erben <ferben@ed.ac.uk>
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
 See the full license in the file "LICENSE" in the top level distribution directory
 *************************************************************************************/
/*  END LEGAL */

#ifndef Hadrons_MDistil_Perambulator_hpp_
#define Hadrons_MDistil_Perambulator_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

// These are members of Distillation
#include <Hadrons/Distil.hpp>

BEGIN_HADRONS_NAMESPACE


/******************************************************************************
 *                             Perambulator                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class PerambulatorPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PerambulatorPar,
		                    std::string, eigenPack,
		                    std::string, noise,
                                    std::string, PerambFileName, //stem!!!
                                    std::string, UniqueIdentifier,
                                    bool, multiFile,
                                    int, nvec,
                                    DistilParameters, Distil,
				    std::string, solver);
};

template <typename FImpl>
class TPerambulator: public Module<PerambulatorPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    // constructor
    TPerambulator(const std::string name);
    // destructor
    virtual ~TPerambulator(void);
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
protected:
    // These variables are created in setup() and freed in Cleanup()
    GridCartesian * grid3d; // Owned by me, so I must delete it
    GridCartesian * grid4d; // Owned by environment (so I won't delete it)
protected:
    virtual void Cleanup(void);
private:
        unsigned int Ls_;
};

// Can't name the module Perambulator, because that's what we've called the object
MODULE_REGISTER_TMP(Peramb, TPerambulator<FIMPL>, MDistil);

/******************************************************************************
 *                 TPerambulator implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPerambulator<FImpl>::TPerambulator(const std::string name)
: grid3d{nullptr}, grid4d{nullptr}, Module<PerambulatorPar>(name)
{}

// destructor
template <typename FImpl>
TPerambulator<FImpl>::~TPerambulator(void)
{
  Cleanup();
};

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPerambulator<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    in.push_back(par().eigenPack);
    in.push_back(par().solver);
    in.push_back(par().noise);
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TPerambulator<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(),getName() + "_unsmeared_sink"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPerambulator<FImpl>::setup(void)
{
  Cleanup();
  grid4d = env().getGrid();
  grid3d = MakeLowerDimGrid(grid4d);
  const int Nt{grid4d->GlobalDimensions()[Tdir]};
  const int nvec{par().nvec};
  const DistilParameters & Distil{par().Distil};
  const int LI{Distil.LI};
  const int nnoise{Distil.nnoise};
  const int Nt_inv{Distil.Nt_inv}; // TODO: PROBABLY BETTER: if (full_tdil) Nt_inv=1; else Nt_inv = TI;
  const int Ns{Distil.Ns};
  std::array<std::string,6> sIndexNames{"Nt", "nvec", "LI", "nnoise", "Nt_inv", "SI"};
  
  envCreate(Perambulator<SpinVector COMMA 6 COMMA sizeof(Real)>, getName(), 1,
            sIndexNames,Nt,nvec,Distil.LI,Distil.nnoise,Distil.Nt_inv,Distil.SI);
  envCreate(std::vector<FermionField>, getName() + "_unsmeared_sink", 1,
            nnoise*LI*Ns*Nt_inv, envGetGrid(FermionField));
  
  envTmpLat(LatticeSpinColourVector, "dist_source");
  envTmpLat(LatticeSpinColourVector, "tmp2");
  envTmpLat(LatticeSpinColourVector, "result");
  envTmpLat(LatticeColourVector, "result_nospin");
  envTmp(LatticeSpinColourVector, "tmp3d",1,LatticeSpinColourVector(grid3d));
  envTmp(LatticeColourVector, "tmp3d_nospin",1,LatticeColourVector(grid3d));
  envTmp(LatticeColourVector, "result_3d",1,LatticeColourVector(grid3d));
  envTmp(LatticeColourVector, "evec3d",1,LatticeColourVector(grid3d));
  
  Ls_ = env().getObjectLs(par().solver);
  envTmpLat(FermionField, "v4dtmp");
  envTmpLat(FermionField, "v5dtmp", Ls_);
  envTmpLat(FermionField, "v5dtmp_sol", Ls_);
}

// clean up any temporaries created by setup (that aren't stored in the environment)
template <typename FImpl>
void TPerambulator<FImpl>::Cleanup(void)
{
  if( grid3d != nullptr ) {
    delete grid3d;
    grid3d = nullptr;
  }
  grid4d = nullptr;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPerambulator<FImpl>::execute(void)
{
  const int Nt{grid4d->GlobalDimensions()[Tdir]};
    const int nvec{par().nvec};
    const DistilParameters & Distil{par().Distil};
    const int LI{Distil.LI};
    const int SI{Distil.SI};
    const int TI{Distil.TI};
    const int nnoise{Distil.nnoise};
    const int tsrc{Distil.tsrc};
    const int Ns{Distil.Ns};
  const bool full_tdil{TI==Nt};
  const bool exact_distillation{full_tdil && LI==nvec};
  const int Nt_inv{full_tdil ? 1 : TI}; // TODO: PROBABLY BETTER: if (full_tdil) Nt_inv=1; else Nt_inv = TI;

    auto &solver=envGet(Solver, par().solver);
    auto &mat = solver.getFMat();
    envGetTmp(FermionField, v4dtmp);
    envGetTmp(FermionField, v5dtmp);
    envGetTmp(FermionField, v5dtmp_sol);

    const std::string &UniqueIdentifier{par().UniqueIdentifier};

    //auto &noise = envGet(std::vector<Complex>, par().noise);
    auto &noise = envGet(NoiseTensor, par().noise);
    auto &perambulator = envGet(Perambulator<SpinVector COMMA 6 COMMA sizeof(Real)>, getName());
    auto &epack = envGet(Grid::Hadrons::EigenPack<LatticeColourVector>, par().eigenPack);
    auto &unsmeared_sink = envGet(std::vector<FermionField>, getName() + "_unsmeared_sink");


    // Load perambulator if it exists on disk instead of creating it
    // Not sure this is how we want it - rather specify an input flag 'read' 
    // and assert that the file is there.

  envGetTmp(LatticeSpinColourVector, dist_source);
  envGetTmp(LatticeSpinColourVector, tmp2);
  envGetTmp(LatticeSpinColourVector, result);
  envGetTmp(LatticeColourVector, result_nospin);
  envGetTmp(LatticeSpinColourVector, tmp3d);
  envGetTmp(LatticeColourVector, tmp3d_nospin);
  envGetTmp(LatticeColourVector, result_3d);
  envGetTmp(LatticeColourVector, evec3d);

    const int Ntlocal{grid4d->LocalDimensions()[3]};
    const int Ntfirst{grid4d->LocalStarts()[3]};

    {

    int t_inv;
    for (int inoise = 0; inoise < nnoise; inoise++) {
      for (int dk = 0; dk < LI; dk++) {
        for (int dt = 0; dt < Nt_inv; dt++) {
          for (int ds = 0; ds < SI; ds++) {
            std::cout <<  "LapH source vector from noise " << inoise << " and dilution component (d_k,d_t,d_alpha) : (" << dk << ","<< dt << "," << ds << ")" << std::endl;
            dist_source = zero;
            tmp3d_nospin = zero;
            evec3d = zero;
            for (int it = dt; it < Nt; it += TI){
              if (full_tdil) t_inv = tsrc; else t_inv = it;
              if( t_inv >= Ntfirst && t_inv < Ntfirst + Ntlocal ) {
                for (int ik = dk; ik < nvec; ik += LI){
                  for (int is = ds; is < Ns; is += SI){ 
                    ExtractSliceLocal(evec3d,epack.evec[ik],0,t_inv,3);
                    //tmp3d_nospin = evec3d * noise[inoise + nnoise*(t_inv + Nt*(ik+nvec*is))];
                    tmp3d_nospin = evec3d * noise(inoise, t_inv, ik, is);
                    tmp3d=zero;
                    pokeSpin(tmp3d,tmp3d_nospin,is);
                    tmp2=zero;
                    InsertSliceLocal(tmp3d,tmp2,0,t_inv-Ntfirst,Grid::QCD::Tdir);
                    dist_source += tmp2;
                  }
                }
              }
            }
            result=zero;
	    v4dtmp = dist_source;
	    if (Ls_ == 1){
	      solver(result, v4dtmp);
	    } else {
	       mat.ImportPhysicalFermionSource(v4dtmp, v5dtmp);
	       solver(v5dtmp_sol, v5dtmp);
	       mat.ExportPhysicalFermionSolution(v5dtmp_sol, v4dtmp);
	       result = v4dtmp;
	    }
            if ((1)) // comment out if unsmeared sink is too large???
              unsmeared_sink[inoise+nnoise*(dk+LI*(dt+Nt_inv*ds))] = result;
            for (int is = 0; is < Ns; is++) {
              result_nospin = peekSpin(result,is);
              for (int t = Ntfirst; t < Ntfirst + Ntlocal; t++) {
                ExtractSliceLocal(result_3d,result_nospin,0,t-Ntfirst,Grid::QCD::Tdir);
                for (int ivec = 0; ivec < nvec; ivec++) {
                  ExtractSliceLocal(evec3d,epack.evec[ivec],0,t,3);
                  pokeSpin(perambulator(t, ivec, dk, inoise,dt,ds),static_cast<Complex>(innerProduct(evec3d, result_3d)),is);
                }
              }
            }
          }
        }
      }
    }
    }
  std::cout <<  "perambulator done" << std::endl;
  perambulator.SliceShare( grid3d, grid4d );

  if(grid4d->IsBoss()) {
    std::string sPerambName{par().PerambFileName};
    if( sPerambName.length() == 0 )
      sPerambName = getName();
    sPerambName.append( "." );
    sPerambName.append( std::to_string(vm().getTrajectory()));
    //perambulator.WriteBinary(sPerambName);
    perambulator.write(sPerambName.c_str());
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_Perambulator_hpp_
