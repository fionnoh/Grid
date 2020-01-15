#ifndef Hadrons_MNoise_SparseGridSpinColorDiagonal_hpp_
#define Hadrons_MNoise_SparseGridSpinColorDiagonal_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SparseGridSpinColorDiagonal                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

class SparseGridSpinColorDiagonalPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SparseGridSpinColorDiagonalPar,
                                    unsigned int, nsrc,
                                    unsigned int, nsparse);
};

template <typename FImpl>
class TSparseGridSpinColorDiagonal: public Module<SparseGridSpinColorDiagonalPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSparseGridSpinColorDiagonal(const std::string name);
    // destructor
    virtual ~TSparseGridSpinColorDiagonal(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SparseGridSpinColorDiagonal, TSparseGridSpinColorDiagonal<FIMPL>, MNoise);
MODULE_REGISTER_TMP(ZSparseGridSpinColorDiagonal, TSparseGridSpinColorDiagonal<ZFIMPL>, MNoise);

/******************************************************************************
 *                 TSparseGridSpinColorDiagonal implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSparseGridSpinColorDiagonal<FImpl>::TSparseGridSpinColorDiagonal(const std::string name)
: Module<SparseGridSpinColorDiagonalPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSparseGridSpinColorDiagonal<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSparseGridSpinColorDiagonal<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSparseGridSpinColorDiagonal<FImpl>::setup(void)
{
    envCreateDerived(DilutedNoise<FImpl>, 
                     SparseGridSpinColorDiagonalNoise<FImpl>,
                     getName(), 1, envGetGrid(FermionField), par().nsrc, par().nsparse);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSparseGridSpinColorDiagonal<FImpl>::execute(void)
{
    auto &noise = envGet(DilutedNoise<FImpl>, getName());
    LOG(Message) << "Generating sparsened, spin-color diagonal noise with nSparse = "
                    << par().nsparse << std::endl;
    noise.generateNoise(rng4d());
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_SparseGridSpinColorDiagonal_hpp_
