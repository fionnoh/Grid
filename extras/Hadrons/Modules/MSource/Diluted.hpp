#ifndef Hadrons_MSource_Diluted_hpp_
#define Hadrons_MSource_Diluted_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Diluted                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class DilutedPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DilutedPar,
                                    unsigned int, tStep,
                                    unsigned int, nHits);
};

template <typename FImpl>
class TDiluted: public Module<DilutedPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
  public:
    // constructor
    TDiluted(const std::string name);
    // destructor
    virtual ~TDiluted(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
  private:
    bool hasT_{false};
    std::string tName_;
};

MODULE_REGISTER_TMP(Diluted, TDiluted<FIMPL>, MSource);
MODULE_REGISTER_TMP(ScalarDiluted, TDiluted<ScalarImplCR>, MSource);

/******************************************************************************
 *                       TDiluted template implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDiluted<FImpl>::TDiluted(const std::string name)
: Module<DilutedPar>(name),
tName_(name + "_t")
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDiluted<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template <typename FImpl>
std::vector<std::string> TDiluted<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDiluted<FImpl>::setup(void)
{
    int nHits = par().nHits;
    int tStep = par().tStep;
    int Nt = env().getDim(Tp);
    if (!(Nt % tStep == 0))
    {
        HADRONS_ERROR(Argument, "Nt is not divisible by tStep");
    }
    int nH = nHits * (Nt / tStep);
    envCreate(std::vector<PropagatorField>, getName(), 1, nH, PropagatorField(env().getGrid()));
    envCacheLat(Lattice<iScalar<vInteger>>, tName_);
    envTmpLat(LatticeComplex, "eta");
    envTmpLat(LatticeComplex, "etaT");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDiluted<FImpl>::execute(void)
{
    auto &src = envGet(std::vector<PropagatorField>, getName());
    auto &t = envGet(Lattice<iScalar<vInteger>>, tName_);
    Complex shift(1., 1.);

    if (!hasT_)
    {
        LatticeCoordinate(t, Tp);
        hasT_ = true;
    }
    envGetTmp(LatticeComplex, eta);
    envGetTmp(LatticeComplex, etaT);

    int nHits = par().nHits;
    int tStep = par().tStep;
    int Nt = env().getDim(Tp);
    int nSrc;
    for (unsigned int i = 0; i < nHits; i++)
    {
        bernoulli(*env().get4dRng(), eta);
        eta = (2. * eta - shift) * (1. / ::sqrt(2.));
        for (unsigned int T = 0; T < Nt; T += tStep)
        {
            nSrc = (T/tStep) + i*(Nt/tStep);
            etaT = 1.;
            etaT = where((t >= T) and (t <= (T + tStep - 1)), etaT, 0. * etaT);
            src[nSrc] = 1.;
            src[nSrc] = src[nSrc] * etaT;
            src[nSrc] = src[nSrc] * eta;
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Diluted_hpp_
