#ifndef Hadrons_MSolver_A2AVectors_hpp_
#define Hadrons_MSolver_A2AVectors_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/Solver.hpp>
#include <Grid/Hadrons/EigenPack.hpp>
#include <Grid/Hadrons/AllToAllVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2AVectors                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsPar,
                                  bool, return5D,
                                  int, nLow,
                                  std::string, sources,
                                  std::string, action,
                                  std::string, eigenPack,
                                  std::string, solver);
};

template <typename FImpl, int nBasis>
class TA2AVectors : public Module<A2AVectorsPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);

    typedef FermionEigenPack<FImpl> EPack;
    typedef CoarseFermionEigenPack<FImpl, nBasis> CoarseEPack;

    typedef A2AModesSchurDiagTwo<typename FImpl::FermionField, FMat, Solver> A2ABase;

  public:
    // constructor
    TA2AVectors(const std::string name);
    // destructor
    virtual ~TA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);

  private:
    unsigned int Ls_;
    std::string className_;
};

MODULE_REGISTER_TMP(A2AVectors, ARG(TA2AVectors<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectors, ARG(TA2AVectors<ZFIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);

/******************************************************************************
 *                 TA2AVectors implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
TA2AVectors<FImpl, nBasis>::TA2AVectors(const std::string name)
: Module<A2AVectorsPar>(name)
, className_ (name + "_class")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
std::vector<std::string> TA2AVectors<FImpl, nBasis>::getInput(void)
{
    int nLow = par().nLow;
    std::string subString = "";
    if (nLow > 0) subString = "_subtract";
    std::vector<std::string> in = {par().solver + subString, par().sources};

    return in;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TA2AVectors<FImpl, nBasis>::getReference(void)
{
    std::vector<std::string> ref = {par().action};

    if (!par().eigenPack.empty())
    {
        ref.push_back(par().eigenPack);
    }

    return ref;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TA2AVectors<FImpl, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {getName(), className_,
                                    getName() + "_w_high_5d", getName() + "_v_high_5d",
                                    getName() + "_w_high_4d", getName() + "_v_high_4d"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TA2AVectors<FImpl, nBasis>::setup(void)
{
    int nLow = par().nLow;
    auto &propSrcs = envGet(std::vector<PropagatorField>, par().sources);
    int nSrc = propSrcs.size();
    int Nc = FImpl::Dimension;

    // Number of high modes
    int nHigh = nSrc * Nc * Ns;
    bool return5D = par().return5D;
    int Ls;

    std::string subString = "";
    if (nLow > 0) subString = "_subtract";
    auto &solver = envGet(Solver, par().solver + subString);
    Ls = env().getObjectLs(par().solver + subString);

    auto &action = envGet(FMat, par().action);

    envTmpLat(FermionField, "fermSrc", Ls);
    envTmpLat(FermionField, "unphysFermSrc", Ls);
    envTmpLat(FermionField, "tmp");

    std::vector<FermionField> *evec;
    const std::vector<RealD> *eval;

    if (nLow > 0)
    {
        // Low modes
        auto &epack = envGet(EPack, par().eigenPack);

        LOG(Message) << "Creating a2a vectors " << getName() <<
                     " using eigenpack '" << par().eigenPack << "' ("
                     << epack.evec.size() << " modes)" <<
                     " and " << nHigh << " high modes." << std::endl;
        evec = &epack.evec;
        eval = &epack.eval;
    }
    else
    {
        LOG(Message) << "Creating a2a vectors " << getName() <<
                     " using " << nHigh << " high modes only." << std::endl;
    }

    int size_5d = 1;
    if (return5D) size_5d = nHigh;
    envCreate(std::vector<FermionField>, getName() + "_w_high_5d", Ls, size_5d, FermionField(env().getGrid(Ls)));
    envCreate(std::vector<FermionField>, getName() + "_v_high_5d", Ls, size_5d, FermionField(env().getGrid(Ls)));
    envCreate(std::vector<FermionField>, getName() + "_w_high_4d", 1, nHigh, FermionField(env().getGrid(1)));
    envCreate(std::vector<FermionField>, getName() + "_v_high_4d", 1, nHigh, FermionField(env().getGrid(1)));

    auto &wHigh5D = envGet(std::vector<FermionField>, getName() + "_w_high_5d");
    auto &vHigh5D = envGet(std::vector<FermionField>, getName() + "_v_high_5d");
    auto &wHigh4D = envGet(std::vector<FermionField>, getName() + "_w_high_4d");
    auto &vHigh4D = envGet(std::vector<FermionField>, getName() + "_v_high_4d");

    envCreate(A2ABase, className_, Ls,
                     evec, eval,
                     action,
                     solver,
                     wHigh5D, vHigh5D,
                     wHigh4D, vHigh4D,
                     nLow, nHigh,
                     return5D);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TA2AVectors<FImpl, nBasis>::execute(void)
{
    auto &action = envGet(FMat, par().action);

    int Nc = FImpl::Dimension;
    int Ls_;
    int nLow = par().nLow;

    std::string subString = "";
    if (nLow > 0) subString = "_subtract";
    Ls_ = env().getObjectLs(par().solver + subString);

    auto &A2AClass = envGet(A2ABase, className_);

    // High modes
    std::string sources = par().sources;
    auto &propSrcs = envGet(std::vector<PropagatorField>, sources);
    int nSrc = propSrcs.size();

    envGetTmp(FermionField, fermSrc);
    envGetTmp(FermionField, unphysFermSrc);
    envGetTmp(FermionField, tmp);

    int nCount = 0;
    for (unsigned int s = 0; s < Ns; ++s)
        for (unsigned int c = 0; c < Nc; ++c)
            for (unsigned int h = 0; h < nSrc; ++h)
            {
                LOG(Message) << "A2A src for s = " << s << " , c = " << c << ", h = " << h << std::endl;
                // source conversion for 4D sources
                if (!env().isObject5d(sources))
                {
                    if (Ls_ == 1)
                    {
                        PropToFerm<FImpl>(fermSrc, propSrcs[h], s, c);
                        tmp = fermSrc;
                    }
                    else
                    {
                        PropToFerm<FImpl>(tmp, propSrcs[h], s, c);
                        LOG(Message) << "ferm_tmp = " << tmp << std::endl;
                        action.ImportPhysicalFermionSource(tmp, fermSrc);
                        action.ImportUnphysicalFermion(tmp, unphysFermSrc);
                    }
                }
                // source conversion for 5D sources
                else
                {
                    if (Ls_ != env().getObjectLs(sources))
                    {
                        HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
                    }
                    else
                    {
                        PropToFerm<FImpl>(fermSrc, propSrcs[h], s, c);
                        action.ExportPhysicalFermionSolution(fermSrc, tmp);
                        unphysFermSrc = fermSrc;
                    }
                }
                LOG(Message) << "A2AClass.highModes Ncount = " << nCount << std::endl;
                A2AClass.highModes(fermSrc, unphysFermSrc, tmp, nCount);
                nCount++;
            }
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectors_hpp_
