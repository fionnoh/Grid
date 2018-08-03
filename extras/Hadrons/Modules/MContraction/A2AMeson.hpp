#ifndef Hadrons_MContraction_A2AMeson_hpp_
#define Hadrons_MContraction_A2AMeson_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/AllToAllVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2AMeson                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;

class A2AMesonPar : Serializable
{
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMesonPar,
                                    std::string, A2A1,
                                    std::string, A2A2,
                                    std::string, gammas,
                                    std::string, output);
};

template <typename FImpl>
class TA2AMeson : public Module<A2AMesonPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
    SOLVER_TYPE_ALIASES(FImpl, );

    typedef A2AModesSchurDiagTwo<typename FImpl::FermionField, FMat, Solver> A2ABase;

    class Result : Serializable
    {
      public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma_snk,
                                        Gamma::Algebra, gamma_src,
                                        std::vector<Complex>, corr);
    };

  public:
    // constructor
    TA2AMeson(const std::string name);
    // destructor
    virtual ~TA2AMeson(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void parseGammaString(std::vector<GammaPair> &gammaList);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER(A2AMeson, ARG(TA2AMeson<FIMPL>), MContraction);
MODULE_REGISTER(ZA2AMeson, ARG(TA2AMeson<ZFIMPL>), MContraction);

/******************************************************************************
*                  TA2AMeson implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AMeson<FImpl>::TA2AMeson(const std::string name)
    : Module<A2AMesonPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AMeson<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().A2A1 + "_class", par().A2A2 + "_class"};
    in.push_back(par().A2A1 + "_w_high_4d");
    in.push_back(par().A2A2 + "_v_high_4d");

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AMeson<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

template <typename FImpl>
void TA2AMeson<FImpl>::parseGammaString(std::vector<GammaPair> &gammaList)
{
    gammaList.clear();
    // Parse individual contractions from input string.
    gammaList = strToVec<GammaPair>(par().gammas);
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMeson<FImpl>::setup(void)
{
    int Ls_ = env().getObjectLs(par().A2A1 + "_class");
    auto &A2AClass1 = envGet(A2ABase, par().A2A1 + "_class");
    auto &A2AClass2 = envGet(A2ABase, par().A2A2 + "_class");

    int N1 = A2AClass1.nLow + A2AClass1.nHigh;
    int N2 = A2AClass2.nLow + A2AClass2.nHigh;
    int nt = env().getDim(Tp);

    envTmp(std::vector<FermionField>, "w1", 1, N1, FermionField(env().getGrid(1)));
    envTmp(std::vector<FermionField>, "v1", 1, N2, FermionField(env().getGrid(1)));
    envTmp(std::vector<FermionField>, "w2", 1, N1, FermionField(env().getGrid(1)));
    envTmp(std::vector<FermionField>, "v2", 1, N2, FermionField(env().getGrid(1)));
    envTmpLat(FermionField, "vTmp5D", Ls_);
    envTmpLat(FermionField, "wTmp5D", Ls_);

    envTmp(std::vector<ComplexD>, "MF1", 1, nt);
    envTmp(std::vector<ComplexD>, "MF2", 1, nt);
    envTmp(std::vector<ComplexD>, "tmp", 1, nt);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMeson<FImpl>::execute(void)
{
    LOG(Message) << "Computing A2A meson contractions" << std::endl;

    Result result;
    Gamma g5(Gamma::Algebra::Gamma5);
    std::vector<GammaPair> gammaList;
    int nt = env().getDim(Tp);

    parseGammaString(gammaList);

    result.gamma_snk = gammaList[0].first;
    result.gamma_src = gammaList[0].second;
    result.corr.resize(nt);

    auto &A2AClass1 = envGet(A2ABase, par().A2A1 + "_class");
    auto &A2AClass2 = envGet(A2ABase, par().A2A2 + "_class");
    int N = A2AClass1.nLow + A2AClass1.nHigh;
    LOG(Message) << "N for A2A cont: " << N << std::endl;

    envGetTmp(std::vector<ComplexD>, MF1);
    envGetTmp(std::vector<ComplexD>, MF2);
    envGetTmp(std::vector<ComplexD>, tmp);

    for (unsigned int t = 0; t < nt; ++t)
    {
        tmp[t] = TensorRemove(MF1[t] * MF2[t] * 0.0);
    }

    Gamma gSnk(gammaList[0].first);
    Gamma gSrc(gammaList[0].second);

    envGetTmp(std::vector<FermionField>, w1);
    envGetTmp(std::vector<FermionField>, v1);
    envGetTmp(std::vector<FermionField>, w2);
    envGetTmp(std::vector<FermionField>, v2);
    envGetTmp(FermionField, vTmp5D);
    envGetTmp(FermionField, wTmp5D);

    LOG(Message) << "Finding v and w vectors for N =  " << N << std::endl;
    for (int i = 0; i < N; i++)
    {
        A2AClass1.wReturn(i, wTmp5D, w1[i]);
        A2AClass1.vReturn(i, vTmp5D, v1[i]);
        A2AClass2.wReturn(i, wTmp5D, w2[i]);
        A2AClass2.vReturn(i, vTmp5D, v2[i]);
    }
    LOG(Message) << "Found v and w vectors for N =  " << N << std::endl;
    for (unsigned int i = 0; i < N; i++)
    {
        v1[i] = gSrc * v1[i];
        v2[i] = gSnk * v2[i];
    }
    int ty;
    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < N; j++)
        {
            mySliceInnerProductVector(MF1, w2[i], v1[j], Tp);
            mySliceInnerProductVector(MF2, w1[j], v2[i], Tp);
            for (unsigned int t = 0; t < nt; ++t)
            {
                for (unsigned int tx = 0; tx < nt; tx++)
                {
                    ty = (tx + t) % nt;
                    tmp[t] += TensorRemove((MF1[tx]) * (MF2[ty]));
                }
            }
        }
        if (i % 10 == 0)
        {
            LOG(Message) << "MF for i = " << i << " of " << N << std::endl;
        }
    }
    double NTinv = 1.0 / static_cast<double>(nt);
    for (unsigned int t = 0; t < nt; ++t)
    {
        result.corr[t] = NTinv * tmp[t];
    }

    saveResult(par().output, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AMeson_hpp_
