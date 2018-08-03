#ifndef Hadrons_MContraction_MesonFieldGamma_hpp_
#define Hadrons_MContraction_MesonFieldGamma_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>
#include <Grid/Hadrons/AllToAllVectors.hpp>
#include <Grid/Hadrons/AllToAllReduction.hpp>
#include <Grid/Grid_Eigen_Dense.h>
#include <fstream>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         MesonFieldGamma                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class MesonFieldPar : Serializable
{
  public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFieldPar,
                                    int, nBlock,
                                    std::string, A2A1,
                                    std::string, A2A2,
                                    std::string, gammas,
                                    std::string, output);
};

template <typename FImpl>
class TMesonFieldGamma : public Module<MesonFieldPar>
{
  public:
    FERM_TYPE_ALIASES(FImpl, );
    SOLVER_TYPE_ALIASES(FImpl, );

    typedef A2AModesSchurDiagTwo<typename FImpl::FermionField, FMat, Solver> A2ABase;

    class Result : Serializable
    {
      public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma,
                                        std::vector<std::vector<std::vector<ComplexD>>>, MesonField);
    };

  public:
    // constructor
    TMesonFieldGamma(const std::string name);
    // destructor
    virtual ~TMesonFieldGamma(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void parseGammaString(std::vector<Gamma::Algebra> &gammaList);
    virtual void vectorOfWs(std::vector<FermionField> &w, int i, int nBlock, FermionField &wTmp5D, std::vector<FermionField> &wVec);
    virtual void vectorOfVs(std::vector<FermionField> &v, int j, int nBlock, FermionField &vTmp5D, std::vector<FermionField> &vVec);
    virtual void gammaMult(std::vector<FermionField> &v, Gamma gamma);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER(MesonFieldGamma, ARG(TMesonFieldGamma<FIMPL>), MContraction);
MODULE_REGISTER(ZMesonFieldGamma, ARG(TMesonFieldGamma<ZFIMPL>), MContraction);

/******************************************************************************
*                  TMesonFieldGamma implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMesonFieldGamma<FImpl>::TMesonFieldGamma(const std::string name)
    : Module<MesonFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMesonFieldGamma<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().A2A1 + "_class", par().A2A2 + "_class"};
    in.push_back(par().A2A1 + "_w_high_4d");
    in.push_back(par().A2A2 + "_v_high_4d");
    return in;
}

template <typename FImpl>
std::vector<std::string> TMesonFieldGamma<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

template <typename FImpl>
void TMesonFieldGamma<FImpl>::parseGammaString(std::vector<Gamma::Algebra> &gammaList)
{
    gammaList.clear();
    // Determine gamma matrices to insert at source/sink.
    if (par().gammas.compare("all") == 0)
    {
        // Do all contractions.
        for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
        {
            gammaList.push_back(((Gamma::Algebra)i));
        }
    }
    else
    {
        // Parse individual contractions from input string.
        gammaList = strToVec<Gamma::Algebra>(par().gammas);
    }
}

template <typename FImpl>
void TMesonFieldGamma<FImpl>::vectorOfWs(std::vector<FermionField> &w, int i, int nBlock, FermionField &wTmp5D, std::vector<FermionField> &wVec)
{
    for (unsigned int ni = 0; ni < nBlock; ni++)
    {
        wVec[ni] = w[i + ni];
    }
}

template <typename FImpl>
void TMesonFieldGamma<FImpl>::vectorOfVs(std::vector<FermionField> &v, int j, int nBlock, FermionField &vTmp5D, std::vector<FermionField> &vVec)
{
    for (unsigned int nj = 0; nj < nBlock; nj++)
    {
        vVec[nj] = v[j+nj];
    }
}

template <typename FImpl>
void TMesonFieldGamma<FImpl>::gammaMult(std::vector<FermionField> &v, Gamma gamma)
{
    int nBlock = v.size();
    for (unsigned int nj = 0; nj < nBlock; nj++)
    {
        v[nj] = gamma * v[nj];
    }
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMesonFieldGamma<FImpl>::setup(void)
{
    int nt = env().getDim(Tp);
    auto &A2A1 = envGet(A2ABase, par().A2A1 + "_class");
    auto &A2A2 = envGet(A2ABase, par().A2A2 + "_class");
    int N1 = A2A1.nLow + A2A1.nHigh;
    int N2 = A2A2.nLow + A2A2.nHigh;
    int nBlock = par().nBlock;

    int Ls_ = env().getObjectLs(par().A2A1 + "_class");

    envTmpLat(FermionField, "vTmp5D", Ls_);
    envTmpLat(FermionField, "wTmp5D", Ls_);

    envTmp(std::vector<FermionField>, "w", 1, N1, FermionField(env().getGrid(1)));
    envTmp(std::vector<FermionField>, "v", 1, N2, FermionField(env().getGrid(1)));

    envTmp(Eigen::MatrixXcd, "MF", 1, Eigen::MatrixXcd::Zero(nt, N1 * N2));

    envTmp(std::vector<FermionField>, "wBlock", 1, nBlock, FermionField(env().getGrid(1)));
    envTmp(std::vector<FermionField>, "vBlock", 1, nBlock, FermionField(env().getGrid(1)));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMesonFieldGamma<FImpl>::execute(void)
{
    LOG(Message) << "Computing A2A meson field for gamma = " << par().gammas << ", taking w from " << par().A2A1 << " and v from " << par().A2A2 << std::endl;

    int nt = env().getDim(Tp);
    int nBlock = par().nBlock;

    std::vector<Result> result;
    std::vector<Gamma::Algebra> gammaResultList;
    std::vector<Gamma> gammaList;

    parseGammaString(gammaResultList);
    result.resize(gammaResultList.size());

    Gamma g5(Gamma::Algebra::Gamma5);
    gammaList.resize(gammaResultList.size(), g5);

    auto &A2A1 = envGet(A2ABase, par().A2A1 + "_class");
    auto &A2A2 = envGet(A2ABase, par().A2A2 + "_class");
    int N1 = A2A1.nLow + A2A1.nHigh;
    int N2 = A2A2.nLow + A2A2.nHigh;

    for (unsigned int i = 0; i < result.size(); ++i)
    {
        result[i].gamma = gammaResultList[i];
        result[i].MesonField.resize(N1, std::vector<std::vector<ComplexD>>(N2, std::vector<ComplexD>(nt)));

        Gamma gamma(gammaResultList[i]);
        gammaList[i] = gamma;
    }

    envGetTmp(FermionField, vTmp5D);
    envGetTmp(FermionField, wTmp5D);

    envGetTmp(std::vector<FermionField>, v);
    envGetTmp(std::vector<FermionField>, w);
    LOG(Message) << "Finding w vectors for N1 =  " << N1 << std::endl;
    for (int i = 0; i < N1; i++)
    {
        A2A1.wReturn(i, wTmp5D, w[i]);
    }
    LOG(Message) << "Finding v vectors for N2 =  " << N2 << std::endl;
    for (int i = 0; i < N2; i++)
    {
        A2A2.vReturn(i, vTmp5D, v[i]);
    }
    LOG(Message) << "Found v and w vectors" << std::endl;

    std::vector<std::vector<ComplexD>> MesonField_ij;
    LOG(Message) << "Before blocked MFs, nBlock = " << nBlock << std::endl;
    envGetTmp(std::vector<FermionField>, vBlock);
    envGetTmp(std::vector<FermionField>, wBlock);
    MesonField_ij.resize(nBlock * nBlock, std::vector<ComplexD>(nt));

    envGetTmp(Eigen::MatrixXcd, MF);

    LOG(Message) << "Before blocked MFs, nBlock = " << nBlock << std::endl;
    for (unsigned int i = 0; i < N1; i += nBlock)
    {
        vectorOfWs(w, i, nBlock, wTmp5D, wBlock);
        for (unsigned int j = 0; j < N2; j += nBlock)
        {
            vectorOfVs(v, j, nBlock, vTmp5D, vBlock);
            for (unsigned int k = 0; k < result.size(); k++)
            {
                gammaMult(vBlock, gammaList[k]);
                sliceInnerProductMesonField(MesonField_ij, wBlock, vBlock, Tp);
                for (unsigned int ni = 0; ni < nBlock; ni++)
                {
                    for (unsigned int nj = 0; nj < nBlock; nj++)
                    {
                        MF.col((i + ni) + (j + nj) * N1) = Eigen::VectorXcd::Map(&MesonField_ij[nj * nBlock + ni][0], MesonField_ij[nj * nBlock + ni].size());
                    }
                }
            }
        }
        if (i % 10 == 0)
        {
            LOG(Message) << "MF for i = " << i << " of " << N1 << std::endl;
        }
    }
    LOG(Message) << "Before Global sum, nBlock = " << nBlock << std::endl;
    vBlock[0]._grid->GlobalSumVector(MF.data(), MF.size());
    LOG(Message) << "After Global sum, nBlock = " << nBlock << std::endl;
    for (unsigned int i = 0; i < N1; i++)
    {
        for (unsigned int j = 0; j < N2; j++)
        {
            for (unsigned int k = 0; k < result.size(); k++)
            {
                for (unsigned int t = 0; t < nt; t++)
                {
                    result[k].MesonField[i][j][t] = MF.col(i + N1 * j)[t];
                }
            }
        }
    }
    saveResult(par().output, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MesonFieldGm_hpp_
