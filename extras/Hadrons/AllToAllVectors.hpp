#ifndef A2A_Vectors_hpp_
#define A2A_Vectors_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Environment.hpp>
#include <Grid/Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

////////////////////////////////
// A2A Modes
////////////////////////////////

template <class Field, class Matrix, class Solver>
class A2AModesSchurDiagTwo
{
  private:
    const std::vector<Field> *evec;
    const std::vector<RealD> *eval;
    Matrix &action;
    Solver &solver;
    std::vector<Field> wHigh5D, vHigh5D, wHigh4D, vHigh4D;
    const bool return5D;
  public:
    const int nLow, nHigh;

  public:
    A2AModesSchurDiagTwo(const std::vector<Field> *_evec, const std::vector<RealD> *_eval,
                         Matrix &_action,
                         Solver &_solver,
                         std::vector<Field> _wHigh5D, std::vector<Field> _vHigh5D,
                         std::vector<Field> _wHigh4D, std::vector<Field> _vHigh4D,
                         const int _nLow, const int _nHigh,
                         const bool _return5D)
                        : evec(_evec), eval(_eval),
                        action(_action),
                        solver(_solver),
                        wHigh5D(_wHigh5D), vHigh5D(_vHigh5D),
                        wHigh4D(_wHigh4D), vHigh4D(_vHigh4D),
                        nLow(_nLow), nHigh(_nHigh),
                        return5D(_return5D){};

    void highModes(Field &source5D, Field &wSource5D, Field &source4D, int i)
    {
        int i5D;
        LOG(Message) << "A2A high modes for i = " << i << std::endl;
        i5D = 0;
        if (return5D) i5D = i;
        this->vHighMode(action, solver, source5D, vHigh5D[i5D], vHigh4D[i]);
        this->wHighMode(wSource5D, source4D, wHigh5D[i5D], wHigh4D[i]);
    }

    void vReturn(int i, Field &vOut5D, Field &vOut4D)
    {
        if (i < nLow)
        {
            this->vLowMode(action, evec->at(i), eval->at(i), vOut5D, vOut4D);
        }
        else
        {
            vOut4D = vHigh4D[i - nLow];
            if (!(return5D)) i = nLow;
            vOut5D = vHigh5D[i - nLow];
        }
    }
    void wReturn(int i, Field &wOut5D, Field &wOut4D)
    {
        if (i < nLow)
        {
            this->wLowMode(action, evec->at(i), eval->at(i), wOut5D, wOut4D);
        }
        else
        {
            wOut4D = wHigh4D[i - nLow];
            if (!(return5D)) i = nLow;
            wOut5D = wHigh5D[i - nLow];
        }
    }

    void vLowMode(Matrix &action, const Field &evec, const RealD &eval, Field &vOut5D, Field &vOut4D)
    {
        GridBase *grid = action.RedBlackGrid();
        Field srcO(grid);
        Field solE(grid);
        Field solO(grid);
        Field tmp(grid);

        srcO = evec;
        srcO.checkerboard = Odd;
        pickCheckerboard(Even, solE, vOut5D);
        pickCheckerboard(Odd, solO, vOut5D);

        /////////////////////////////////////////////////////
        // v_ie = -(1/eval_i) * MeeInv Meo MooInv evec_i
        /////////////////////////////////////////////////////
        action.MooeeInv(srcO, tmp);
        assert(tmp.checkerboard == Odd);
        action.Meooe(tmp, solE);
        assert(solE.checkerboard == Even);
        action.MooeeInv(solE, tmp);
        assert(tmp.checkerboard == Even);
        solE = (-1.0 / eval) * tmp;
        assert(solE.checkerboard == Even);

        /////////////////////////////////////////////////////
        // v_io = (1/eval_i) * MooInv evec_i
        /////////////////////////////////////////////////////
        action.MooeeInv(srcO, tmp);
        assert(tmp.checkerboard == Odd);
        solO = (1.0 / eval) * tmp;
        assert(solO.checkerboard == Odd);

        setCheckerboard(vOut5D, solE);
        assert(solE.checkerboard == Even);
        setCheckerboard(vOut5D, solO);
        assert(solO.checkerboard == Odd);

        action.ExportPhysicalFermionSolution(vOut5D, vOut4D);
    }

    void wLowMode(Matrix &action, const Field &evec, const RealD &eval, Field &wOut5D, Field &wOut4D)
    {
        GridBase *grid = action.RedBlackGrid();
        SchurDiagTwoOperator<Matrix, Field> _HermOpEO(action);

        Field srcO(grid);
        Field solE(grid);
        Field solO(grid);
        Field tmp(grid);

        GridBase *fgrid = action.Grid();
        Field tmp_wout(fgrid);

        srcO = evec;
        srcO.checkerboard = Odd;
        pickCheckerboard(Even, solE, tmp_wout);
        pickCheckerboard(Odd, solO, tmp_wout);

        /////////////////////////////////////////////////////
        // w_ie = - MeeInvDag MoeDag Doo evec_i
        /////////////////////////////////////////////////////
        _HermOpEO.Mpc(srcO, tmp);
        assert(tmp.checkerboard == Odd);
        action.MeooeDag(tmp, solE);
        assert(solE.checkerboard == Even);
        action.MooeeInvDag(solE, tmp);
        assert(tmp.checkerboard == Even);
        solE = (-1.0) * tmp;

        /////////////////////////////////////////////////////
        // w_io = Doo evec_i
        /////////////////////////////////////////////////////
        _HermOpEO.Mpc(srcO, solO);
        assert(solO.checkerboard == Odd);

        setCheckerboard(tmp_wout, solE);
        assert(solE.checkerboard == Even);
        setCheckerboard(tmp_wout, solO);
        assert(solO.checkerboard == Odd);

        action.DminusDag(tmp_wout, wOut5D);

        action.ExportPhysicalFermionSource(wOut5D, wOut4D);
    }

    void vHighMode(Matrix &action, Solver &solver, const Field &source, Field &vOut5D, Field &vOut4D)
    {
        GridBase *fgrid = action.Grid();
        solver(vOut5D, source); // Note: solver is solver(out, in)
        action.ExportPhysicalFermionSolution(vOut5D, vOut4D);
    }

    void wHighMode(const Field &wSource5D, const Field &source4D, Field &wOut5D, Field &wOut4D)
    {
        wOut5D = wSource5D;
        wOut4D = source4D;
    }
};

// TODO: A2A for coarse eigenvectors

// template <class FineField, class CoarseField, class Matrix, class Solver>
// class A2ALMSchurDiagTwoCoarse : public A2AModesSchurDiagTwo<FineField, Matrix, Solver>
// {
//   private:
//     const std::vector<FineField> &subspace;
//     const std::vector<CoarseField> &evec_coarse;
//     const std::vector<RealD> &eval_coarse;
//     Matrix &action;

//   public:
//     A2ALMSchurDiagTwoCoarse(const std::vector<FineField> &_subspace, const std::vector<CoarseField> &_evec_coarse, const std::vector<RealD> &_eval_coarse, Matrix &_action)
//         : subspace(_subspace), evec_coarse(_evec_coarse), eval_coarse(_eval_coarse), action(_action){};

//     void operator()(int i, FineField &vout, FineField &wout)
//     {
//         FineField prom_evec(subspace[0]._grid);
//         blockPromote(evec_coarse[i], prom_evec, subspace);
//         this->vLowMode(action, prom_evec, eval_coarse[i], vout);
//         this->wLowMode(action, prom_evec, eval_coarse[i], wout);
//     }
// };

END_HADRONS_NAMESPACE

#endif // A2A_Vectors_hpp_