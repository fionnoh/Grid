#ifndef Hadrons_MUtilities_PropMult_hpp_
#define Hadrons_MUtilities_PropMult_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         PropMult                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class PropMultPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PropMultPar,
                                    std::string, q,
                                    Gamma::Algebra, gamma);
};

template <typename FImpl>
class TPropMult: public Module<PropMultPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // constructor
    TPropMult(const std::string name);
    // destructor
    virtual ~TPropMult(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
};

MODULE_REGISTER_TMP(PropMult, TPropMult<FIMPL>, MUtilities);
MODULE_REGISTER_TMP(ZPropMult, TPropMult<ZFIMPL>, MUtilities);

/******************************************************************************
 *                 TPropMult implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPropMult<FImpl>::TPropMult(const std::string name)
: Module<PropMultPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPropMult<FImpl>::getInput(void)
{
    std::vector<std::string> in={par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TPropMult<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropMult<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().q);
    if (Ls_ > 1)
    {
        envCreateLat(PropagatorField, getName(), Ls_);
    }
    else
    {
        envCreateLat(PropagatorField, getName());
    }
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropMult<FImpl>::execute(void)
{
    auto  &prop = envGet(PropagatorField, getName());
    auto  &q  = envGet(PropagatorField, par().q);
    Gamma gamma(par().gamma);

    prop = q*gamma;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_PropMult_hpp_
