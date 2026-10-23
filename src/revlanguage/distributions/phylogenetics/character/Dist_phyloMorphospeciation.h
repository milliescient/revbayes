#ifndef Dist_phyloMorphospeciation_H
#define Dist_phyloMorphospeciation_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "PhyloMorphospeciationProcess.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {

    /**
     * RevLanguage wrapper of dnPhyloMorphospeciation, morphological characters whose changes are
     * the speciation events of a dnFBDSP tree.
     *
     *     morph ~ dnPhyloMorphospeciation( tree=bd )
     *     morph.clamp( matrix )
     */
    class Dist_phyloMorphospeciation : public TypedDistribution< AbstractHomologousDiscreteCharacterData > {

    public:
        Dist_phyloMorphospeciation( void );

        Dist_phyloMorphospeciation*                     clone(void) const;
        static const std::string&                       getClassType(void);
        static const TypeSpec&                          getClassTypeSpec(void);
        std::string                                     getDistributionFunctionName(void) const;
        MethodTable                                     getDistributionMethods(void) const;
        const TypeSpec&                                 getTypeSpec(void) const;
        const MemberRules&                              getParameterRules(void) const;

        RevBayesCore::PhyloMorphospeciationProcess*     createDistribution(void) const;

    protected:
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);

        RevPtr<const RevVariable>                       tree;
        RevPtr<const RevVariable>                       num_states;
        RevPtr<const RevVariable>                       epsilon;
        RevPtr<const RevVariable>                       epsilon_a;
    };

}

#endif
