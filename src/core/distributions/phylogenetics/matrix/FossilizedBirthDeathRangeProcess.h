#ifndef FossilizedBirthDeathRangeProcess_H
#define FossilizedBirthDeathRangeProcess_H

#include <set>
#include <utility>
#include "AbstractFossilizedBirthDeathRangeProcess.h"

#include "MatrixReal.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    /**
     * @brief Piecewise-constant fossilized birth-death range process on a matrix of ranges (Rev: dnFBDRP).
     *
     * Each row holds a taxon's death, youngest appearance, oldest appearance and birth.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class FossilizedBirthDeathRangeProcess : public TypedDistribution<MatrixReal>, public AbstractFossilizedBirthDeathRangeProcess {
        
    public:
        FossilizedBirthDeathRangeProcess (const DagNode *speciation,
										  const DagNode *extinction,
										  const DagNode *psi,
										  const TypedDagNode<double>* rho,
										  const TypedDagNode<RbVector<double> > *times,
										  const std::string &condition,
										  const std::vector<Taxon> &taxa,
										  bool complete_record,
										  bool coexistence,
										  bool require_parent,
                                          const TypedDagNode<double>* origin = NULL,
                                            double present = 0.0,
                                          const TypedDagNode<RbVector<double> >* psi_mixture_rates = NULL,
                                          const TypedDagNode<Simplex>* psi_mixture_weights = NULL,
                                          const TypedDagNode<RbVector<double> >* speciation_timeline = NULL,
                                          const TypedDagNode<RbVector<double> >* extinction_timeline = NULL,
                                          const TypedDagNode<RbVector<double> >* sampling_timeline = NULL);
        
        // public member functions
        void                                            setMcmcMode(bool tf) override;

        FossilizedBirthDeathRangeProcess*               clone(void) const override;                                         //!< Create an independent clone

        void                                            setValue(MatrixReal *v, bool force = false) override;
        std::string                                     getHiddenStateString(void) const override;                          //!< The matrix with the statuses and psi mixture components it does not carry.
        void                                            setHiddenStateFromString(const std::string &s) override;

    protected:
        void                                            updateRanges() override;

        // Parameter management functions
        double                                          computeLnProbability(void) override;                                //!< Compute the log-transformed probability of the current value.

        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;  //!< Swap a parameter

        std::vector<std::pair<size_t,double> >          stored_repairs;                                                     //!< Elements the repair wrote that no move stored, undone on restore.
        void                                            repairRanges(void) override;                                        //!< Hold tau_K equal to tau_1 below two occurrences.
        void                                            repairRanges(const std::set<size_t> &touched);                      //!< Follow whichever of the two columns the move wrote.
        void                                            keepSpecialization(const DagNode *toucher) override;
        void                                            restoreSpecialization(const DagNode *toucher) override;
        void                                            touchSpecialization(const DagNode *toucher, bool touchAll) override;

        bool                                            coexistence;                                        //!< Weight each birth by how many lineages were alive to produce it.
        bool                                            require_parent;                                     //!< Reject a birth with none alive. Off only where ranges were left out.

    private:
        
        // helper functions
        void                                            updateGamma(bool force = false);                                    //!< Number of species alive at time t.
        void                                            redrawValue(void) override;

        std::vector<size_t>                             gamma_i;            //!< The number of coexisting lineages at each taxon birth time
        std::vector<std::vector<bool> >                 gamma_links;        //!< Which taxa coexist at each birth time
        std::vector<bool>                               dirty_gamma;        //!< Indicates whether gamma needs updating
    };
}

#endif
