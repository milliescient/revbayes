#ifndef FossilOccurrenceDistribution_H
#define FossilOccurrenceDistribution_H

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "TypedDistribution.h"
#include "RbVector.h"
#include "Taxon.h"

namespace RevBayesCore {

    class AbstractFossilizedBirthDeathRangeProcess;

    /**
     * @brief The fossil record conditioned on a range process (Rev: dnFossilOccurrences).
     *
     * The ranges, psi, the occurrences and the reporting model are all held by the range process.
     */
    class FossilOccurrenceDistribution : public TypedDistribution< RbVector<Taxon> > {

    public:
        FossilOccurrenceDistribution(const DagNode *ranges, AbstractFossilizedBirthDeathRangeProcess::Reporting reporting, const std::vector<double> &bins = std::vector<double>());
        virtual ~FossilOccurrenceDistribution() {}

        FossilOccurrenceDistribution* clone(void) const override;
        double                      computeLnProbability(void) override;
        void                        redrawValue(void) override;
        void                        setValue(RbVector<Taxon> *v, bool force=false) override;

    protected:
        void                        swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;

        void                        resolveRanges(void);
        static AbstractFossilizedBirthDeathRangeProcess* rangesOf(const DagNode *n);

    private:
        const DagNode*                                   ranges_node;
        AbstractFossilizedBirthDeathRangeProcess*        ranges;
        AbstractFossilizedBirthDeathRangeProcess::Reporting reporting;
        std::vector<Taxon>                               taxa;
        std::vector<double>                              bins;                   //!< Interval boundaries for simulated occurrences, empty for exact ages.
    };
}

#endif
