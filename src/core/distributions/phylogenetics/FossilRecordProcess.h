#ifndef FossilRecordProcess_H
#define FossilRecordProcess_H

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "TypedDistribution.h"
#include "RbVector.h"
#include "Taxon.h"

namespace RevBayesCore {

    class AbstractFossilizedBirthDeathRangeProcess;

    /**
     * @brief The fossil record conditioned on a range process (Rev: dnFossilRecord).
     *
     * The ranges, appearances, psi and the occurrences are all read from the range process,
     * which also holds the reporting model.
     */
    class FossilRecordProcess : public TypedDistribution< RbVector<Taxon> > {

    public:
        FossilRecordProcess(const DagNode *ranges, AbstractFossilizedBirthDeathRangeProcess::Reporting reporting);
        virtual ~FossilRecordProcess() {}

        FossilRecordProcess*        clone(void) const override;
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
    };
}

#endif
