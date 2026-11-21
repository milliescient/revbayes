#ifndef PhyloCTMCSiteHomogeneousConditional_H
#define PhyloCTMCSiteHomogeneousConditional_H

#include <bitset>
#include <cmath>
#include "PhyloCTMCSiteHomogeneous.h"
#include "DistributionBinomial.h"
#include "DistributionNegativeBinomial.h"
#include "RlUserInterface.h"
#include "RateMatrix_JC.h"
#include "RbMathCombinatorialFunctions.h"

// A, B and C per singleton count, packed three to a slot
#define RB_FAST_A(u,j) ((u)[(j)*3    ])
#define RB_FAST_B(u,j) ((u)[(j)*3 + 1])
#define RB_FAST_C(u,j) ((u)[(j)*3 + 2])


namespace RevBayesCore {

    struct AscertainmentBias {
        
        enum Coding { ALL                 = 0x00,
                      VARIABLE            = 0x03,
                      INFORMATIVE         = 0x0F,
                      ORDERED             = 0x10    // uninformative patterns are those of ordered (Wagner) parsimony
                    };
    };

    template<class charType>
    class PhyloCTMCSiteHomogeneousConditional : public PhyloCTMCSiteHomogeneous<charType> {

    public:
        PhyloCTMCSiteHomogeneousConditional(const TypedDagNode< Tree > *t, size_t nChars, bool c, size_t nSites, bool amb, AscertainmentBias::Coding cod = AscertainmentBias::ALL, bool internal = false, bool gapmatch = true);
        PhyloCTMCSiteHomogeneousConditional(const PhyloCTMCSiteHomogeneousConditional &n);
        virtual                                            ~PhyloCTMCSiteHomogeneousConditional(void);                                                                   //!< Virtual destructor

        // public member functions
        PhyloCTMCSiteHomogeneousConditional*                clone(void) const;                                                                        //!< Create an independent clone
        void                                                setValue(AbstractHomologousDiscreteCharacterData *v, bool f=false);
        virtual void                                        redrawValue(void);

    protected:

        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r);
        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r, size_t m);
        virtual void                                        computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                        computeTipLikelihood(const TopologyNode &node, size_t nIdx);

        virtual void                                        computeRootCorrection(size_t root, size_t l, size_t r);
        virtual void                                        computeRootCorrection(size_t root, size_t l, size_t r, size_t m);
        virtual void                                        computeInternalNodeCorrection(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                        computeInternalNodeCorrection(const TopologyNode &n, size_t nIdx, size_t l, size_t r, size_t m);
        virtual void                                        computeTipCorrection(const TopologyNode &node, size_t nIdx);

        virtual void                                        resizeLikelihoodVectors(void);

        vector<double>::const_iterator                      getCorrectionLikelihoodsForNode(int node_index) const;
        vector<double>::iterator                            getMutableCorrectionLikelihoodsForNode(int node_index);

        bool                                                warned;

        int                                                 coding;
        bool                                                ordered;
        size_t                                              N;

        size_t                                              numCorrectionMasks;
        size_t                                              activeCorrectionOffset;
        size_t                                              correctionNodeOffset;
        size_t                                              correctionMixtureOffset;
        size_t                                              correctionMaskOffset;
        size_t                                              correctionOffset;
        size_t                                              numCorrectionPatterns;

        std::vector<double>                                 correctionLikelihoods;

        std::vector<std::vector<bool> >                     correctionMaskMatrix;
        std::vector<size_t>                                 correctionMaskCounts;
        std::vector<size_t>                                 maskObservationCounts;
        // one restricted entry per pattern needing it, count zero so it is never a divisor
        size_t                                              numRealMasks;
        std::vector<size_t>                                 restrictedPatternOfEntry;
        std::vector<size_t>                                 restrictedEntryOfPattern;
        std::vector<double>                                 perMaskMixtureCorrections;

        virtual double                                      sumRootLikelihood( void );
        virtual void                                        computeRootLikelihoods( std::vector<double> &rv ) const;
        double                                              compatibleExcludedProb( size_t mask ) const;
        double                                              excludedProbForMixture( std::vector<double>::const_iterator p_node, size_t mask, size_t mixture ) const;
        void                                                buildRestrictedEntries( void );
        virtual bool                                        isSitePatternCompatible( std::map<size_t, size_t> );
        virtual bool                                        isSitePatternCompatible( std::map<RbBitSet, size_t> );
        std::vector<size_t>                                 getIncludedSiteIndices( void );
        void                                                updateCorrections( const TopologyNode& node, size_t nodeIndex );

        // exchangeable states: the probability depends on the size of S alone
        bool                                                useFastInformative(void) const;
        mutable int                                         fastInformativeCache = -1;     //!< -1 unknown, else 0 or 1
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        virtual double                                      computeLnProbability(void);
        void                                                fastMerge(std::vector<double>::const_iterator l, std::vector<double>::const_iterator r, double *d) const;
        void                                                computeTipCorrectionFast(const TopologyNode &node, size_t nodeIndex);
        void                                                computeInternalNodeCorrectionFast(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right);
        void                                                computeInternalNodeCorrectionFast(const TopologyNode &node, size_t nodeIndex, size_t left, size_t right, size_t middle);
        void                                                computeRootCorrectionFast(size_t root, size_t left, size_t right);
        void                                                computeRootCorrectionFast(size_t root, size_t left, size_t right, size_t middle);

    private:
        virtual void                                        simulateConditional( const TopologyNode &node, std::vector<charType> &taxa, size_t rateIndex, std::map<size_t, size_t>& charCounts);
    };

}

#include <climits>

template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::PhyloCTMCSiteHomogeneousConditional(const TypedDagNode<Tree> *t, size_t nChars, bool c, size_t nSites, bool amb, AscertainmentBias::Coding ty, bool internal, bool gapmatch) :
    PhyloCTMCSiteHomogeneous<charType>(  t, nChars, c, nSites, amb, internal, gapmatch ), warned(false), coding(ty & ~AscertainmentBias::ORDERED), ordered(ty & AscertainmentBias::ORDERED), N(nSites), numCorrectionMasks(0), numRealMasks(0)
{
    if (coding != AscertainmentBias::ALL)
    {
        numCorrectionMasks      = 1;
        correctionMaskMatrix    = std::vector<std::vector<bool> >(1, std::vector<bool>(this->num_nodes,0) );
        
        // number of correction patterns per character state
        if (coding & (AscertainmentBias::INFORMATIVE ^ AscertainmentBias::VARIABLE))
        {
            numCorrectionPatterns = std::pow(2.0f, float(this->num_chars - 1));
        }
        else
        {
            numCorrectionPatterns = 1;
        }

        correctionOffset        = this->num_chars*numCorrectionPatterns;
        correctionMaskOffset    = correctionOffset*this->num_chars;
        correctionMixtureOffset = numCorrectionMasks*correctionMaskOffset;
        correctionNodeOffset    = this->num_site_mixtures*correctionMixtureOffset;
        activeCorrectionOffset  = this->num_nodes*correctionNodeOffset;

        correctionLikelihoods = std::vector<double>(activeCorrectionOffset*2, 0.0);

        perMaskMixtureCorrections = std::vector<double>(this->num_site_mixtures*numCorrectionMasks, 0.0);
    }
}


template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::~PhyloCTMCSiteHomogeneousConditional( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!

}


template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>* RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::clone( void ) const {

    return new PhyloCTMCSiteHomogeneousConditional<charType>( *this );
}

template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::PhyloCTMCSiteHomogeneousConditional(const PhyloCTMCSiteHomogeneousConditional &n) :
    PhyloCTMCSiteHomogeneous< charType>( n ),
    warned(n.warned),
    coding(n.coding),
    ordered(n.ordered),
    N(n.N),
    numCorrectionMasks(n.numCorrectionMasks),
    activeCorrectionOffset(n.activeCorrectionOffset),
    correctionNodeOffset(n.correctionNodeOffset),
    correctionMixtureOffset(n.correctionMixtureOffset),
    correctionMaskOffset(n.correctionMaskOffset),
    correctionOffset(n.correctionOffset),
    numCorrectionPatterns(n.numCorrectionPatterns),
    correctionLikelihoods(n.correctionLikelihoods),
    correctionMaskMatrix(n.correctionMaskMatrix),
    correctionMaskCounts(n.correctionMaskCounts),
    maskObservationCounts(n.maskObservationCounts),
numRealMasks(n.numRealMasks),
restrictedPatternOfEntry(n.restrictedPatternOfEntry),
restrictedEntryOfPattern(n.restrictedEntryOfPattern),
    perMaskMixtureCorrections(n.perMaskMixtureCorrections)
{

}

template<class charType>
std::vector<size_t> RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::getIncludedSiteIndices( void )
{
    correctionMaskCounts.clear();
    maskObservationCounts.clear();
    correctionMaskMatrix.clear();

    // find the unique site patterns and compute their respective frequencies
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    size_t tips = this->tau->getValue().getNumberOfTips();

    // create a vector with the correct site indices
    // some of the sites may have been excluded
    std::vector<size_t> siteIndices;
    size_t siteIndex = 0;
    size_t incompatible = 0;

    std::map<std::string, size_t> maskIndices;

    if (coding != AscertainmentBias::ALL)
    {
        // if we are using a correction, add a correction mask with 0 gaps.
        // it is required when simulating data, but not necessarily used
        // in computing the likelihood (e.g. if all sites have at least one gap)

        std::string gapless(tips, ' ');
        std::vector<bool> mask(tips, false);

        maskIndices[gapless] = 0;

        correctionMaskCounts.push_back(0);
        maskObservationCounts.push_back(tips);
        correctionMaskMatrix.push_back(mask);
    }

    this->num_sites = this->value->getNumberOfIncludedCharacters();

    for (size_t i = 0; i < this->num_sites; ++i)
    {
        while ( this->value->isCharacterExcluded(siteIndex) )
        {
            siteIndex++;
            if ( siteIndex >= this->value->getNumberOfCharacters()  )
            {
                throw RbException( "The character matrix cannot set to this variable because it does not have enough included characters." );
            }
        }


        std::map<size_t,size_t> charCounts;
        std::map<RbBitSet,size_t> bitCounts;
        size_t numGap = 0;

        std::string mask = "";
        std::vector<bool> maskData;
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                AbstractDiscreteTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );
                DiscreteCharacterState &c = taxon.getCharacter(siteIndex);

                bool gap = c.isGapState();
                // if we treat unknown characters as gaps and this is an unknown character then we change it
                // because we might then have a pattern more
                if ( this->treatAmbiguousAsGaps && (c.isAmbiguous() || c.isMissingState()) )
                {
                    gap = true;
                }
                else if ( this->treatUnknownAsGap && (c.getNumberOfStates() == c.getNumberObservedStates() || c.isMissingState()) )
                {
                    gap = true;
                }

                if (gap)
                {
                    if (coding != AscertainmentBias::ALL)
                        mask += "-";
                    numGap++;
                }
                else
                {
                    if (coding != AscertainmentBias::ALL)
                        mask += " ";

                    bitCounts[c.getState()]++;
                }

                if (coding != AscertainmentBias::ALL)
                    maskData.push_back(gap);
            }
        }

        if ( !isSitePatternCompatible(bitCounts) )
        {
            incompatible++;
        }
        else
        {
            siteIndices.push_back(siteIndex);

            if (coding != AscertainmentBias::ALL)
            {
                // increase the count for this mask
                std::map<std::string, size_t>::iterator it = maskIndices.find(mask);
                if (it != maskIndices.end())
                {
                    correctionMaskCounts[it->second]++;
                }
                else
                {
                    maskIndices[mask] = correctionMaskCounts.size();

                    correctionMaskCounts.push_back(1);
                    maskObservationCounts.push_back(tips - numGap);
                    correctionMaskMatrix.push_back(maskData);
                }
            }
        }

        siteIndex++;
    }

    // test if there were additional sites that we did not use
    while ( siteIndex < this->value->getNumberOfCharacters() )
    {
        if ( !this->value->isCharacterExcluded(siteIndex)  )
        {
            throw RbException( "The character matrix cannot set to this variable because it has too many included characters." );
        }
        siteIndex++;
    }

    // warn if we have to exclude some incompatible characters.
    if (incompatible > 0)
    {
        if (!warned)
        {
            std::stringstream ss;
            ss << "WARNING: There are " << incompatible << " characters incompatible with the specified coding bias. These characters will be excluded.";
            RBOUT(ss.str());

            warned = true;
        }

        // resize our datset to account for the newly excluded characters
        this->num_sites = siteIndices.size();
        this->site_pattern.clear();
        this->site_pattern.resize(this->num_sites);
    }

    // readjust the number of correction sites to account for masked sites
    if (coding != AscertainmentBias::ALL)
    {
        numCorrectionMasks = correctionMaskCounts.size();
        numRealMasks       = numCorrectionMasks;
    }

    return siteIndices;
}


template<class charType>
inline vector<double>::iterator RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::getMutableCorrectionLikelihoodsForNode( int node_index )
{
    return correctionLikelihoods.begin() + this->activeLikelihood[node_index]*this->activeCorrectionOffset + node_index*correctionNodeOffset;
}

template<class charType>
inline vector<double>::const_iterator RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::getCorrectionLikelihoodsForNode( int node_index ) const
{
    return correctionLikelihoods.begin() + this->activeLikelihood[node_index]*this->activeCorrectionOffset + node_index*correctionNodeOffset;
}

template<class charType>
bool RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::isSitePatternCompatible( std::map<size_t, size_t> charCounts )
{

    size_t num_informative = 0;
    for (std::map<size_t, size_t>::iterator it = charCounts.begin(); it != charCounts.end(); it++)
    {
        if (it->second > 1)
        {
            num_informative++;
        }
    }

    bool informative = num_informative > 1;

    // ordered: informative if some threshold between adjacent states has two taxa on each side
    if (ordered)
    {
        informative = false;
        for (size_t j = 1; j < this->num_chars && informative == false; j++)
        {
            size_t below = 0;
            size_t above = 0;
            for (std::map<size_t, size_t>::iterator it = charCounts.begin(); it != charCounts.end(); it++)
            {
                (it->first < j ? below : above) += it->second;
            }
            informative = (below > 1 && above > 1);
        }
    }

    if ( (charCounts.size() == 1 && (coding == AscertainmentBias::VARIABLE)) || (!informative && (coding == AscertainmentBias::INFORMATIVE)) )
    {
        return false;
    }

    return true;
}


// Thanks to MH
template<class charType>
bool RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::isSitePatternCompatible( std::map<RbBitSet, size_t> charCounts )
{
    // A character could have been recorded if any resolution of its ambiguous cells satisfies the
    // coding, so each test below asks whether such a resolution exists rather than reading the
    // cells as though each distinct set were a state of its own.
    size_t k = this->num_chars;

    // a cell admitting every state carries nothing, so it can neither make a character vary nor
    // supply one of the taxa an informative pattern needs
    std::map<RbBitSet, size_t> scored;
    for (std::map<RbBitSet, size_t>::iterator it = charCounts.begin(); it != charCounts.end(); it++)
        if ( it->first.count() > 0 && it->first.count() < k )
            scored[it->first] = it->second;

    size_t total = 0;
    for (std::map<RbBitSet, size_t>::iterator it = scored.begin(); it != scored.end(); it++)
        total += it->second;

    std::vector<size_t> admit(k, 0);
    for (std::map<RbBitSet, size_t>::iterator it = scored.begin(); it != scored.end(); it++)
        for (size_t i = 0; i < k && i < it->first.size(); i++)
            if ( it->first.test(i) )
                admit[i] += it->second;

    if ( coding == AscertainmentBias::VARIABLE )
    {
        // every resolution is constant only when every cell is the same single state
        return !( total > 0 && scored.size() == 1 && scored.begin()->first.count() == 1 );
    }

    if ( total < 4 )
        return false;

    // Two cells must take one state and two others a second state. Writing A and B for the cells
    // admitting each, such an assignment exists exactly when |A| >= 2, |B| >= 2 and |A or B| >= 4.
    if ( ordered )
    {
        for (size_t j = 1; j < k; j++)
        {
            size_t below = 0, above = 0;
            for (std::map<RbBitSet, size_t>::iterator it = scored.begin(); it != scored.end(); it++)
            {
                bool any_below = false, any_above = false;
                for (size_t i = 0; i < k && i < it->first.size(); i++)
                    if ( it->first.test(i) )
                        ( i < j ? any_below : any_above ) = true;
                if (any_below) below += it->second;
                if (any_above) above += it->second;
            }
            if ( below >= 2 && above >= 2 )
                return true;
        }
        return false;
    }

    for (size_t x = 0; x < k; x++)
    {
        for (size_t y = x + 1; y < k; y++)
        {
            if ( admit[x] < 2 || admit[y] < 2 )
                continue;
            size_t either = 0;
            for (std::map<RbBitSet, size_t>::iterator it = scored.begin(); it != scored.end(); it++)
                if ( it->first.test(x) || it->first.test(y) )
                    either += it->second;
            if ( either >= 4 )
                return true;
        }
    }
    return false;
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeRootLikelihood( size_t root, size_t left, size_t right)
{

    PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihood(root, left, right);

    if (coding != AscertainmentBias::ALL)
    {
        computeRootCorrection(root, left, right);
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle)
{

    PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihood(root, left, right, middle);

    if (coding != AscertainmentBias::ALL)
    {
        computeRootCorrection(root, left, right, middle);
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{

    PhyloCTMCSiteHomogeneous<charType>::computeInternalNodeLikelihood(node, node_index, left, right);

    if (coding != AscertainmentBias::ALL)
    {
        computeInternalNodeCorrection(node, node_index, left, right);
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeTipLikelihood(const TopologyNode &node, size_t node_index)
{

    PhyloCTMCSiteHomogeneous<charType>::computeTipLikelihood(node, node_index);

    if (coding != AscertainmentBias::ALL)
    {
        computeTipCorrection(node, node_index);
    }
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeTipCorrection(const TopologyNode &node, size_t node_index)
{
    if ( useFastInformative() )
    {
        computeTipCorrectionFast(node, node_index);
        return;
    }

    // ensure we actually have a matrix to use
    this->updateTransitionProbabilityMatrix(node_index);

    std::vector<double>::iterator p_node = getMutableCorrectionLikelihoodsForNode(node_index);
    
    size_t data_tip_index = this->taxon_name_2_tip_index_map[ node.getName() ];
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // const TransitionProbabilityMatrix&    pij = this->transition_prob_matrices[mixture];
        const TransitionProbabilityMatrix&    pij = this->pmatrices[node_index][mixture];

        // iterate over correction masks
        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            bool gap = correctionMaskMatrix[mask][data_tip_index];
            
            // a restricted entry admits only the states its pattern records
            RbBitSet admits;
            if ( mask >= numRealMasks && mask - numRealMasks < restrictedPatternOfEntry.size() )
                admits = this->ambiguous_char_matrix[data_tip_index][ restrictedPatternOfEntry[mask - numRealMasks] ];

            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < this->num_chars; a++)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;

                // iterate over constant/autapomorphic states
                for (size_t i = 0; i < this->num_chars; i++)
                {
                    if (i > 0 && numCorrectionPatterns == 1)
                        break;

                    // get bit pattern for this constant/autapomorphic state
                    size_t c = ( i == 0 ? 0 : std::pow(2.0f, float(i - 1)) );

                    std::vector<double>::iterator         uc = u  + c*this->num_chars;

                    // if we have a gap, then fill with ones for constant patterns
                    // or zero for other patterns
                    if (gap)
                    {
                        std::fill(uc, uc + this->num_chars, i == 0);
                    }
                    // if this is a constant pattern, fill with transition prob from ci to a
                    else if (c == 0)
                    {
                        // iterate over initial states
                        for (size_t ci = 0; ci < this->num_chars; ci++)
                        {
                            uc[ci] = ( admits.empty() || admits.test(a) ) ? pij[ci][a] : 0.0;
                        }
                    }
                    // if this is an autapomorphic state cj, fill with transition prob from ci to cj
                    else
                    {
                        // get autamorphic state
                        size_t cj = (i-1 >= a ? i : i-1);

                        // iterate over initial states
                        for (size_t ci = 0; ci < this->num_chars; ci++)
                        {
                            uc[ci] = ( admits.empty() || admits.test(cj) ) ? pij[ci][cj] : 0.0;
                        }
                    }
                }
            }
        }
    }
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeInternalNodeCorrection(const TopologyNode &node, size_t node_index, size_t left, size_t right, size_t middle)
{
    if ( useFastInformative() )
    {
        computeInternalNodeCorrectionFast(node, node_index, left, right, middle);
        return;
    }

    // ensure we actually have a matrix to use
    this->updateTransitionProbabilityMatrix(node_index);

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    std::vector<double>::const_iterator   p_left   = getCorrectionLikelihoodsForNode(left);
    std::vector<double>::const_iterator   p_right  = getCorrectionLikelihoodsForNode(right);
    std::vector<double>::const_iterator   p_middle = getCorrectionLikelihoodsForNode(middle);
    std::vector<double>::iterator         p_node   = getMutableCorrectionLikelihoodsForNode(node_index);
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // const TransitionProbabilityMatrix&    pij = this->transition_prob_matrices[mixture];
        const TransitionProbabilityMatrix&    pij = this->pmatrices[node_index][mixture];

        // iterate over correction masks
        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < this->num_chars; a++)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;
                std::vector<double>::const_iterator         u_l = p_left   + offset;
                std::vector<double>::const_iterator         u_r = p_right  + offset;
                std::vector<double>::const_iterator         u_m = p_middle + offset;

                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    std::vector<double>::iterator         uc = u  + c*this->num_chars;

                    std::fill(uc, uc + this->num_chars, 0.0);

                    // iterate over partitions of c
                    for (size_t p1 = 0; p1 <= c; p1++)
                    {
                        if ( (p1 | c) == c)
                        {
                            size_t p_tmp = p1 ^ c;

                            // iterate over partitions of p_tmp
                            for (size_t p2 = 0; p2 <= p_tmp; p2++)
                            {
                                if ( (p2 | p_tmp) == p_tmp)
                                {
                                    size_t p3 = p2 ^ p_tmp;

                                    std::vector<double>::const_iterator         lc = u_l  + p1*this->num_chars;
                                    std::vector<double>::const_iterator         rc = u_r  + p2*this->num_chars;
                                    std::vector<double>::const_iterator         mc = u_m  + p3*this->num_chars;

                                    // iterate over initial states
                                    for (size_t ci = 0; ci < this->num_chars; ci++)
                                    {
                                        // iterate over ending states
                                        for (size_t cj = 0; cj < this->num_chars; cj++)
                                        {
                                            uc[ci] += pij[ci][cj] * lc[cj] * rc[cj] * mc[cj];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeInternalNodeCorrection(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{
    if ( useFastInformative() )
    {
        computeInternalNodeCorrectionFast(node, node_index, left, right);
        return;
    }

    // ensure we actually have a matrix to use
    this->updateTransitionProbabilityMatrix(node_index);

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    std::vector<double>::const_iterator   p_left  = getCorrectionLikelihoodsForNode(left);
    std::vector<double>::const_iterator   p_right = getCorrectionLikelihoodsForNode(right);
    std::vector<double>::iterator         p_node  = getMutableCorrectionLikelihoodsForNode(node_index);
    
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // const TransitionProbabilityMatrix&    pij = this->transition_prob_matrices[mixture];
        const TransitionProbabilityMatrix&    pij = this->pmatrices[node_index][mixture];

        // iterate over correction masks
        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < this->num_chars; a++)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;
                std::vector<double>::const_iterator         u_l = p_left   + offset;
                std::vector<double>::const_iterator         u_r = p_right  + offset;

                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    std::vector<double>::iterator         uc = u  + c*this->num_chars;

                    std::fill(uc, uc + this->num_chars, 0.0);

                    // iterate over partitions of c
                    for (size_t p1 = 0; p1 <= c; p1++)
                    {
                        if ( (p1 | c) == c)
                        {
                            size_t p2 = p1 ^ c;

                            std::vector<double>::const_iterator         lc = u_l  + p1*this->num_chars;
                            std::vector<double>::const_iterator         rc = u_r  + p2*this->num_chars;

                            // iterate over initial states
                            for (size_t ci = 0; ci < this->num_chars; ci++)
                            {
                                // iterate over ending states
                                for (size_t cj = 0; cj < this->num_chars; cj++)
                                {
                                    uc[ci] += pij[ci][cj] * lc[cj] * rc[cj];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeRootCorrection( size_t root, size_t left, size_t right, size_t middle)
{
    if ( useFastInformative() )
    {
        computeRootCorrectionFast(root, left, right, middle);
        return;
    }

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    std::vector<double>::iterator         p_node   = getMutableCorrectionLikelihoodsForNode(root);
    std::vector<double>::const_iterator   p_left   = getCorrectionLikelihoodsForNode(left);
    std::vector<double>::const_iterator   p_right  = getCorrectionLikelihoodsForNode(right);
    std::vector<double>::const_iterator   p_middle = getCorrectionLikelihoodsForNode(middle);

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    this->getRootFrequencies(ff);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        const std::vector<double> &f = ff[mixture % ff.size()];

        // iterate over correction masks
        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < this->num_chars; a++)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;
                std::vector<double>::const_iterator         u_l = p_left   + offset;
                std::vector<double>::const_iterator         u_r = p_right  + offset;
                std::vector<double>::const_iterator         u_m = p_middle + offset;

                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    std::vector<double>::iterator         uc = u  + c*this->num_chars;

                    std::fill(uc, uc + this->num_chars, 0.0);

                    // iterate over partitions of c
                    for (size_t p1 = 0; p1 <= c; p1++)
                    {
                        if ( (p1 | c) == c)
                        {
                            size_t p_tmp = p1 ^ c;

                            // iterate over partitions of p_tmp
                            for (size_t p2 = 0; p2 <= p_tmp; p2++)
                            {
                                if ( (p2 | p_tmp) == p_tmp)
                                {
                                    size_t p3 = p2 ^ p_tmp;

                                    std::vector<double>::const_iterator         lc = u_l  + p1*this->num_chars;
                                    std::vector<double>::const_iterator         rc = u_r  + p2*this->num_chars;
                                    std::vector<double>::const_iterator         mc = u_m  + p3*this->num_chars;

                                    // iterate over initial states
                                    for (size_t ci = 0; ci < this->num_chars; ci++)
                                    {
                                        uc[ci] += f[ci] * lc[ci] * rc[ci] * mc[ci];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeRootCorrection( size_t root, size_t left, size_t right)
{
    if ( useFastInformative() )
    {
        computeRootCorrectionFast(root, left, right);
        return;
    }

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    std::vector<double>::iterator         p_node  = getMutableCorrectionLikelihoodsForNode(root);
    std::vector<double>::const_iterator   p_left  = getCorrectionLikelihoodsForNode(left);
    std::vector<double>::const_iterator   p_right = getCorrectionLikelihoodsForNode(right);

    // get the root frequencies
    std::vector<std::vector<double> > ff;
    this->getRootFrequencies(ff);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        const std::vector<double> &f = ff[mixture % ff.size()];

        // iterate over correction masks
        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            // iterate over ancestral (non-autapomorphic) states
            for (size_t a = 0; a < this->num_chars; a++)
            {
                size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset + a*correctionOffset;

                std::vector<double>::iterator                 u = p_node   + offset;
                std::vector<double>::const_iterator         u_l = p_left   + offset;
                std::vector<double>::const_iterator         u_r = p_right  + offset;

                // iterate over combinations of autapomorphic states
                for (size_t c = 0; c < numCorrectionPatterns; c++)
                {
                    std::vector<double>::iterator         uc = u  + c*this->num_chars;

                    std::fill(uc, uc + this->num_chars, 0.0);
                    
                    // iterate over partitions of c
                    for (size_t p1 = 0; p1 <= c; p1++)
                    {
                        if ( (p1 | c) == c)
                        {
                            size_t p2 = p1 ^ c;

                            std::vector<double>::const_iterator         lc = u_l  + p1*this->num_chars;
                            std::vector<double>::const_iterator         rc = u_r  + p2*this->num_chars;

                            // iterate over initial states
                            for (size_t ci = 0; ci < this->num_chars; ci++)
                            {
                                uc[ci] += f[ci] * lc[ci] * rc[ci];
                            }
                        }
                    }
                }
            }
        }
    }
}

template<class charType>
double RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::excludedProbForMixture( std::vector<double>::const_iterator p_node, size_t mask, size_t mixture ) const
{
    double prob = 0.0;
    
    if ( useFastInformative() )
    {
        size_t k = this->num_chars;
        size_t nObs = maskObservationCounts[mask];
        std::vector<double>::const_iterator u = p_node + mixture*correctionMixtureOffset + mask*correctionMaskOffset;

        for (size_t j = 0; j < k; j++)
        {
            // weight 1/nObs when every observed taxon is in a different state
            double w = (nObs > j + 1 ? 1.0 : (nObs == j + 1 ? 1.0/nObs : 0.0));
            if (w == 0.0)
                continue;

            double tot = RB_FAST_A(u,j) + j*RB_FAST_B(u,j) + (k - 1 - j)*RB_FAST_C(u,j);
            prob += k * RbMath::choose(k-1, j) * w * tot;
        }
    }
    else
    // iterate over ancestral (non-autapomorphic) states
    for (size_t a = 0; a < this->num_chars; a++)
    {
        size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset + a*correctionOffset;

        std::vector<double>::const_iterator             u = p_node   + offset;

        // iterate over combinations of autapomorphic states
        for (size_t c = 0; c < numCorrectionPatterns; c++)
        {
            std::vector<double>::const_iterator         uc = u  + c*this->num_chars;

            if (c == 0 && coding == AscertainmentBias::VARIABLE)
            {
                // iterate over initial states
                for (size_t ci = 0; ci < this->num_chars; ci++)
                {
                    prob += uc[ci];
                }
            }
            else if (coding == AscertainmentBias::INFORMATIVE)
            {
                // get number of autapomorphic states
                size_t popcount = std::bitset<sizeof(size_t)*CHAR_BIT>(c).count();

               /* if the number of observations is <= the number of autapomorphies then
                * 1. don't count patterns with zero prob (num obs < num auto)
                * 2. don't double-count patterns (num obs == num auto)
                *
                * if num obs == num auto + 1
                * then count the all-autapomorphies pattern only once, when a is its
                * smallest state (bits of c below a are the states below a)
                */
                // ordered: a pattern is uninformative only with at most one autapomorphic
                // state below a and one above, so the all-autapomorphies pattern with one
                // state on each side of a is counted only for that a
                size_t lower = std::bitset<sizeof(size_t)*CHAR_BIT>(c & ((size_t(1) << a) - 1)).count();
                if ( ordered && (lower > 1 || popcount - lower > 1) )
                {
                    continue;
                }

                if ( (maskObservationCounts[mask] == popcount + 1 && (lower == 0 || (ordered && popcount == 2))) ||
                     maskObservationCounts[mask] > popcount + 1)
                {
                    // iterate over initial states
                    for (size_t ci = 0; ci < this->num_chars; ci++)
                    {
                        prob += uc[ci];
                    }
                }
            }
        }
    }
    assert(std::isnan(prob) || (0.0 <= prob and prob <= 1.00000000001));
    prob = std::min(prob, 1.0);
    return prob;
}

template<class charType>
double RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::compatibleExcludedProb( size_t mask ) const
{
    const TopologyNode &root = this->tau->getValue().getRoot();
    std::vector<double>::const_iterator p_node = getCorrectionLikelihoodsForNode( root.getIndex() );
    std::vector<double> mixtureProbs = this->getMixtureProbs();
    double total = 0.0;
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        total += excludedProbForMixture(p_node, mask, mixture) * mixtureProbs[mixture];
    return total;
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeRootLikelihoods( std::vector<double> &rv ) const
{
    PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihoods( rv );

    if ( restrictedEntryOfPattern.empty() || this->getPInv() > 0.0 )
        return;

    for (size_t p = 0; p < rv.size() && p < restrictedEntryOfPattern.size(); ++p)
    {
        size_t entry = restrictedEntryOfPattern[p];
        if (entry == 0)
            continue;

        double removed = compatibleExcludedProb(entry);
        size_t count = this->pattern_counts[p];
        if (removed <= 0.0 || count == 0)
            continue;

        double adjusted = std::exp(rv[p] / double(count)) - removed;
        if (adjusted <= 0.0)
            continue;
        rv[p] = std::log(adjusted) * double(count);
    }
}

template<class charType>
double RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::sumRootLikelihood( void )
{
    double sumPartialProbs = PhyloCTMCSiteHomogeneous<charType>::sumRootLikelihood();
    if (coding == AscertainmentBias::ALL)
        return sumPartialProbs;

    // get the root node
    const TopologyNode &root = this->tau->getValue().getRoot();

    // get the index of the root node
    size_t node_index = root.getIndex();
    
    std::vector<double>::const_iterator p_node = getCorrectionLikelihoodsForNode(node_index);
    
    std::vector<double> perMaskCorrections = std::vector<double>(numCorrectionMasks, 0.0);
    
    std::vector<double> mixtureProbs = this->getMixtureProbs();
    
    // iterate over each correction mask
    for (size_t mask = 0; mask < numCorrectionMasks; mask++)
    {
        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {
            double prob = excludedProbForMixture(p_node, mask, mixture);

//            perMaskCorrections[mask] += prob * 0.25;
//            perMaskMixtureCorrections[mask*this->num_site_mixtures + mixture] = (1.0 - prob) * 0.25;
            perMaskCorrections[mask] += prob * mixtureProbs[mixture];
            perMaskMixtureCorrections[mask*this->num_site_mixtures + mixture] = (1.0 - prob) * mixtureProbs[mixture];
        
            // add corrections for invariant sites
            double prob_invariant = this->getPInv();
            if (prob_invariant > 0.0)
            {
                perMaskCorrections[mask] *= (1.0 - prob_invariant);
                perMaskMixtureCorrections[mask*this->num_site_mixtures + mixture] *= (1.0 - prob_invariant);
            }
            
        }

        // add corrections for invariant sites
        double prob_invariant = this->getPInv();
        if (prob_invariant > 0.0)
        {
            if ( coding & AscertainmentBias::VARIABLE)
            {
                perMaskCorrections[mask] += prob_invariant;
            }
        }

        // normalize the log-probability
//        perMaskCorrections[mask] /= this->num_site_mixtures;
        assert(std::isnan(perMaskCorrections[mask]) || (0.0 <= perMaskCorrections[mask] and perMaskCorrections[mask] <= 1.00000000001));
        perMaskCorrections[mask] = std::min(perMaskCorrections[mask], 1.0);

        perMaskCorrections[mask] = log1p(-perMaskCorrections[mask]);

        // apply the correction for this correction mask
        sumPartialProbs -= perMaskCorrections[mask]*correctionMaskCounts[mask];
    }

    return sumPartialProbs;
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::buildRestrictedEntries( void )
{
    restrictedPatternOfEntry.clear();
    restrictedEntryOfPattern.clear();

    if (coding == AscertainmentBias::ALL || numRealMasks == 0)
        return;

    if (correctionMaskCounts.size() > numRealMasks)
    {
        correctionMaskCounts.resize(numRealMasks);
        maskObservationCounts.resize(numRealMasks);
        correctionMaskMatrix.resize(numRealMasks);
    }
    numCorrectionMasks = numRealMasks;

    size_t tips = this->gap_matrix.size();
    if (this->using_ambiguous_characters == false || tips == 0 || this->pattern_block_size == 0)
        return;

    restrictedEntryOfPattern.assign(this->pattern_block_size, 0);

    for (size_t p = 0; p < this->pattern_block_size; ++p)
    {
        bool needed = false;
        if (coding == AscertainmentBias::VARIABLE)
        {
            // compress() has already found the states every scored cell admits
            needed = ( p < this->site_invariant.size() && this->site_invariant[p] );
        }
        else
        {
            for (size_t t = 0; t < tips && needed == false; ++t)
            {
                if ( this->gap_matrix[t][p] )
                    continue;
                size_t n = this->ambiguous_char_matrix[t][p].count();
                if ( n > 1 && n < this->num_chars )
                    needed = true;
            }
        }
        if (needed == false)
            continue;

        std::vector<bool> maskData(this->num_nodes, false);
        size_t scored = 0;
        for (size_t t = 0; t < tips; ++t)
        {
            bool g = this->gap_matrix[t][p]
                     || this->ambiguous_char_matrix[t][p].count() >= this->num_chars;
            maskData[t] = g;
            if (g == false)
                ++scored;
        }

        restrictedEntryOfPattern[p] = correctionMaskCounts.size();
        restrictedPatternOfEntry.push_back(p);
        correctionMaskCounts.push_back(0);
        maskObservationCounts.push_back(scored);
        correctionMaskMatrix.push_back(maskData);
    }
    numCorrectionMasks = correctionMaskCounts.size();
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::resizeLikelihoodVectors( void )
{
    buildRestrictedEntries();


    RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::resizeLikelihoodVectors();
    
    if (coding != AscertainmentBias::ALL)
    {
        correctionMixtureOffset = numCorrectionMasks*correctionMaskOffset;
        correctionNodeOffset    = this->num_site_mixtures*correctionMixtureOffset;
        activeCorrectionOffset  = this->num_nodes*correctionNodeOffset;

        // this buffer holds 2^k entries per state per node, so it is rebuilt only when its
        // size actually changes
        if (correctionLikelihoods.size() != activeCorrectionOffset*2)
            correctionLikelihoods.assign(activeCorrectionOffset*2, 0.0);

        if (perMaskMixtureCorrections.size() != this->num_site_mixtures*numCorrectionMasks)
            perMaskMixtureCorrections.assign(this->num_site_mixtures*numCorrectionMasks, 0.0);
    }
}


// do not cache an answer formed before the rate matrix is attached
template<class charType>
double RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeLnProbability( void )
{
    fastInformativeCache = -1;
    return PhyloCTMCSiteHomogeneous<charType>::computeLnProbability();
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::touchSpecialization( const DagNode *toucher, bool touchAll )
{
    fastInformativeCache = -1;
    PhyloCTMCSiteHomogeneous<charType>::touchSpecialization(toucher, touchAll);
}


// cached until a touch
template<class charType>
bool RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::useFastInformative( void ) const
{
    // the entries are built after this is first called, and the fast path cannot fill them
    if ( restrictedPatternOfEntry.empty() == false )
        return false;

    if (fastInformativeCache >= 0)
        return fastInformativeCache == 1;

    // no answer while both matrix slots are empty
    if (this->homogeneous_rate_matrix == NULL && this->heterogeneous_rate_matrices == NULL)
        return false;

    fastInformativeCache = 0;

    // setting this falls back to the general recursion, which must give the same answer
    static const bool disabled = (getenv("RB_NO_FAST_INFORMATIVE") != NULL);
    if (disabled)
        return false;

    if (coding != AscertainmentBias::INFORMATIVE || ordered || this->num_chars < 3)
        return false;

    // a Jukes-Cantor generator is the Mk model, whose states are exchangeable. Rev hands a single
    // Q to the vector form, so both have to be looked at
    if (this->homogeneous_rate_matrix != NULL)
    {
        if (dynamic_cast<const RateMatrix_JC *>( &this->homogeneous_rate_matrix->getValue() ) == NULL)
            return false;
    }
    else if (this->heterogeneous_rate_matrices != NULL)
    {
        const RbVector<RateGenerator> &qs = this->heterogeneous_rate_matrices->getValue();
        if (qs.size() == 0)
            return false;
        for (size_t i = 0; i < qs.size(); i++)
            if (dynamic_cast<const RateMatrix_JC *>( &qs[i] ) == NULL)
                return false;
    }
    else
        return false;

    // root frequencies may still be given separately, and must be flat
    std::vector<std::vector<double> > ff;
    this->getRootFrequencies(ff);
    for (size_t i = 0; i < ff.size(); i++)
        for (size_t j = 0; j < ff[i].size(); j++)
            if (fabs(ff[i][j] - 1.0/this->num_chars) > 1e-12)
                return false;

    fastInformativeCache = 1;
    return true;
}


// disjoint singleton sets, so the counts convolve binomially
template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::fastMerge(std::vector<double>::const_iterator l, std::vector<double>::const_iterator r, double *d) const
{
    size_t k = this->num_chars;
    std::fill(d, d + 3*k, 0.0);

    for (size_t j = 0; j < k; j++)
    {
        for (size_t i = 0; i <= j; i++)
        {
            double bin = RbMath::choose(j, i);
            RB_FAST_A(d,j) += bin * RB_FAST_A(l,i) * RB_FAST_A(r,j-i);
            RB_FAST_C(d,j) += bin * RB_FAST_C(l,i) * RB_FAST_C(r,j-i);
            if (j >= 1)
            {
                if (i >= 1)
                    RB_FAST_B(d,j) += RbMath::choose(j-1, i-1) * RB_FAST_B(l,i) * RB_FAST_C(r,j-i);
                if (j - i >= 1)
                    RB_FAST_B(d,j) += RbMath::choose(j-1, i)   * RB_FAST_C(l,i) * RB_FAST_B(r,j-i);
            }
        }
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeTipCorrectionFast(const TopologyNode &node, size_t node_index)
{
    this->updateTransitionProbabilityMatrix(node_index);

    std::vector<double>::iterator p_node = getMutableCorrectionLikelihoodsForNode(node_index);
    size_t k = this->num_chars;
    size_t data_tip_index = this->taxon_name_2_tip_index_map[ node.getName() ];

    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        const TransitionProbabilityMatrix &pij = this->pmatrices[node_index][mixture];
        double same = pij[0][0], diff = pij[1][0];

        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            bool gap = correctionMaskMatrix[mask][data_tip_index];
            std::vector<double>::iterator u = p_node + mixture*correctionMixtureOffset + mask*correctionMaskOffset;
            std::fill(u, u + 3*k, 0.0);

            if (gap)
            {
                // scored in no taxon here, so it uses no singleton and permits any state
                RB_FAST_A(u,0) = 1.0;
                RB_FAST_C(u,0) = 1.0;
            }
            else
            {
                // the tip shows the background, or it is the one taxon showing a singleton
                RB_FAST_A(u,0) = same;
                RB_FAST_C(u,0) = diff;
                RB_FAST_A(u,1) = diff;
                RB_FAST_B(u,1) = same;
                RB_FAST_C(u,1) = diff;
            }
        }
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeInternalNodeCorrectionFast(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{
    this->updateTransitionProbabilityMatrix(node_index);

    std::vector<double>::const_iterator p_left  = getCorrectionLikelihoodsForNode(left);
    std::vector<double>::const_iterator p_right = getCorrectionLikelihoodsForNode(right);
    std::vector<double>::iterator       p_node  = getMutableCorrectionLikelihoodsForNode(node_index);
    size_t k = this->num_chars;
    std::vector<double> d(3*k);

    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        const TransitionProbabilityMatrix &pij = this->pmatrices[node_index][mixture];
        double same = pij[0][0], diff = pij[1][0];

        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset;
            fastMerge(p_left + offset, p_right + offset, &d[0]);

            std::vector<double>::iterator u = p_node + offset;
            for (size_t j = 0; j < k; j++)
            {
                double A = RB_FAST_A(d,j), B = RB_FAST_B(d,j), C = RB_FAST_C(d,j);
                double nOther = (k - 1 > j ? k - 1 - j : 0);
                RB_FAST_A(u,j) = same*A + j*diff*B + nOther*diff*C;
                RB_FAST_B(u,j) = diff*A + same*B + (j >= 1 ? (j-1) : 0)*diff*B + nOther*diff*C;
                RB_FAST_C(u,j) = diff*A + j*diff*B + same*C + (nOther >= 1 ? nOther - 1 : 0)*diff*C;
            }
        }
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeInternalNodeCorrectionFast(const TopologyNode &node, size_t node_index, size_t left, size_t right, size_t middle)
{
    this->updateTransitionProbabilityMatrix(node_index);

    std::vector<double>::const_iterator p_left   = getCorrectionLikelihoodsForNode(left);
    std::vector<double>::const_iterator p_right  = getCorrectionLikelihoodsForNode(right);
    std::vector<double>::const_iterator p_middle = getCorrectionLikelihoodsForNode(middle);
    std::vector<double>::iterator       p_node   = getMutableCorrectionLikelihoodsForNode(node_index);
    size_t k = this->num_chars;
    std::vector<double> lr(3*k), d(3*k);

    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        const TransitionProbabilityMatrix &pij = this->pmatrices[node_index][mixture];
        double same = pij[0][0], diff = pij[1][0];

        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset;
            fastMerge(p_left + offset, p_right + offset, &lr[0]);
            fastMerge(lr.cbegin(), p_middle + offset, &d[0]);

            std::vector<double>::iterator u = p_node + offset;
            for (size_t j = 0; j < k; j++)
            {
                double A = RB_FAST_A(d,j), B = RB_FAST_B(d,j), C = RB_FAST_C(d,j);
                double nOther = (k - 1 > j ? k - 1 - j : 0);
                RB_FAST_A(u,j) = same*A + j*diff*B + nOther*diff*C;
                RB_FAST_B(u,j) = diff*A + same*B + (j >= 1 ? (j-1) : 0)*diff*B + nOther*diff*C;
                RB_FAST_C(u,j) = diff*A + j*diff*B + same*C + (nOther >= 1 ? nOther - 1 : 0)*diff*C;
            }
        }
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeRootCorrectionFast(size_t root, size_t left, size_t right)
{
    std::vector<double>::const_iterator p_left  = getCorrectionLikelihoodsForNode(left);
    std::vector<double>::const_iterator p_right = getCorrectionLikelihoodsForNode(right);
    std::vector<double>::iterator       p_node  = getMutableCorrectionLikelihoodsForNode(root);
    size_t k = this->num_chars;
    std::vector<double> d(3*k);

    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset;
            fastMerge(p_left + offset, p_right + offset, &d[0]);

            std::vector<double>::iterator u = p_node + offset;
            for (size_t i = 0; i < 3*k; i++)
                u[i] = d[i] / k;
        }
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::computeRootCorrectionFast(size_t root, size_t left, size_t right, size_t middle)
{
    std::vector<double>::const_iterator p_left   = getCorrectionLikelihoodsForNode(left);
    std::vector<double>::const_iterator p_right  = getCorrectionLikelihoodsForNode(right);
    std::vector<double>::const_iterator p_middle = getCorrectionLikelihoodsForNode(middle);
    std::vector<double>::iterator       p_node   = getMutableCorrectionLikelihoodsForNode(root);
    size_t k = this->num_chars;
    std::vector<double> lr(3*k), d(3*k);

    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        for (size_t mask = 0; mask < numCorrectionMasks; mask++)
        {
            size_t offset = mixture*correctionMixtureOffset + mask*correctionMaskOffset;
            fastMerge(p_left + offset, p_right + offset, &lr[0]);
            fastMerge(lr.cbegin(), p_middle + offset, &d[0]);

            std::vector<double>::iterator u = p_node + offset;
            for (size_t i = 0; i < 3*k; i++)
                u[i] = d[i] / k;
        }
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::updateCorrections( const TopologyNode& node, size_t nodeIndex )
{

    // only necessary if the node is actually dirty
    //if ( dirtyNodes[nodeIndex] )
    //{
    //    dirtyNodes[nodeIndex] = false;

        // start by filling the likelihood vector for the children of the root
        if ( node.isTip() )
        {
            // this->updateTransitionProbabilities(nodeIndex );
            computeTipCorrection( node, nodeIndex );
        }
        else if ( node.getNumberOfChildren() == 2 ) // rooted trees have two children for the root
        {
            const TopologyNode &left = node.getChild(0);
            size_t leftIndex = left.getIndex();
            updateCorrections( left, leftIndex );
            const TopologyNode &right = node.getChild(1);
            size_t rightIndex = right.getIndex();
            updateCorrections( right, rightIndex );

            if (node.isRoot())
                computeRootCorrection( nodeIndex, leftIndex, rightIndex );
            else
            {
                // this->updateTransitionProbabilities( nodeIndex );
                computeInternalNodeCorrection( node, nodeIndex, leftIndex, rightIndex );
            }

        }
        else if ( node.getNumberOfChildren() == 3 ) // unrooted trees have three children for the root
        {
            const TopologyNode &left = node.getChild(0);
            size_t leftIndex = left.getIndex();
            updateCorrections( left, leftIndex );
            const TopologyNode &right = node.getChild(1);
            size_t rightIndex = right.getIndex();
            updateCorrections( right, rightIndex );
            const TopologyNode &middle = node.getChild(2);
            size_t middleIndex = middle.getIndex();
            updateCorrections( middle, middleIndex );

            if (node.isRoot())
                computeRootCorrection( nodeIndex, leftIndex, rightIndex, middleIndex );
            else
            {
                // this->updateTransitionProbabilities( nodeIndex );
                computeInternalNodeCorrection( node, nodeIndex, leftIndex, rightIndex, middleIndex );
            }

        }
        else
        {
            throw RbException("Node has an unexpected number of children. Only 2 or 3 are allowed.");
        }
    //}
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::redrawValue( void ) {

    if (coding == AscertainmentBias::ALL)
    {
        PhyloCTMCSiteHomogeneous<charType>::redrawValue();
        
        return;
    }

    // delete the old value first
    delete this->value;

    // create a new character data object
    this->value = new HomologousDiscreteCharacterData<charType>();

    size_t num_tips = this->tau->getValue().getNumberOfTips();
    size_t num_nodes = this->tau->getValue().getNumberOfNodes();

    RandomNumberGenerator* rng = GLOBAL_RNG;

    const TopologyNode &root = this->tau->getValue().getRoot();
    size_t root_index = this->tau->getValue().getRoot().getIndex();

    // for safety, we resize the likelihood vectors here
    resizeLikelihoodVectors();
    
    updateCorrections(root, root_index);

    std::vector< DiscreteTaxonData<charType> > taxa = std::vector< DiscreteTaxonData<charType> >(num_tips, DiscreteTaxonData<charType>( Taxon("") ) );

    /*// first sample a total number of characters (M) from the marginal posterior:
    // M - N | N ~ NegBinomial(N+1, exp(lnCorrection) )
    double M_minus_N = RbStatistics::NegativeBinomial::rv(N + 1, exp(perMaskCorrections[0]), *rng);

    // then sample the observed number of characters (num_sites) from the likelihood:
    // num_sites | M ~ Binomial(M, exp(lnCorrection) )
    this->num_sites = RbStatistics::Binomial::rv( M_minus_N + N, exp(perMaskCorrections[0]), *rng);*/

    // sample the rate categories in proportion to the total probability (correction) for each mixture.
    double total = 0.0;
    std::vector<double> mixtureProbs = this->getMixtureProbs();
    for ( size_t i = 0; i < this->num_site_mixtures; ++i )
        total += mixtureProbs[i];

    std::vector<size_t> perSiteRates;
    for ( size_t i = 0; i < this->num_sites; ++i )
    {
        // draw the state
        double u = rng->uniform01()*total;
        size_t rateIndex = 0;

        double tmp = 0.0;
        while(tmp < u)
        {
            tmp += mixtureProbs[rateIndex];
            if (tmp < u)
                rateIndex++;
        }

        perSiteRates.push_back( rateIndex );
    }
    
    // get the root frequencies
    std::vector<std::vector<double> > ff;
    this->getRootFrequencies(ff);

    // then sample site-patterns using rejection sampling,
    // rejecting those that match the unobservable ones.
    for ( size_t i = 0; i < this->num_sites; i++ )
    {
        size_t rateIndex = perSiteRates[i];

        const std::vector< double > &freqs = ff[perSiteRates[i] % ff.size()];

        std::vector<charType> siteData(num_nodes, charType(this->num_chars));

        // create the character
        charType &c = siteData[root_index];
        c.setToFirstState();
        // draw the state
        double u = rng->uniform01();
        size_t stateIndex = 0;
        while ( true )
        {
            u -= freqs[stateIndex];
            ++stateIndex;

            if ( u > 0.0 && stateIndex < this->num_chars)
            {
                ++c;
            }
            else
            {
                break;
            }

        }

        // recursively simulate the sequences
        std::map<size_t, size_t> charCounts;
        simulateConditional( root, siteData, rateIndex, charCounts);

        if ( !isSitePatternCompatible(charCounts) )
        {
            i--;
            continue;
        }

        // add the taxon data to the character data
        for (size_t t = 0; t < num_tips; ++t)
        {
            taxa[t].addCharacter(siteData[t]);
        }
    }

    // add the taxon data to the character data
    for (size_t i = 0; i < this->tau->getValue().getNumberOfTips(); ++i)
    {
        taxa[i].setTaxon( this->tau->getValue().getNode(i).getTaxon() );
        this->value->addTaxonData( taxa[i] );
    }

    // compress the data and initialize internal variables
    this->reInitialized();

    this->markAllPartialLikelihoodsDirty();
    this->pmatrices.mark_all_dirty();

    // flip the active likelihood pointers
    for (size_t index = 0; index < this->changed_nodes.size(); ++index)
    {
        if ( this->changed_nodes[index] == false )
        {
            this->activeLikelihood[index] = (this->activeLikelihood[index] == 0 ? 1 : 0);
            this->changed_nodes[index] = true;
        }
    }

    updateCorrections(root, root_index);
}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::simulateConditional( const TopologyNode &node, std::vector<charType> &data, size_t rateIndex, std::map<size_t, size_t>& charCounts) {

    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();

    // get the sequence of this node
    size_t node_index = node.getIndex();
    charType &parentState = data[ node_index ];

    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (auto child: node.getChildren())
    {
        // update the transition probability matrix
        this->updateTransitionProbabilityMatrix( child->getIndex() );

        std::uint64_t cp = parentState.getStateIndex();

        const double *freqs = this->pmatrices[child->getIndex()][ rateIndex ][ cp ];

        // create the character
        charType &c = data[ child->getIndex() ];
        c.setToFirstState();
        // draw the state
        double u = rng->uniform01();
        size_t stateIndex = 0;
        while ( true )
        {
            u -= *freqs;
            ++stateIndex;

            if ( u > 0.0 && stateIndex < this->num_chars)
            {
                ++c;
                ++freqs;
            }
            else
            {
                break;
            }

        }

        if (child->isTip())
            charCounts[c.getStateIndex()]++;
        else
            simulateConditional( *child, data, rateIndex, charCounts);
    }

}

template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneousConditional<charType>::setValue(AbstractHomologousDiscreteCharacterData *v, bool force)
{
    warned = false;

    // delegate to the parent class
    PhyloCTMCSiteHomogeneous<charType>::setValue(v, force);
}

#endif
