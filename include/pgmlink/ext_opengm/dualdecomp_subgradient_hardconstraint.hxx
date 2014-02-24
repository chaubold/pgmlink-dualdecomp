#ifndef OPENGM_DUALDECOMP_SUBGRADIENT_HARDCONSTRAINT_HXX
#define OPENGM_DUALDECOMP_SUBGRADIENT_HARDCONSTRAINT_HXX

#pragma once

#include <vector>
#include <list>
#include <typeinfo>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include <opengm/inference/inference.hxx>
#include <opengm/inference/visitors/visitor.hxx>
#include <opengm/graphicalmodel/space/discretespace.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/inference/dualdecomposition/dualdecomposition_base.hxx>
#include <opengm/inference/dualdecomposition/dualdecomposition_subgradient.hxx>

#include <boost/function.hpp>
#include "../hardconstraintchecker.h"
#include "../log.h"

namespace opengm {

/// Inference based on dual decomposition using sub-gradient descent
/// Tailored to suit our tracking needs -> copies hard constraints to subproblems
template<class GM, class INF, class DUALBLOCK >
class DualDecompositionSubGradientWithHardConstraints
        : public Inference<GM,typename INF::AccumulationType>,  public DualDecompositionBase<GM, DUALBLOCK >
{
public:
    typedef GM                                                 GmType;
    typedef GM                                                 GraphicalModelType;
    typedef typename INF::AccumulationType                     AccumulationType;
    OPENGM_GM_TYPE_TYPEDEFS;
    typedef VerboseVisitor<DualDecompositionSubGradientWithHardConstraints<GM, INF,DUALBLOCK> > VerboseVisitorType;
    typedef TimingVisitor<DualDecompositionSubGradientWithHardConstraints<GM, INF,DUALBLOCK> >  TimingVisitorType;
    typedef EmptyVisitor<DualDecompositionSubGradientWithHardConstraints<GM, INF,DUALBLOCK> >   EmptyVisitorType;

    typedef INF                                                InfType;
    typedef DUALBLOCK                                          DualBlockType;
    typedef DualDecompositionBase<GmType, DualBlockType>       DDBaseType;

    typedef typename DualBlockType::DualVariableType           DualVariableType;
    typedef typename DDBaseType::SubGmType                     SubGmType;
    typedef typename DualBlockType::SubFactorType              SubFactorType;
    typedef typename DualBlockType::SubFactorListType          SubFactorListType;
    typedef typename DDBaseType::SubVariableType               SubVariableType;
    typedef typename DDBaseType::SubVariableListType           SubVariableListType;

    typedef boost::function<void(const SubGmType&, size_t, InfType&)> HardConstraintConfigurator;

    class Parameter : public DualDecompositionBaseParameter{
    public:
        /// Parameter for Subproblems
        typename InfType::Parameter subPara_;
        bool useAdaptiveStepsize_;
        bool useProjectedAdaptiveStepsize_;
        Parameter() : useAdaptiveStepsize_(false), useProjectedAdaptiveStepsize_(false){};
    };

    using  DualDecompositionBase<GmType, DualBlockType >::gm_;
    using  DualDecompositionBase<GmType, DualBlockType >::subGm_;
    using  DualDecompositionBase<GmType, DualBlockType >::dualBlocks_;
    using  DualDecompositionBase<GmType, DualBlockType >::numDualsOvercomplete_;
    using  DualDecompositionBase<GmType, DualBlockType >::numDualsMinimal_;
    using  DualDecompositionBase<GmType, DualBlockType >::modelWithSameVariables_;

    DualDecompositionSubGradientWithHardConstraints(const GmType&, const HardConstraintConfigurator& hcc, pgmlink::HardConstraintChecker* hard_constraint_checker);
    DualDecompositionSubGradientWithHardConstraints(const GmType&, const Parameter&, const HardConstraintConfigurator& hcc, pgmlink::HardConstraintChecker* hard_constraint_checker);
    virtual std::string name() const {return "DualDecompositionSubGradientWithHardConstraints";};
    virtual const GmType& graphicalModel() const {return gm_;};
    virtual InferenceTermination infer();
    template <class VISITOR> InferenceTermination infer(VISITOR&);
    virtual ValueType bound() const;
    virtual ValueType value() const;
    virtual InferenceTermination arg(std::vector<LabelType>&, const size_t = 1)const;
    virtual DualDecompositionBaseParameter& parameter();

private:
    virtual void allocate();
    void dualStep(const size_t);
    template <class T_IndexType, class T_LabelType>
    void getPartialSubGradient(const size_t, const std::vector<T_IndexType>&, std::vector<T_LabelType>&)const;
    double euclideanProjectedSubGradientNorm();

protected:
    // overrride parent method
    template<class ACC> void getBounds(const std::vector<std::vector<LabelType> >&, const std::vector<SubVariableListType>&, ValueType&, ValueType&, std::vector<LabelType>&);

private:
    // Members
    std::vector<std::vector<LabelType> >  subStates_;

    Accumulation<ValueType,LabelType,AccumulationType> acUpperBound_;
    Accumulation<ValueType,LabelType,AccumulationType> acNegLowerBound_;
    ValueType upperBound_;
    ValueType lowerBound_;
    std::vector<ValueType> values_;

    Parameter              para_;
    std::vector<ValueType> mem_;

    opengm::Timer primalTimer_;
    opengm::Timer dualTimer_;
    double primalTime_;
    double dualTime_;
    HardConstraintConfigurator hardConstraintConfigurator_;
    pgmlink::HardConstraintChecker* hard_constraint_checker_;
};

//**********************************************************************************
template<class GM, class INF, class DUALBLOCK>
DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::DualDecompositionSubGradientWithHardConstraints(const GmType& gm, const HardConstraintConfigurator& hcc, pgmlink::HardConstraintChecker* hard_constraint_checker)
    : DualDecompositionBase<GmType, DualBlockType >(gm),
      hardConstraintConfigurator_(hcc),
      hard_constraint_checker_(hard_constraint_checker)
{
    this->init(para_);
    subStates_.resize(subGm_.size());
    for(size_t i=0; i<subGm_.size(); ++i)
        subStates_[i].resize(subGm_[i].numberOfVariables());
}

template<class GM, class INF, class DUALBLOCK>
DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::DualDecompositionSubGradientWithHardConstraints(const GmType& gm, const Parameter& para, const HardConstraintConfigurator& hcc, pgmlink::HardConstraintChecker *hard_constraint_checker)
    :   DualDecompositionBase<GmType, DualBlockType >(gm),para_(para),
      hardConstraintConfigurator_(hcc),
      hard_constraint_checker_(hard_constraint_checker)
{
    this->init(para_);
    subStates_.resize(subGm_.size());
    for(size_t i=0; i<subGm_.size(); ++i)
        subStates_[i].resize(subGm_[i].numberOfVariables());
}


////////////////////////////////////////////////////////////////////

template <class GM, class INF, class DUALBLOCK>
void DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::allocate()
{
    if(DDIsView<DualVariableType>::isView()){
        mem_.resize(numDualsOvercomplete_,0.0);
    }
    else
        mem_.resize(1,0.0);
    //std::cout << mem_.size() <<std::flush;
    ValueType *data = &mem_[0];
    for(typename std::vector<DualBlockType>::iterator it=dualBlocks_.begin(); it!=dualBlocks_.end(); ++it){
        for(size_t i=0; i<(*it).duals_.size(); ++i){
            DualVariableType& dv = (*it).duals_[i];
            DualVarAssign(dv, data);
            if(DDIsView<DualVariableType>::isView())
                data += dv.size();
        }
    }
}
template <class GM, class INF, class DUALBLOCK>
DualDecompositionBaseParameter& DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::parameter()
{
    return para_;
}

/////////////////////////
template<class GM, class INF, class DUALBLOCK>
InferenceTermination DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::
infer()
{
    EmptyVisitorType visitor;
    return infer(visitor);
}
template<class GM, class INF, class DUALBLOCK>
template<class VISITOR>
InferenceTermination DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::
infer(VISITOR& visitor)
{
    std::cout.precision(15);
    //visitor.startInference();
    visitor.begin(*this,this->value(),this->bound());

    srand (time(NULL));

    opengm::Timer overallTimer;
    overallTimer.tic();

    LOG(pgmlink::logINFO) << "Number of Subproblems: " << subGm_.size();
    for(size_t iteration=0; iteration<para_.maximalNumberOfIterations_; ++iteration){
        LOG(pgmlink::logDEBUG) << "Iteration: " << iteration;
        // Solve Subproblems
        primalTime_=0;
        primalTimer_.tic();
        //omp_set_num_threads(para_.numberOfThreads_);
        const std::vector<SubVariableListType>& subVariableLists = para_.decomposition_.getVariableLists();

        //#pragma omp parallel for
        for(size_t subModelId=0; subModelId<subGm_.size(); ++subModelId){
            InfType inf(subGm_[subModelId],para_.subPara_);
            hardConstraintConfigurator_(subGm_[subModelId], subModelId, inf);
            inf.setStartingPoint(subStates_[subModelId].begin());
            inf.infer();
            inf.arg(subStates_[subModelId]);

            std::stringstream s;
            s << "Submodel " << subModelId << " - Found solution:\n";

            for(size_t varId=0; varId<gm_.numberOfVariables(); ++varId){
                for(typename SubVariableListType::const_iterator its = subVariableLists[varId].begin();
                    its!=subVariableLists[varId].end();++its){
                    if(subModelId != (*its).subModelId_)
                    {
                        continue;
                    }

                    const size_t& subVariableId = (*its).subVariableId_;
                    s << "(" << varId << ")=" << subStates_[subModelId][subVariableId] << " ";
                }
            }

            LOG(pgmlink::logDEBUG2) << s.str();
        }
        primalTimer_.toc();

        dualTimer_.tic();
        // Calculate lower-bound
        std::vector<LabelType> temp;
        std::vector<LabelType> temp2;
        (*this).template getBounds<AccumulationType>(subStates_, subVariableLists, lowerBound_, upperBound_, temp);
        acNegLowerBound_(-lowerBound_,temp2);
        acUpperBound_(upperBound_, temp);

        //dualStep
        double stepsize = 1.0 / (iteration + 1.0);
        if(para_.useAdaptiveStepsize_){
            stepsize = para_.stepsizeStride_ * fabs(acUpperBound_.value() - lowerBound_)
                    /(*this).subGradientNorm(1);
        }
        else if(para_.useProjectedAdaptiveStepsize_){
            double subgradientNorm = euclideanProjectedSubGradientNorm();
            stepsize = para_.stepsizeStride_ * fabs(acUpperBound_.value() - lowerBound_)
                    /subgradientNorm/subgradientNorm;
        }
        else{
            double primalDualGap = fabs(acUpperBound_.value() + acNegLowerBound_.value());
            double subgradientNorm = euclideanProjectedSubGradientNorm();//(*this).subGradientNorm(1);
            stepsize = para_.getStepsize(iteration,primalDualGap,subgradientNorm);
            //            std::cout << stepsize << std::endl;
        }

        if(typeid(AccumulationType) == typeid(opengm::Minimizer))
            stepsize *=1;
        else
            stepsize *= -1;

        std::vector<size_t> s;
        for(typename std::vector<DualBlockType>::const_iterator it = dualBlocks_.begin(); it != dualBlocks_.end(); ++it){
            const size_t numDuals = (*it).duals_.size();
            typename SubFactorListType::const_iterator lit = (*((*it).subFactorList_)).begin();
            s.resize((*lit).subIndices_.size());

            std::vector<double> weights;
            weights.push_back(1.0 * rand() / RAND_MAX); // random number between 0 and 1

            for(size_t i = 1; i < numDuals; ++i)
            {
                weights.push_back((1.0 - weights[0]) / (numDuals-1)*numDuals); // random weight for the other duals
            }

            weights[0] *= numDuals; // scale all weights such that their sum equals the number of duals

            for(size_t i=0; i<numDuals; ++i){
                getPartialSubGradient<size_t>((*lit).subModelId_, (*lit).subIndices_, s);

                ++lit;
                (*it).duals_[i](s.begin()) += stepsize;
                for(size_t j=0; j<numDuals; ++j){
                    (*it).duals_[j](s.begin()) -= weights[j] * stepsize/numDuals;
                }
            }
            (*it).test();
        }
        dualTimer_.toc();

        primalTime_ = primalTimer_.elapsedTime();
        dualTime_   = dualTimer_.elapsedTime();
        visitor((*this), upperBound_, lowerBound_);
        //visitor((*this), lowerBound_, -acNegLowerBound_.value(), upperBound_, acUpperBound_.value(), primalTime_, dualTime_);

        LOG(pgmlink::logINFO) << "************************************\n";
        LOG(pgmlink::logINFO) << "Current gap: " << fabs(acUpperBound_.value() + acNegLowerBound_.value()) << "\nSubgradient (" << s.size() << "): ";

        std::stringstream str;
        for(std::vector<size_t>::iterator it = s.begin(); it != s.end(); ++it)
        {
            str << *it << " ";
        }
        str << "\nBounds: [" << lowerBound_ << ", " << upperBound_ << "]" << std::endl;
        LOG(pgmlink::logINFO) << str.str();

        LOG(pgmlink::logINFO) << "************************************" << std::endl;

        // Test for Convergence
        ValueType o;
        AccumulationType::iop(0.0001,-0.0001,o);
        OPENGM_ASSERT(AccumulationType::bop(lowerBound_, upperBound_+o));
        OPENGM_ASSERT(AccumulationType::bop(-acNegLowerBound_.value(), acUpperBound_.value()+o));

        if(   fabs(acUpperBound_.value() + acNegLowerBound_.value())                       <= para_.minimalAbsAccuracy_
              || fabs((acUpperBound_.value()+ acNegLowerBound_.value())/acUpperBound_.value()) <= para_.minimalRelAccuracy_){
            visitor.end((*this), acUpperBound_.value(), -acNegLowerBound_.value());
            LOG(pgmlink::logINFO) << "Inference took: " << overallTimer.elapsedTime();
            return NORMAL;
        }
    }
    visitor.end((*this), acUpperBound_.value(), -acNegLowerBound_.value());
    LOG(pgmlink::logINFO) << "Inference took: " << overallTimer.elapsedTime();
    return NORMAL;
}


template<class GM, class INF, class DUALBLOCK>
InferenceTermination DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::
arg(std::vector<LabelType>& conf, const size_t n)const
{
    if(n!=1){
        return UNKNOWN;
    }
    else{
        acUpperBound_.state(conf);
        return NORMAL;
    }
}

template<class GM, class INF, class DUALBLOCK>
typename GM::ValueType DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::value() const
{
    return acUpperBound_.value();
}

template<class GM, class INF, class DUALBLOCK>
typename GM::ValueType DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::bound() const
{
    return -acNegLowerBound_.value();
}


///////////////////////////////////////////////////////////////

template <class GM, class INF, class DUALBLOCK>
void DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::dualStep(const size_t iteration)
{

}

template <class GM, class INF, class DUALBLOCK>
template <class T_IndexType, class T_LabelType>
inline void DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::getPartialSubGradient
(
        const size_t                             subModelId,
        const std::vector<T_IndexType>&    subIndices,
        std::vector<T_LabelType> &                 s
        )const
{
    OPENGM_ASSERT(subIndices.size() == s.size());
    for(size_t n=0; n<s.size(); ++n){
        s[n] = subStates_[subModelId][subIndices[n]];
    }
}

template <class GM, class INF, class DUALBLOCK>
double DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::euclideanProjectedSubGradientNorm()
{
    double norm = 0;
    std::vector<LabelType> s;
    marray::Marray<double> M;
    typename std::vector<DUALBLOCK>::const_iterator it;
    typename SubFactorListType::const_iterator                  lit;
    for(it = dualBlocks_.begin(); it != dualBlocks_.end(); ++it){
        const size_t numDuals = (*it).duals_.size();
        marray::Marray<double> M( (*it).duals_[0].shapeBegin(), (*it).duals_[0].shapeEnd() ,0.0);
        lit = (*((*it).subFactorList_)).begin();
        s.resize((*lit).subIndices_.size());
        for(size_t i=0; i<numDuals; ++i){
            getPartialSubGradient((*lit).subModelId_, (*lit).subIndices_, s);
            ++lit;
            M(s.begin()) += 1;
        }
        for(size_t i=0; i<M.size(); ++i){
            norm  += M(i)            * pow(1.0 - M(i)/numDuals,2);
            norm  += (numDuals-M(i)) * pow(      M(i)/numDuals,2);
        }
    }
    return sqrt(norm);
}

template<class GM, class INF, class DUALBLOCK>
template <class ACC>
void DualDecompositionSubGradientWithHardConstraints<GM,INF,DUALBLOCK>::getBounds
(
   const std::vector<std::vector<LabelType> >& subStates,
   const std::vector<SubVariableListType>& subVariableLists,
   ValueType& lowerBound,
   ValueType& upperBound,
   std::vector<LabelType> & upperState
   )
{
   // Calculate lower-bound
   lowerBound=0;
   for(size_t subModelId=0; subModelId<subGm_.size(); ++subModelId){
       ValueType subModelEnergy = subGm_[subModelId].evaluate(subStates[subModelId]);
       LOG(pgmlink::logDEBUG1) << "SubModel " << subModelId << " has energy: " << subModelEnergy;
      lowerBound += subModelEnergy;
   }

   // Calculate upper-bound
   Accumulation<ValueType,LabelType,ACC> ac;

   // Set modelWithSameVariables_
   if(modelWithSameVariables_[0] == Tribool::Maybe){
      for(size_t varId=0; varId<gm_.numberOfVariables(); ++varId){
         for(typename SubVariableListType::const_iterator its = subVariableLists[varId].begin();
             its!=subVariableLists[varId].end();++its){
            const size_t& subModelId    = (*its).subModelId_;
            const size_t& subVariableId = (*its).subVariableId_;
            if(subVariableId != varId){
               modelWithSameVariables_[subModelId] = Tribool::False;
            }
         }
      }
      for(size_t subModelId=0; subModelId<subGm_.size(); ++subModelId){
         if(gm_.numberOfVariables() != subGm_[subModelId].numberOfVariables()){
            modelWithSameVariables_[subModelId] = Tribool::False;
         }
         if(modelWithSameVariables_[subModelId] == Tribool::Maybe){
            modelWithSameVariables_[subModelId] = Tribool::True;
         }
      }
   }

   // Build Primal-Candidates
   std::vector<std::vector<LabelType> > args(subGm_.size());
   bool somethingToFill = false;
   for(size_t subModelId=0; subModelId<subGm_.size(); ++subModelId){
      if(modelWithSameVariables_[subModelId] == Tribool::False){
         args[subModelId].assign(gm_.numberOfVariables(),std::numeric_limits<LabelType>::max());
         somethingToFill = true;
      }
      else{
          std::stringstream s;
          s << "Accumulating primal candidates - evaluating GM with submodel " << subModelId << " solution:" << std::endl;
          for(size_t i = 0; i < subStates[subModelId].size(); ++i)
          {
              s << subStates[subModelId][i] << " ";
          }
          LOG(pgmlink::logDEBUG2) << s.str();
         ac(gm_.evaluate(subStates[subModelId]),subStates[subModelId]);
      }
   }

   int min_num_violated = std::numeric_limits<int>::max();

   if(somethingToFill){
      for(size_t varId=0; varId<gm_.numberOfVariables(); ++varId){
         for(typename SubVariableListType::const_iterator its = subVariableLists[varId].begin();
             its!=subVariableLists[varId].end();++its){
            const size_t& subModelId    = (*its).subModelId_;
            const size_t& subVariableId = (*its).subVariableId_;
            if(modelWithSameVariables_[subModelId] == Tribool::False){
               args[subModelId][varId] = subStates[subModelId][subVariableId];
            }
         }
         for(size_t subModelId=0; subModelId<subGm_.size(); ++subModelId){
            if(modelWithSameVariables_[subModelId] == Tribool::False &&
               args[subModelId][varId] == std::numeric_limits<LabelType>::max())
            {
               const size_t& aSubModelId    = subVariableLists[varId].front().subModelId_;
               const size_t& aSubVariableId = subVariableLists[varId].front().subVariableId_;
               args[subModelId][varId] = subStates[aSubModelId][aSubVariableId];
            }
         }
      }
      for(size_t subModelId=0; subModelId<subGm_.size(); ++subModelId){
         if(modelWithSameVariables_[subModelId] == Tribool::False){
             ValueType value = gm_.evaluate(args[subModelId]);
             int num_violated_constraints = hard_constraint_checker_->check_configuration(args[subModelId]);
             value += 1000000 * num_violated_constraints;

             min_num_violated = std::min(min_num_violated, num_violated_constraints);

             std::stringstream s;
             s << "Accumulating energy - evaluating GM with submodel " << subModelId << " solution: " << value << std::endl;
             for(size_t i = 0; i < args[subModelId].size(); ++i)
             {
                 s << args[subModelId][i] << " ";
             }
             LOG(pgmlink::logDEBUG2) << s.str();
            ac(value,args[subModelId]);
         }
      }
   }

   LOG(pgmlink::logINFO) << "Min number violated constraints: " << min_num_violated;

   upperBound = ac.value();
   ac.state(upperState);
}


}

#endif // OPENGM_DUALDECOMP_SUBGRADIENT_HARDCONSTRAINT_HXX
