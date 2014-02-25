#pragma once
#ifndef DUALDECOMP_BUNDLE_HARDCONSTRAINT_HXX
#define DUALDECOMP_BUNDLE_HARDCONSTRAINT_HXX

#include <vector>
#include <list>
#include <typeinfo>
#include <opengm/inference/inference.hxx>
#include <opengm/inference/visitors/visitor.hxx>
#include <opengm/graphicalmodel/space/discretespace.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/inference/dualdecomposition/dualdecomposition_base.hxx>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

// include conic bundle
#include <CBSolver.hxx>

// pgmlink includes
#include <boost/function.hpp>
#include "../hardconstraintchecker.h"
#include "../log.h"

namespace opengm {

/// Inference based on dual decomposition using bundle methods
template<class GM, class INF, class DUALBLOCK >
class DualDecompositionBundleWithHardConstraints
        : public Inference<GM,typename INF::AccumulationType>,
        public DualDecompositionBase<GM, DUALBLOCK>,
        public ConicBundle::FunctionOracle
{
public:
    typedef GM                                                 GmType;
    typedef GM                                                 GraphicalModelType;
    typedef typename INF::AccumulationType                     AccumulationType;
    OPENGM_GM_TYPE_TYPEDEFS;
    typedef VerboseVisitor<DualDecompositionBundleWithHardConstraints<GM, INF,DUALBLOCK> > VerboseVisitorType;
    typedef TimingVisitor<DualDecompositionBundleWithHardConstraints<GM, INF,DUALBLOCK> >  TimingVisitorType;
    typedef EmptyVisitor<DualDecompositionBundleWithHardConstraints<GM, INF,DUALBLOCK> >   EmptyVisitorType;
    typedef INF                                                InfType;
    typedef DUALBLOCK                                          DualBlockType;
    typedef typename DualBlockType::DualVariableType           DualVariableType;
    typedef DualDecompositionBase<GmType, DualBlockType>       DDBaseType;

    typedef typename DDBaseType::SubGmType                     SubGmType;
    typedef typename DualBlockType::SubFactorType              SubFactorType;
    typedef typename DualBlockType::SubFactorListType          SubFactorListType;
    typedef typename DDBaseType::SubVariableType               SubVariableType;
    typedef typename DDBaseType::SubVariableListType           SubVariableListType;

    typedef boost::function<void(const SubGmType&, size_t, InfType&)> HardConstraintConfigurator;

    class Parameter : public DualDecompositionBaseParameter{
    public:
        /// The relative accuracy which have to be garantee to stop with an approximative solution (set 0 for optimality)
        double minimalRelAccuracy_;
        /// Parameter for Subproblems
        typename InfType::Parameter subPara_;
        /// Relative Precision of dual bound
        double relativeDualBoundPrecision_;
        /// Maximal size of bundle
        size_t maxBundlesize_;
        /// Some variables will be fixed automatically to the center value if their bounds are strongly active (i.e., the corresponding multipliers are big).
        bool activeBoundFixing_;
        /// Lower bound on the weight for the quadratic term of the augmented subproblem.
        double minDualWeight_;
        /// Upper bound on the weight for the quadratic term of the augmented subproblem.
        double maxDualWeight_;
        /// Use a special solver that only employs a minimal bundle consisting of just one new and one aggregate gradient so that there is no real bundle available.
        bool noBundle_;
        /// Uses heuristic for stepsize/trustregion-radius
        bool useHeuristicStepsize_;

        Parameter()
            : relativeDualBoundPrecision_(0.0),
              maxBundlesize_(100),
              activeBoundFixing_(false),
              minDualWeight_(-1),
              maxDualWeight_(-1),
              noBundle_(false),
              useHeuristicStepsize_(true)
        {};
    };

    using  DualDecompositionBase<GmType, DualBlockType >::gm_;
    using  DualDecompositionBase<GmType, DualBlockType >::subGm_;
    using  DualDecompositionBase<GmType, DualBlockType >::dualBlocks_;
    using  DualDecompositionBase<GmType, DualBlockType >::numDualsOvercomplete_;
    using  DualDecompositionBase<GmType, DualBlockType >::numDualsMinimal_;
    using  DualDecompositionBase<GmType, DualBlockType >::modelWithSameVariables_;

    ~DualDecompositionBundleWithHardConstraints();
    DualDecompositionBundleWithHardConstraints(const GmType&,
                                               const HardConstraintConfigurator& hcc,
                                               pgmlink::HardConstraintChecker* hard_constraint_checker);
    DualDecompositionBundleWithHardConstraints(const GmType&,
                                               const Parameter&,
                                               const HardConstraintConfigurator& hcc,
                                               pgmlink::HardConstraintChecker* hard_constraint_checker);
    virtual std::string name() const {return "DualDecompositionBundleWithHardConstraints";};
    virtual const GmType& graphicalModel() const {return gm_;};
    virtual InferenceTermination infer();
    template<class VisitorType>
    InferenceTermination infer(VisitorType&);
    virtual ValueType bound() const;
    virtual ValueType value() const;
    virtual InferenceTermination arg(std::vector<LabelType>&, const size_t = 1)const;
    virtual int evaluate(const ConicBundle::DVector&, double, double&, ConicBundle::DVector&, std::vector<ConicBundle::DVector>&,
                         std::vector<ConicBundle::PrimalData*>&, ConicBundle::PrimalExtender*&);
    virtual DualDecompositionBaseParameter& parameter();


protected:
    // overrride parent method
    template<class ACC> void getBounds(const std::vector<std::vector<LabelType> >&, const std::vector<SubVariableListType>&, ValueType&, ValueType&, std::vector<LabelType>&);

private:
    virtual void allocate();
    int dualStep(const size_t iteration);
    template <class T_IndexType, class T_LabelType>
    void getPartialSubGradient(const size_t, const std::vector<T_IndexType>&, std::vector<T_LabelType>&)const;
    double euclideanSubGradientNorm();

    // Members
    std::vector<std::vector<LabelType> >  subStates_;
    ConicBundle::CBSolver* solver;
    size_t nullStepCounter_;

    Accumulation<ValueType,LabelType,AccumulationType> acUpperBound_;
    Accumulation<ValueType,LabelType,AccumulationType> acNegLowerBound_;
    ValueType upperBound_;
    ValueType lowerBound_;

    Parameter              para_;
    std::vector<ValueType> mem_;
    std::vector<ValueType> mem2_;

    opengm::Timer primalTimer_;
    opengm::Timer dualTimer_;
    double primalTime_;
    double dualTime_;

    HardConstraintConfigurator hardConstraintConfigurator_;
    pgmlink::HardConstraintChecker* hard_constraint_checker_;
};

//**********************************************************************************
template<class GM, class INF, class DUALBLOCK>
DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::~DualDecompositionBundleWithHardConstraints()
{
    delete solver;
}

template<class GM, class INF, class DUALBLOCK>
DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::DualDecompositionBundleWithHardConstraints(const GmType& gm,
                                                                                                         const HardConstraintConfigurator &hcc,
                                                                                                         pgmlink::HardConstraintChecker *hard_constraint_checker)
    : DualDecompositionBase<GmType, DualBlockType >(gm),
      hardConstraintConfigurator_(hcc),
      hard_constraint_checker_(hard_constraint_checker)
{
    this->init(para_);
    subStates_.resize(subGm_.size());
    for(size_t i=0; i<subGm_.size(); ++i)
        subStates_[i].resize(subGm_[i].numberOfVariables());

    solver = new ConicBundle::CBSolver(para_.noBundle_);
    // Setup bundle-solver
    solver->init_problem(numDualsMinimal_);
    solver->add_function(*this);
    solver->set_out(&std::cout,0);//1=output

    solver->set_max_bundlesize(*this,para_.maxBundlesize_);
    //solver->set_eval_limit(1000);
    //solver->set_inner_update_limit(1);
    solver->set_active_bounds_fixing(para_.activeBoundFixing_);
    solver->set_term_relprec(para_.relativeDualBoundPrecision_);
    solver->set_min_weight(para_.minDualWeight_);
    solver->set_max_weight(para_.maxDualWeight_);
    nullStepCounter_ =0;
}

template<class GM, class INF, class DUALBLOCK>
DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::DualDecompositionBundleWithHardConstraints(const GmType& gm,
                                                                                                         const Parameter& para,
                                                                                                         const HardConstraintConfigurator &hcc,
                                                                                                         pgmlink::HardConstraintChecker *hard_constraint_checker)
    :  DualDecompositionBase<GmType, DualBlockType >(gm),
      hardConstraintConfigurator_(hcc),
      hard_constraint_checker_(hard_constraint_checker)
{
    para_ = para;
    this->init(para_);

    subStates_.resize(subGm_.size());
    for(size_t i=0; i<subGm_.size(); ++i)
        subStates_[i].resize(subGm_[i].numberOfVariables());

    solver = new ConicBundle::CBSolver(para_.noBundle_);
    // Setup bundle-solver
    solver->init_problem(numDualsMinimal_);
    solver->add_function(*this);
    solver->set_out(&std::cout,0);//1=output
    solver->set_max_bundlesize(*this,para_.maxBundlesize_);
    //solver->set_eval_limit(1000);
    //solver->set_inner_update_limit(1);
    solver->set_active_bounds_fixing(para.activeBoundFixing_);
    solver->set_term_relprec(para_.relativeDualBoundPrecision_);
    solver->set_min_weight(para_.minDualWeight_);
    solver->set_max_weight(para_.maxDualWeight_);
    nullStepCounter_ =0;
}


////////////////////////////////////////////////////////////////////

template <class GM, class INF, class DUALBLOCK>
void DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::allocate()
{
    mem_.resize(numDualsOvercomplete_,0);
    mem2_.resize(numDualsOvercomplete_,0);
    ValueType *data1Front = &mem_[0];
    ValueType *data1Back  = &mem_[numDualsOvercomplete_];
    ValueType *data2Front = &mem2_[0];
    ValueType *data2Back  = &mem2_[numDualsOvercomplete_];
    for(typename std::vector<DualBlockType>::iterator it=dualBlocks_.begin(); it!=dualBlocks_.end(); ++it){
        for(size_t i=0; i<(*it).duals_.size(); ++i){
            DualVariableType& dv1 = (*it).duals_[i];
            DualVariableType& dv2 = (*it).duals2_[i];
            if(i+1==(*it).duals_.size()){
                data1Back -= dv1.size();
                data2Back -= dv2.size();
                dv1.assign( dv1.shapeBegin(),dv1.shapeEnd(),data1Back);
                dv2.assign( dv2.shapeBegin(),dv2.shapeEnd(),data2Back);
            }
            else{
                dv1.assign( dv1.shapeBegin(),dv1.shapeEnd(),data1Front);
                dv2.assign( dv2.shapeBegin(),dv2.shapeEnd(),data2Front);
                data1Front += dv1.size();
                data2Front += dv2.size();
            }
        }
    }
    OPENGM_ASSERT(data1Front ==  data1Back );
    OPENGM_ASSERT(data2Front ==  data2Back );
    OPENGM_ASSERT(data1Front ==  &mem_[0]+numDualsMinimal_);
    OPENGM_ASSERT(data2Front ==  &mem2_[0]+numDualsMinimal_ );
}

template <class GM, class INF, class DUALBLOCK>
DualDecompositionBaseParameter& DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::parameter()
{
    return para_;
}

/////////////////////////

template<class GM, class INF, class DUALBLOCK>
InferenceTermination DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::
infer()
{
    EmptyVisitorType visitor;
    return infer(visitor);
}

template<class GM, class INF, class DUALBLOCK>
template<class VisitorType>
InferenceTermination DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::
infer(VisitorType& visitor)
{
    std::cout.precision(15);
    visitor.begin(*this,this->value(),this->bound());
    for(size_t iteration=0; iteration<para_.maximalNumberOfIterations_; ++iteration){
        // Dual Step
        dualTimer_.tic();
        int ret;
        if(dualBlocks_.size() == 0){
            // Solve subproblems
            for(size_t subModelId=0; subModelId<subGm_.size(); ++subModelId){
                InfType inf(subGm_[subModelId],para_.subPara_);
                inf.infer();
                inf.arg(subStates_[subModelId]);
            }

            // Calculate lower-bound
            std::vector<LabelType> temp;
            std::vector<LabelType> temp2;
            const std::vector<SubVariableListType>& subVariableLists = para_.decomposition_.getVariableLists();
            (*this).template getBounds<AccumulationType>(subStates_, subVariableLists, lowerBound_, upperBound_, temp);
            acNegLowerBound_(-lowerBound_,temp2);
            acUpperBound_(upperBound_, temp);
            ret=128;
        }
        else{
            ret = dualStep(iteration);
        }
        dualTimer_.toc();
        dualTime_ = dualTimer_.elapsedTime() - primalTime_;
        std::cout.precision(15);
        visitor((*this), upperBound_, lowerBound_);
        //visitor((*this),lowerBound_,-acNegLowerBound_.value(), upperBound_, acUpperBound_.value(), primalTime_, dualTime_);

        dualTime_  = 0;
        primalTime_ = 0;


        // Test for Convergence
        ValueType o;
        AccumulationType::iop(0.0001,-0.0001,o);
        OPENGM_ASSERT(AccumulationType::bop(lowerBound_, upperBound_+o));
        OPENGM_ASSERT(AccumulationType::bop(-acNegLowerBound_.value(), acUpperBound_.value()+o));

        if(   fabs(acUpperBound_.value() + acNegLowerBound_.value())                       <= para_.minimalAbsAccuracy_
              || fabs((acUpperBound_.value()+ acNegLowerBound_.value())/acUpperBound_.value()) <= para_.minimalRelAccuracy_
              || ret ==1){
            visitor.end((*this), acUpperBound_.value(), -acNegLowerBound_.value());
            return NORMAL;
        }
        if(ret>0){
            break;
        }
    }
    visitor.end((*this), acUpperBound_.value(), -acNegLowerBound_.value());
    return NORMAL;
}

template<class GM, class INF, class DUALBLOCK>
InferenceTermination DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::
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
typename GM::ValueType DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::value() const
{
    return acUpperBound_.value();
}

template<class GM, class INF, class DUALBLOCK>
typename GM::ValueType DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::bound() const
{
    return -acNegLowerBound_.value();
}


///////////////////////////////////////////////////////////////

template <class GM, class INF, class DUALBLOCK>
int DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::dualStep(const size_t iteration)
{
    int retval;
    if(para_.useHeuristicStepsize_){
        solver->set_min_weight(para_.minDualWeight_);
        solver->set_max_weight(para_.maxDualWeight_);
    }
    else if(iteration == 0){
        solver->set_min_weight(100);
        solver->set_max_weight(100);
    }
    else{
        if(solver->get_objval() == solver->get_candidate_value() || iteration==1){
            //Serious Step
            double primalDualGap   = fabs(acUpperBound_.value() + acNegLowerBound_.value());

            double subgradientNorm =  (*this).euclideanSubGradientNorm();
            double stepsize = (primalDualGap)/subgradientNorm * para_.stepsizeStride_;

            if(para_.minDualWeight_>=0)
                stepsize = std::min(1/para_.minDualWeight_, stepsize);
            if(para_.maxDualWeight_>=0)
                stepsize = std::max(1/para_.maxDualWeight_, stepsize);

            double t   = 1/stepsize;
            solver->set_next_weight(t);
            solver->set_min_weight(t);
            solver->set_max_weight(t);
            nullStepCounter_ =0;
        }
        else{
            // Null Step
            ++nullStepCounter_;
        }
    }

    retval=solver->do_descent_step(1);

    if (retval){
        std::cout<<"descent_step returned "<<retval<<std::endl;
    }
    //std::cout << solver->get_last_weight() << std::endl;
    return solver->termination_code();
}

template <class GM, class INF, class DUALBLOCK>
int DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::evaluate
(
        const ConicBundle::DVector&            dual, // Argument/Lagrange multipliers
        double                                 relprec,
        double&                                objective_value,
        ConicBundle::DVector&                  cut_vals,
        std::vector<ConicBundle::DVector>&     subgradients,
        std::vector<ConicBundle::PrimalData*>& primal_solutions,
        ConicBundle::PrimalExtender*&          primal_extender
        )
{
    typename std::vector<DualBlockType>::iterator it;
    typename SubFactorListType::const_iterator lit;

    for(size_t i=0; i<numDualsMinimal_; ++i){
        mem_[i] = dual[i];
    }
    for(it = dualBlocks_.begin(); it != dualBlocks_.end(); ++it){
        const size_t numDuals = (*it).duals_.size();
        (*it).duals_[numDuals-1] = -(*it).duals_[0];
        for(size_t i=1; i<numDuals-1;++i){
            (*it).duals_[numDuals-1] -= (*it).duals_[i];
        }
    }
    // Solve Subproblems
    objective_value=0;
    primalTimer_.tic();

    //#ifdef WITH_OPENMP
    //      omp_set_num_threads(para_.numberOfThreads_);
    //#pragma omp parallel for
    //#endif
    for(size_t subModelId=0; subModelId<subGm_.size(); ++subModelId){
        InfType inf(subGm_[subModelId],para_.subPara_);
        inf.infer();
        inf.arg(subStates_[subModelId]);
    }
    primalTimer_.toc();
    primalTime_ +=  primalTimer_.elapsedTime();

    // Calculate lower-bound
    std::vector<LabelType> temp;
    std::vector<LabelType> temp2;
    const std::vector<SubVariableListType>& subVariableLists = para_.decomposition_.getVariableLists();
    (*this).template getBounds<AccumulationType>(subStates_, subVariableLists, lowerBound_, upperBound_, temp);
    acNegLowerBound_(-lowerBound_,temp2);
    acUpperBound_(upperBound_, temp);
    objective_value = -lowerBound_;

    // Store subgradient
    std::vector<size_t> s;
    mem2_.assign(mem2_.size(),0);
    for(it = dualBlocks_.begin(); it != dualBlocks_.end(); ++it){
        const size_t numDuals = (*it).duals_.size();
        lit = (*((*it).subFactorList_)).begin();
        s.resize((*lit).subIndices_.size());
        for(size_t i=0; i<numDuals; ++i){
            getPartialSubGradient((*lit).subModelId_, (*lit).subIndices_, s);
            ++lit;
            (*it).duals2_[i](s.begin()) += -1.0;
        }
        for(size_t i=0; i<numDuals-1; ++i){
            (*it).duals2_[i] -=  (*it).duals2_[numDuals-1] ;
        }
    }

    //construct first subgradient and objective value
    ConicBundle::PrimalDVector h(numDualsMinimal_,0);
    cut_vals.push_back(objective_value);
    for(size_t i=0; i<numDualsMinimal_; ++i){
        h[i] = mem2_[i];
    }
    subgradients.push_back(h);
    //  primal_solutions.push_back(h.clone_primal_data());
    return 0;
}



template <class GM, class INF, class DUALBLOCK>
template <class T_IndexType, class T_LabelType>
inline void DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::getPartialSubGradient
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
double DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::euclideanSubGradientNorm()
{
    double norm = 0;
    std::vector<size_t> s,s2;
    typename std::vector<DUALBLOCK>::const_iterator it;
    typename SubFactorListType::const_iterator                  lit;
    for(it = dualBlocks_.begin(); it != dualBlocks_.end(); ++it){
        const size_t numDuals = (*it).duals_.size();
        const SubFactorType& sf = (*((*it).subFactorList_)).back();
        lit = (*((*it).subFactorList_)).begin();
        s.resize((*lit).subIndices_.size());
        s2.resize((*lit).subIndices_.size());
        getPartialSubGradient(sf.subModelId_, sf.subIndices_, s2);
        for(size_t i=0; i<numDuals-1; ++i){
            getPartialSubGradient((*lit).subModelId_, (*lit).subIndices_, s);
            ++lit;
            if(s!=s2)
                norm += 2;
        }
    }
    return sqrt(norm);
}

template<class GM, class INF, class DUALBLOCK>
template <class ACC>
void DualDecompositionBundleWithHardConstraints<GM,INF,DUALBLOCK>::getBounds
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
                value += 360 * num_violated_constraints;

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

#endif // DUALDECOMP_BUNDLE_HARDCONSTRAINT_HXX
