#include <opengm/inference/lpcplex.hxx>
#include <opengm/graphicalmodel/graphicalmodel_factor.hxx>
#include <lemon/graph_to_eps.h>
#include "pgmlink/hypotheses.h"
#include <lemon/connectivity.h>
#include <lemon/adaptors.h>
#include "pgmlink/reasoner_constracking_dd.h"

pgmlink::DualDecompositionConservationTracking::DualDecompositionConservationTracking(unsigned int max_number_objects,
        boost::function<double (const Traxel&, const size_t)> detection,
        boost::function<double (const Traxel&, const size_t)> division,
        boost::function<double (const double)> transition,
        double forbidden_cost,
        double ep_gap,
        bool with_tracklets,
        bool with_divisions,
        boost::function<double (const Traxel&)> disappearance_cost_fn,
        boost::function<double (const Traxel&)> appearance_cost_fn,
        bool with_misdetections_allowed,
        bool with_appearance,
        bool with_disappearance,
        double transition_parameter,
        bool with_constraints,
        size_t timesteps_per_block,
        size_t num_overlapping_timesteps)
    : pgmlink::ConservationTracking(
          max_number_objects,
          detection,
          division,
          transition,
          forbidden_cost,
          ep_gap,
          with_tracklets,
          with_divisions,
          disappearance_cost_fn,
          appearance_cost_fn,
          with_misdetections_allowed,
          with_appearance,
          with_disappearance,
          transition_parameter,
          with_constraints),
      hypotheses_graph_(NULL),
      dd_optimizer_(NULL),
      timesteps_per_block_(timesteps_per_block),
      num_overlapping_timesteps_(num_overlapping_timesteps),
      first_dd_iteration_(true)
{}

pgmlink::DualDecompositionConservationTracking::~DualDecompositionConservationTracking()
{

}

void pgmlink::DualDecompositionConservationTracking::add_constraint(
        pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::InfType* optimizer,
        size_t sub_gm_index,
        std::vector<std::size_t>::iterator ids_begin,
        std::vector<std::size_t>::iterator ids_end,
        std::vector<int>::iterator coeffs_begin,
        int lower, int higher, const char* name)
{
    for(std::vector<std::size_t>::iterator it = ids_begin; it != ids_end; ++it)
    {
        if(*it == (size_t)-1)
        {
            std::stringstream s;
            s << "------ discarding constraint between cplex nodes: ";
            for(std::vector<std::size_t>::iterator i = ids_begin; i != ids_end; ++i)
            {
                s << *i << " ";
            }
            LOG(pgmlink::logDEBUG4) << s.str();
            LOG(pgmlink::logDEBUG4) << "------ " << name;

            return;
        }
    }

    std::stringstream s;
    s << "++++++ adding constraint between cplex nodes: ";
    for(std::vector<std::size_t>::iterator it = ids_begin; it != ids_end; ++it)
    {
        s << *it << " ";
    }
    LOG(pgmlink::logDEBUG4) << s.str();
    LOG(pgmlink::logDEBUG4) << "++++++ " << name << std::endl;

    optimizer->addConstraint(ids_begin, ids_end, coeffs_begin, lower, higher, name);
}

void pgmlink::DualDecompositionConservationTracking::constraint_debug_output(
        const pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::SubGmType& subGM)
{
    lemon::InDegMap<HypothesesGraph> degree(*hypotheses_graph_);

    // compare subGM to hypotheses_graph
    for(size_t i = 0; i < std::max(subGM.numberOfVariables(), pgm_->Model()->numberOfVariables()); i++)
    {
        if(i < subGM.numberOfVariables())
        {
            size_t num_factors_opengm = subGM.numberOfFactors(i);
            std::cout << "\tVariable " << i << " has " << num_factors_opengm
                      << " factors in OpenGM->Decomposition" << std::endl;
        }

        if(i < pgm_->Model()->numberOfVariables())
        {
            size_t num_factors_pgm = pgm_->Model()->numberOfFactors(i);
            std::cout << "\tVariable " << i << " has " << num_factors_pgm
                      << " factors in PGM" << std::endl;
        }
    }

    for(size_t i = 0; i < std::max(subGM.numberOfFactors(), pgm_->Model()->numberOfFactors()); i++)
    {
        if(i < subGM.numberOfFactors())
        {
            std::cout << "\tFactor " << i << " has " << subGM[i].numberOfVariables()
                      << " variables in OpenGM->Decomposition\n\t\t";

            for(size_t v = 0; v < subGM[i].numberOfVariables(); v++)
            {
                std::cout << subGM[i].variableIndex(v) << " ";
            }
            std::cout << std::endl;
        }

        if(i < pgm_->Model()->numberOfFactors())
        {
            std::cout << "\tFactor " << i << " has " << (*(pgm_->Model()))[i].numberOfVariables()
                      << " variables in PGM\n\t\t";

            for(size_t v = 0; v < (*(pgm_->Model()))[i].numberOfVariables(); v++)
            {
                std::cout << (*(pgm_->Model()))[i].variableIndex(v) << " ";
            }

            std::cout << std::endl;
        }
    }

    exit(0);
}

void pgmlink::DualDecompositionConservationTracking::configure_hard_constraints(
        const pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::SubGmType& subGM,
        size_t sub_gm_index,
        pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::InfType& optimizer)
{
    if(!first_dd_iteration_)
    {
        hard_constraint_checker_.disable_adding_constraints();
    }
    else
    {
        first_dd_iteration_ = false;
    }
    assert(hypotheses_graph_ != NULL);

    LOG(logINFO) << "Found subproblem with " << subGM.numberOfVariables() << " variables";
    current_sub_gm_id_ = sub_gm_index;
    current_sub_optimizer_ = &optimizer;

    // add constraint
    pgmlink::ConservationTracking::add_constraints(*hypotheses_graph_, boost::bind(
                                                       &pgmlink::DualDecompositionConservationTracking::add_constraint,
                                                       this, &optimizer, sub_gm_index, _1, _2, _3, _4, _5, _6));

    // constraint_debug_output(subGM);
}

void pgmlink::DualDecompositionConservationTracking::debug_graph_output(GraphicalModelType* model)
{
    // count factor orders
    std::map<size_t, size_t> factorOrders;

    std::cout << "Model has " << model->numberOfFactors() << " factors and "
              << model->numberOfVariables() << " variables" << std::endl;
    for(size_t i = 0; i < model->numberOfFactors(); i++)
    {
        size_t order = (*model)[i].numberOfVariables();
        if(factorOrders.find(order) == factorOrders.end())
        {
            factorOrders[order] = 1;
        }
        else
            factorOrders[order]++;
    }

    std::cout << "Factor Variable Count:" << std::endl;
    for(std::map<size_t, size_t>::iterator it = factorOrders.begin();
        it != factorOrders.end();
        ++it)
    {
        std::cout << "\t[" << it->first << "] = " << it->second << std::endl;
    }

    // count variable order
    std::map<size_t, size_t> variableOrders;

    for(size_t i = 0; i < model->numberOfVariables(); i++)
    {
        size_t order = model->numberOfFactors(i);

        if(variableOrders.find(order) == variableOrders.end())
        {
            variableOrders[order] = 1;
        }
        else
            variableOrders[order]++;
    }

    std::cout << "Variable Factor Count:" << std::endl;
    for(std::map<size_t, size_t>::iterator it = variableOrders.begin();
        it != variableOrders.end();
        ++it)
    {
        std::cout << "\t[" << it->first << "] = " << it->second << std::endl;
    }

    std::cout << "Graph decomposed!" << std::endl;
    exit(0);
}

void pgmlink::DualDecompositionConservationTracking::count_connected_components(GraphicalModelType* model)
{
//    // create a partition of all connected components
//    opengm::Partition<size_t> connectedComponents(model->numberOfVariables());
//    // iterate over all factors
//    for(size_t i = 0; i < model->numberOfFactors(); i++) {
//        // iterate over all connected variables of factor and merge them to one partition
//        const GraphicalModelType::ConstVariableIterator variablesBegin = model->variablesOfFactorBegin(i);
//        const GraphicalModelType::ConstVariableIterator variablesEnd = model->variablesOfFactorEnd(i);
//        OPENGM_ASSERT(variablesBegin != variablesEnd);

//        for(GraphicalModelType::ConstVariableIterator iter = variablesBegin + 1; iter != variablesEnd; iter++) {
//            connectedComponents.merge(*(iter - 1), *iter);
//        }
//    }

//    LOG(pgmlink::logINFO) << "Number of Connected Components: " << connectedComponents.numberOfSets();

    // hypotheses graph connected components
    lemon::Undirector<const HypothesesGraph> undirected_hypotheses_graph(*hypotheses_graph_);
    lemon::Undirector<const HypothesesGraph>::NodeMap<int> connected_components(undirected_hypotheses_graph);

    int num_connected_components = lemon::connectedComponents(undirected_hypotheses_graph, connected_components);
    LOG(pgmlink::logINFO) << "Number of Connected Components: " << num_connected_components;

    {
        std::map<int, int> component_cardinality;
        for(lemon::Undirector<const HypothesesGraph>::NodeIt n(undirected_hypotheses_graph); n != lemon::INVALID; ++n)
        {
            component_cardinality[connected_components[n]]++;
        }

        for(std::map<int, int>::iterator it = component_cardinality.begin(); it != component_cardinality.end(); ++it)
        {
            LOG(pgmlink::logINFO) << "\tFound component " << it->first << " containing " << it->second << " nodes";
        }
    }

    // check for bridge nodes in hypotheses graph
    lemon::Undirector<const HypothesesGraph>::NodeMap<bool> node_filter(undirected_hypotheses_graph);
    lemon::Undirector<const HypothesesGraph>::EdgeMap<bool> edge_filter(undirected_hypotheses_graph);
    lemon::SubGraph< const lemon::Undirector<const HypothesesGraph> > sub_graph(undirected_hypotheses_graph,
                                                                                node_filter, edge_filter);

    for(lemon::Undirector<const HypothesesGraph>::EdgeIt e(undirected_hypotheses_graph); e != lemon::INVALID; ++e)
    {
        sub_graph.enable(e);
    }
    for(lemon::Undirector<const HypothesesGraph>::NodeIt n(undirected_hypotheses_graph); n != lemon::INVALID; ++n)
    {
        sub_graph.enable(n);
    }

    for(lemon::Undirector<const HypothesesGraph>::EdgeIt e(undirected_hypotheses_graph); e != lemon::INVALID; ++e)
    {
        sub_graph.disable(e);
        int new_num_connected_components = lemon::countConnectedComponents(sub_graph);

        if(new_num_connected_components > num_connected_components)
        {
            num_connected_components = new_num_connected_components;
            //LOG(pgmlink::logINFO) << "Increased number of Connected Components to: " << num_connected_components;
        }
        else
        {
            sub_graph.enable(e);
        }
    }

    LOG(pgmlink::logINFO) << "Finished checking for bridges, now we have " << num_connected_components << " components";

    {
        lemon::SubGraph< const lemon::Undirector<const HypothesesGraph> >::NodeMap<int> connected_components(sub_graph);

        num_connected_components = lemon::connectedComponents(sub_graph, connected_components);

        std::map<int, int> component_cardinality;
        for(lemon::SubGraph< const lemon::Undirector<const HypothesesGraph> >::NodeIt n(sub_graph); n != lemon::INVALID; ++n)
        {
            component_cardinality[connected_components[n]]++;
        }

        float average_size = 0;
        for(std::map<int, int>::iterator it = component_cardinality.begin(); it != component_cardinality.end(); ++it)
        {
            LOG(pgmlink::logINFO) << "\tFound component " << it->first << " containing " << it->second << " nodes";
            average_size += it->second;
        }
        average_size /= num_connected_components;

        LOG(pgmlink::logINFO) << "Average component size: " << average_size << " of " << num_connected_components << " components";
    }

    std::exit(0);
}

void pgmlink::DualDecompositionConservationTracking::decompose_graph(
        GraphicalModelType* model,
        DualDecompositionSubGradient::Parameter& dd_parameter)
{
    opengm::GraphicalModelDecomposer<GraphicalModelType> decomposer;
    opengm::GraphicalModelDecomposer<GraphicalModelType>::DecompositionType decomposition(
                model->numberOfVariables(),model->numberOfFactors(),0);

    const size_t num_time_steps = nodes_by_timestep_.size();
    size_t num_sub_models = num_time_steps / timesteps_per_block_;

    if(num_time_steps % timesteps_per_block_ != 0)
    {
        num_sub_models++;
    }

    for(size_t sub_model_id = 0; sub_model_id < num_sub_models; ++sub_model_id)
    {
        decomposition.addSubModel();

        std::vector<size_t> sub_variable_map(model->numberOfVariables(),
                                             std::numeric_limits<std::size_t>::max());

        // using overlap, the size of each submodel becomes
        // timestep_per_block + 2 * (overlapping_frames - 1)
        // but the first and last submodel only need one overlap
        size_t first_timestep = sub_model_id * timesteps_per_block_;
        if(first_timestep > num_overlapping_timesteps_)
        {
            first_timestep -= num_overlapping_timesteps_ - 1;
        }

        size_t last_timestep = std::min((sub_model_id + 1) * timesteps_per_block_
                                        + num_overlapping_timesteps_ - 1, num_time_steps);

        // add all variables to their submodels
        for(size_t timestep = first_timestep; timestep < last_timestep; ++timestep)
        {
            std::vector<size_t>& nodes_at_timestep = nodes_by_timestep_[timestep];

            for(std::vector<size_t>::iterator node = nodes_at_timestep.begin();
                node != nodes_at_timestep.end();
                ++node)
            {
                // only add, if we do not have this node in that submodel yet
                if(sub_variable_map[*node] == std::numeric_limits<std::size_t>::max())
                {
                    sub_variable_map[*node] = decomposition.addSubVariable(sub_model_id, *node);
                    LOG(pgmlink::logDEBUG2) << "Adding node: " << *node << " to submodel: " << sub_model_id;
                }
            }
        }

        // add factors to the submodels in which all linked variables are found
        for(size_t factor_id = 0; factor_id < model->numberOfFactors(); ++factor_id)
        {
            if((*model)[factor_id].numberOfVariables() == 0)
            {
               std::vector<size_t> sub_variable_indices(0);
               decomposition.addSubFactor(sub_model_id, factor_id, sub_variable_indices);
            }
            else
            {
                if((*model)[factor_id].numberOfVariables() == 1 && sub_model_id == 0)
                {
                    LOG(pgmlink::logDEBUG2) << "Found unary for variable: "
                              << (*model)[factor_id].variableIndex(0);
                }

                std::vector<size_t> sub_variable_indices((*model)[factor_id].numberOfVariables());
                bool all_variables_inside_submodel = true;

                std::stringstream factor_output;
                factor_output << "Adding factor: " << factor_id << " to submodel " << sub_model_id << " with factors: ";

                for(size_t i = 0;
                    i < (*model)[factor_id].numberOfVariables() && all_variables_inside_submodel;
                    ++i) {
                    const size_t var_id = (*model)[factor_id].variableIndex(i);
                    sub_variable_indices[i] = sub_variable_map[var_id];
                    factor_output << var_id << " ";

                    bool found = (sub_variable_map[var_id] != std::numeric_limits<std::size_t>::max());

                    all_variables_inside_submodel = all_variables_inside_submodel && found;
                }

                if(all_variables_inside_submodel)
                {
                    decomposition.addSubFactor(sub_model_id, factor_id, sub_variable_indices);
                    LOG(pgmlink::logDEBUG2) << factor_output.str();
                }
            }
        }
    }
    decomposition.reorder();
    LOG(pgmlink::logDEBUG2) << "done reordering, completing... ";
    decomposition.complete();

    //dd_parameter.decomposition_ = decomposer.decomposeIntoClosedBlocks(*model, 2);
    dd_parameter.decomposition_ = decomposition;
    LOG(pgmlink::logINFO) << "Decomposed into " << dd_parameter.decomposition_.numberOfSubModels()
              << " submodels";

    for(unsigned int i = 0; i < dd_parameter.decomposition_.numberOfSubModels(); i++)
    {
        LOG(pgmlink::logINFO) << "\tSubproblem " << i << ": "
                  << dd_parameter.decomposition_.numberOfSubFactors(i) << " factors and "
                  << dd_parameter.decomposition_.numberOfSubVariables(i) << " variables and "
                  << dd_parameter.decomposition_.getEmptyFactorLists().size() << " empty factors";
    }


    if(!dd_parameter.decomposition_.isValid(*model))
    {
        LOG(pgmlink::logERROR) << "Model decomposition is invalid!!!!!!!!!";
    }
    else
    {
        LOG(pgmlink::logDEBUG2) << "Model decomposition valid!";
    }
}

void pgmlink::DualDecompositionConservationTracking::infer()
{
    if (!with_constraints_) {
        opengm::hdf5::save(optimizer_->graphicalModel(), "./conservationTracking.h5", "conservationTracking");
        throw std::runtime_error("GraphicalModel::infer(): inference with soft constraints is not implemented yet. "
                                 "The conservation tracking factor graph has been saved to file");
    }

    DualDecompositionSubGradient::Parameter dd_parameter;
    GraphicalModelType* model = pgm_->Model();

    LOG(pgmlink::logINFO) << "Beginning Graph Decomposition";
    LOG(pgmlink::logDEBUG) << "Original Graph had: " << model->numberOfFactors() << " factors and "
              << model->numberOfVariables() << " variables";

//    count_connected_components(model);
    decompose_graph(model, dd_parameter);

//    debug_graph_output(model);

    dd_parameter.decompositionId_ = DualDecompositionSubGradient::Parameter::MANUAL;
    dd_parameter.maximalDualOrder_ = 1;
    dd_parameter.minimalAbsAccuracy_ = 0.0001;
    dd_parameter.subPara_.verbose_ = true;
    dd_parameter.subPara_.integerConstraint_ = true;
    dd_parameter.subPara_.epGap_ = ep_gap_;

    dd_optimizer_ = new DualDecompositionSubGradient(*model, dd_parameter,
                        boost::bind(&pgmlink::DualDecompositionConservationTracking::configure_hard_constraints,
                                    this, _1, _2, _3),
                                                     &hard_constraint_checker_);

    DualDecompositionSubGradient::VerboseVisitorType visitor;
    opengm::InferenceTermination status = dd_optimizer_->infer(visitor);

    if (status != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): optimizer terminated abnormally");
    }
}

void pgmlink::DualDecompositionConservationTracking::extractSolution(
        std::vector<pgm::OpengmModelDeprecated::ogmInference::LabelType> &solution)
{
    opengm::InferenceTermination status = dd_optimizer_->arg(solution);
    if (status != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): solution extraction terminated abnormally");
    }

    std::stringstream s;
    s << "Dual Decomposition - Found solution: ";
    for(size_t i = 0; i < solution.size(); i++)
    {
        s << solution[i] << " ";
    }
    LOG(pgmlink::logINFO) << s.str();
}

void pgmlink::DualDecompositionConservationTracking::add_constraints(const pgmlink::HypothesesGraph &g)
{
    hypotheses_graph_ = &g;
}

size_t pgmlink::DualDecompositionConservationTracking::cplex_id(size_t opengm_id, size_t state)
{
    // map indices to subproblem indices
    typedef opengm::GraphicalModelDecomposition::SubVariableListType SubVarListType;
    const std::vector<SubVarListType>& sub_variable_list =
            dd_optimizer_->parameter().decomposition_.getVariableLists();

    std::size_t new_id = -1; // init with max id

    for(SubVarListType::const_iterator sub_var_it = sub_variable_list[opengm_id].begin();
        sub_var_it != sub_variable_list[opengm_id].end();
        ++sub_var_it)
    {
        if(sub_var_it->subModelId_ == current_sub_gm_id_)
        {
            // use remapping
            new_id = sub_var_it->subVariableId_;
        }
    }

    if(new_id == (size_t)-1)
    {
        // return an "error" value
        LOG(pgmlink::logDEBUG4) << "OpenGM ID " << opengm_id
                                << " does not exist in subproblem " << current_sub_gm_id_;
        return -1;
    }

    // add constraint if the indices are present in the current subproblem
    return current_sub_optimizer_->lpNodeVi(new_id, state);
}
