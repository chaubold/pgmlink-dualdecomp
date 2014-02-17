#include "pgmlink/hardconstraintchecker.h"
#include "pgmlink/log.h"

namespace pgmlink {

HardConstraintChecker::HardConstraintChecker()
{
}

void HardConstraintChecker::add_incoming_constraint(const IdList &transition_nodes, size_t detection_node)
{
    incoming_constraints_.push_back(std::make_pair(transition_nodes, detection_node));
}

bool HardConstraintChecker::check_incoming_constraints(const Configuration &config)
{
    // check all listed constraints
    for(std::vector< std::pair<IdList, size_t> >::iterator it = incoming_constraints_.begin();
        it != incoming_constraints_.end();
        ++it)
    {
        size_t detection = it->second;
        IdList& transition_nodes = it->first;

        size_t num_incoming = 0;
        size_t num_detections = config[detection];

        // count the number of incoming transitions = sum_j (Y_ij)
        for(IdList::iterator transition_node = transition_nodes.begin();
            transition_node != transition_nodes.end();
            ++transition_node)
        {
            num_incoming += config[*transition_node];
        }

        // check that as many incoming than detections
        if(num_incoming != num_detections)
        {
            LOG(logWARNING) << "Num detections not equal to num incoming at nodes(values): "
                            << detection << "(" << num_detections << ") != "
                            << num_incoming;
            return false;
        }
    }

    return true;
}

void HardConstraintChecker::add_outgoing_constraint(const IdList &transition_nodes,
                                                    size_t detection_node,
                                                    size_t division_node)
{
    outgoing_constraints_.push_back(boost::make_tuple(transition_nodes, detection_node, division_node));
}


bool HardConstraintChecker::check_outgoing_constraints(const Configuration &config)
{
    // check all listed constraints
    for(std::vector< boost::tuple<IdList, size_t, size_t> >::iterator it = outgoing_constraints_.begin();
        it != outgoing_constraints_.end();
        ++it)
    {
        size_t detection = (*it).get<1>();
        size_t division = (*it).get<2>();
        IdList& transition_nodes = (*it).get<0>();

        size_t num_outgoing = 0;
        size_t num_detections = config[detection];
        size_t num_divisions = config[division];

        // we cannot have more divisions than detections
        if(num_divisions > num_detections)
        {
            LOG(logWARNING) << "More divisions than detections at nodes(values): "
                            << detection << "(" << num_detections << ") < "
                            << division << "(" << num_divisions << ")";
            return false;
        }

        // count the number of outgoing transitions = sum_j (Y_ij)
        for(IdList::iterator transition_node = transition_nodes.begin();
            transition_node != transition_nodes.end();
            ++transition_node)
        {
            num_outgoing += config[*transition_node];
        }

        // check that the detections and divisions are equal to the number of outgoings
        if(num_outgoing != num_detections + num_divisions)
        {
            LOG(logWARNING) << "Outgoing Transition (A + D = sum Y) constraint violated by nodes(value): "
                            << detection << "(" << num_detections << ") + "
                            << division << "(" << num_divisions << ") != "
                            << num_outgoing << "(in " << transition_nodes.size() << " nodes)";
            return false;
        }
    }

    return true;
}

void HardConstraintChecker::add_appearance_disappearance_constraint(size_t appearance_node, size_t disappearance_node)
{
    appearance_disappearance_constraints_.push_back(std::make_pair(appearance_node, disappearance_node));
}

bool HardConstraintChecker::check_appearance_disappearance_constraints(const Configuration &config)
{
    // check all listed constraints
    for(std::vector< std::pair<size_t, size_t> >::iterator it = appearance_disappearance_constraints_.begin();
        it != appearance_disappearance_constraints_.end();
        ++it)
    {
        size_t appearance_node = it->first;
        size_t disappearance_node = it->second;

        size_t num_appearances = config[appearance_node];
        size_t num_disappearances = config[disappearance_node];

        // check that either A = V or A = 0 or V = 0
        if(num_appearances != num_disappearances && num_appearances > 0 && num_disappearances > 0)
        {
            LOG(logWARNING) << "Appearance-Disappearance constraint violated by nodes(value): "
                            << appearance_node << "(" << num_appearances << ") != "
                            << disappearance_node << "(" << num_disappearances << ")";
            return false;
        }
    }

    return true;
}

bool HardConstraintChecker::check_configuration(const Configuration &config)
{
    return check_outgoing_constraints(config)
            && check_incoming_constraints(config)
            && check_appearance_disappearance_constraints(config);
}

} // namespace pgmlink
