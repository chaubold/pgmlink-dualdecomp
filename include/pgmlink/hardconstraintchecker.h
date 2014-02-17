#ifndef HARDCONSTRAINTCHECKER_H
#define HARDCONSTRAINTCHECKER_H

#include <iostream>
#include <vector>
#include <boost/tuple/tuple.hpp>

namespace pgmlink {

typedef std::vector<size_t> IdList;
typedef std::vector<size_t> Configuration;

class HardConstraintChecker
{
public:
    HardConstraintChecker();

    void add_incoming_constraint(
            const IdList& transition_nodes,
            size_t detection_node
            );
    void add_outgoing_constraint(
            const IdList& transition_nodes,
            size_t detection_node,
            size_t division_node
            );

    void add_appearance_disappearance_constraint(
            size_t appearance_node,
            size_t disappearance_node);

    int check_configuration(const Configuration& config);
    void disable_adding_constraints();
private:
    int check_outgoing_constraints(const Configuration& config);
    int check_incoming_constraints(const Configuration& config);
    int check_appearance_disappearance_constraints(const Configuration& config);

private:
    std::vector< std::pair<IdList, size_t> > incoming_constraints_;
    std::vector< boost::tuple<IdList, size_t, size_t> > outgoing_constraints_;
    std::vector< std::pair<size_t, size_t> > appearance_disappearance_constraints_;

    bool adding_constraints_enabled_;
};

} // namespace pgmlink
#endif // HARDCONSTRAINTCHECKER_H
