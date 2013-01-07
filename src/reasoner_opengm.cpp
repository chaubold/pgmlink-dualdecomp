#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <opengm/inference/lpcplex.hxx>
#include <opengm/datastructures/marray/marray.hxx>

#include "pgmlink/graphical_model.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_opengm.h"
#include "pgmlink/traxels.h"

//#include <ostream>

using namespace std;

namespace pgmlink {
  Chaingraph::~Chaingraph() {
   if(mrf_ != NULL) {
	delete mrf_;
	mrf_ = NULL;
    }
    if(optimizer_ != NULL) {
	delete optimizer_;
	optimizer_ = NULL;
    }
  }

double Chaingraph::forbidden_cost() const {
    return forbidden_cost_;
}

bool Chaingraph::with_constraints() const {
    return with_constraints_;
}

void Chaingraph::formulate( const HypothesesGraph& hypotheses ) {
    LOG(logDEBUG) << "Chaingraph::formulate: entered";
    reset();
    mrf_ = new OpengmModel<>::ogmGraphicalModel();

    LOG(logDEBUG) << "Chaingraph::formulate: add_detection_nodes";
    add_detection_nodes( hypotheses );
    LOG(logDEBUG) << "Chaingraph::formulate: add_transition_nodes";
    add_transition_nodes( hypotheses );
    LOG(logDEBUG) << "Chaingraph::formulate: add_finite_factors";
    add_finite_factors( hypotheses );

    typedef opengm::LPCplex<OpengmModel<>::ogmGraphicalModel, OpengmModel<>::ogmAccumulator> cplex_optimizer;
    cplex_optimizer::Parameter param;
    param.verbose_ = true;
    param.integerConstraint_ = true;
    param.epGap_ = ep_gap_;
    LOG(logDEBUG) << "Chaingraph::formulate ep_gap = " << param.epGap_;

    optimizer_ = new cplex_optimizer(*mrf_, param);

	if (with_constraints_) {
		LOG(logDEBUG) << "Chaingraph::formulate: add_constraints";
		add_constraints(hypotheses);
	}

	if (fixed_detections_) {
		LOG(logDEBUG) << "Chaingraph::formulate: fix_detections";
		fix_detections(hypotheses, 1);
	}
}



void Chaingraph::infer() {
    opengm::InferenceTermination status = optimizer_->infer();
    if(status != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): optimizer terminated unnormally");
    }
}


void Chaingraph::conclude( HypothesesGraph& g ) {
    // extract solution from optimizer
    vector<OpengmModel<>::ogmInference::LabelType> solution;
    opengm::InferenceTermination status = optimizer_->arg(solution);
    if(status != opengm::NORMAL) {
	throw runtime_error("GraphicalModel::infer(): solution extraction terminated unnormally");
    }

    // add 'active' properties to graph
    g.add(node_active()).add(arc_active());
    property_map<node_active, HypothesesGraph::base_graph>::type& active_nodes = g.get(node_active());
    property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = g.get(arc_active());

    // write state after inference into 'active'-property maps
    for(std::map<HypothesesGraph::Node, size_t>::const_iterator it = node_map_.begin(); it != node_map_.end(); ++it) {
	bool state = false;
	if(solution[it->second] == 1) state = true;
	active_nodes.set(it->first, state);
    }
    for(std::map<HypothesesGraph::Arc, size_t>::const_iterator it = arc_map_.begin(); it != arc_map_.end(); ++it) {
	bool state = false;
	if(solution[it->second] == 1) state = true;
	active_arcs.set(it->first, state);
    }
}

  const OpengmModel<>::ogmGraphicalModel* Chaingraph::get_graphical_model() const {
    return mrf_;
  }

  const std::map<HypothesesGraph::Node, size_t>& Chaingraph::get_node_map() const {
    return node_map_;
  }

  const std::map<HypothesesGraph::Arc, size_t>& Chaingraph::get_arc_map() const {
    return arc_map_;
  }

void Chaingraph::reset() {
    if(mrf_ != NULL) {
	delete mrf_;
	mrf_ = NULL;
    }
    if(optimizer_ != NULL) {
	delete optimizer_;
	optimizer_ = NULL;
    }
    node_map_.clear();
    arc_map_.clear();
}

void Chaingraph::add_detection_nodes( const HypothesesGraph& g) {
    for(HypothesesGraph::NodeIt n(g); n!=lemon::INVALID; ++n) {
	mrf_->addVariable(2);
	node_map_[n] = mrf_->numberOfVariables() - 1; 
    }
}
void Chaingraph::add_transition_nodes( const HypothesesGraph& g) {
    for(HypothesesGraph::ArcIt a(g); a!=lemon::INVALID; ++a) {
	mrf_->addVariable(2);
	arc_map_[a] = mrf_->numberOfVariables() - 1; 
    }
}

void Chaingraph::add_finite_factors( const HypothesesGraph& g) {
  LOG(logDEBUG) << "Chaingraph::add_finite_factors: entered";
  property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());		
  ////
  //// add detection factors
  ////
  LOG(logDEBUG) << "Chaingraph::add_finite_factors: add detection factors";
  for(HypothesesGraph::NodeIt n(g); n!=lemon::INVALID; ++n) {
    size_t vi[] = {node_map_[n]};
    vector<size_t> coords(1,0);
    OpengmExplicitFactor<double> table( vi, vi+1 );
    coords[0] = 0;
    table.set_value( coords, non_detection_(traxel_map[n]) ); 
    coords[0] = 1;
    table.set_value( coords, detection_(traxel_map[n]) );

    table.add_to( *mrf_ );
  }
  
  ////
  //// add transition factors
  ////
  for(HypothesesGraph::NodeIt n(g); n!=lemon::INVALID; ++n) {
    add_outgoing_factor(g, n);
    add_incoming_factor(g, n);
  }
}	    

namespace {
    inline size_t cplex_id(size_t opengm_id) {
	return 2*opengm_id + 1;
    }
}

void Chaingraph::add_constraints( const HypothesesGraph& g ) {
    LOG(logDEBUG) << "Chaingraph::add_constraints: entered";
    typedef opengm::LPCplex<OpengmModel<>::ogmGraphicalModel, OpengmModel<>::ogmAccumulator> cplex;
    ////
    //// outgoing transitions
    ////
    LOG(logDEBUG) << "Chaingraph::add_constraints: outgoing transitions";
    for(HypothesesGraph::NodeIt n(g); n!=lemon::INVALID; ++n) {
	// couple detection and transitions
	for(HypothesesGraph::OutArcIt a(g, n); a!=lemon::INVALID; ++a) {
	    couple(n, a);
	}
	// couple transitions
	vector<size_t> cplex_idxs;
	for(HypothesesGraph::OutArcIt a(g, n); a!=lemon::INVALID; ++a) {
	    cplex_idxs.push_back(cplex_id(arc_map_[a]));
	}
	vector<int> coeffs(cplex_idxs.size(), 1);
	// 0 <= 1*transition + ... + 1*transition <= 2
	dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 2);
    }

    ////
    //// incoming transitions
    ////
    LOG(logDEBUG) << "Chaingraph::add_constraints: incoming transitions";
    for(HypothesesGraph::NodeIt n(g); n!=lemon::INVALID; ++n) {
	// couple detection and transitions
	for(HypothesesGraph::InArcIt a(g, n); a!=lemon::INVALID; ++a) {
	    couple(n, a);
	}
	    
	// couple transitions
	vector<size_t> cplex_idxs;
	for(HypothesesGraph::InArcIt a(g, n); a!=lemon::INVALID; ++a) {
	    cplex_idxs.push_back(cplex_id(arc_map_[a]));
	}
	vector<int> coeffs(cplex_idxs.size(), 1);
	// 0 <= 1*transition + ... + 1*transition <= 1
	dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 1);
    }
}

void Chaingraph::couple(HypothesesGraph::Node& n, HypothesesGraph::Arc& a) {
	    typedef opengm::LPCplex<OpengmModel<>::ogmGraphicalModel, OpengmModel<>::ogmAccumulator> cplex;
    	    vector<size_t> cplex_idxs; 
	    cplex_idxs.push_back(cplex_id(node_map_[n]));
	    cplex_idxs.push_back(cplex_id(arc_map_[a]));
	    vector<int> coeffs;
	    coeffs.push_back(1);
	    coeffs.push_back(-1);
	    // 0 <= 1*detection - 1*transition <= 1
	    dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin() , 0, 1);
}

  void Chaingraph::fix_detections( const HypothesesGraph& g, size_t val ) {
	    typedef opengm::LPCplex<OpengmModel<>::ogmGraphicalModel, OpengmModel<>::ogmAccumulator> cplex;
	    for(HypothesesGraph::NodeIt n(g); n!=lemon::INVALID; ++n) {
	      vector<size_t> cplex_idxs; 
	      cplex_idxs.push_back(cplex_id(node_map_[n]));
	      vector<int> coeffs;
	      coeffs.push_back(1);
	      // val <= 1*detection <= val
	      dynamic_cast<cplex*>(optimizer_)->addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin() , val, val);
	    }
}

  void Chaingraph::add_outgoing_factor( const HypothesesGraph& g, const HypothesesGraph::Node& n ) {
    LOG(logDEBUG) << "Chaingraph::add_outgoing_factor(): entered";
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    // collect and count outgoing arcs
    vector<HypothesesGraph::Arc> arcs; 
    vector<size_t> vi; 		// opengm variable indeces
    vi.push_back(node_map_[n]); // first detection node, remaining will be transition nodes
    int count = 0;
    for(HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
      arcs.push_back(a);
      vi.push_back(arc_map_[a]);
      ++count;
    }

    // construct factor
    if(count == 0) {
      // build value table
      size_t table_dim = 1; 		// only one detection var
      std::vector<size_t> coords;
      OpengmExplicitFactor<double> table( vi );

      // opportunity
      coords = std::vector<size_t>(table_dim, 0); 		// (0)
      table.set_value( coords, opportunity_cost_ );

      // disappearance 
      coords = std::vector<size_t>(table_dim, 0);
      coords[0] = 1; 						// (1)
      table.set_value( coords, disappearance_(traxel_map[n]) );

      table.add_to( *mrf_ );

    } else if(count == 1) {
      // no division possible
      size_t table_dim = 2; 		// detection var + 1 * transition var
      std::vector<size_t> coords;
      OpengmExplicitFactor<double> table( vi, forbidden_cost_ );

      // opportunity configuration
      coords = std::vector<size_t>(table_dim, 0); // (0,0)
      table.set_value( coords, opportunity_cost_ );

      // disappearance configuration
      coords = std::vector<size_t>(table_dim, 0);
      coords[0] = 1; // (1,0)
      table.set_value( coords, disappearance_(traxel_map[n]) );

      // move configurations
      coords = std::vector<size_t>(table_dim, 1);
      // (1,1)
      table.set_value( coords, move_(traxel_map[n], traxel_map[g.target(arcs[0])]) );

      table.add_to( *mrf_ );      

    } else {
      // build value table
      size_t table_dim = count + 1; 		// detection var + n * transition var
      std::vector<size_t> coords;
      OpengmExplicitFactor<double> table( vi, forbidden_cost_ );

      // opportunity configuration
      coords = std::vector<size_t>(table_dim, 0); // (0,0,...,0)
      table.set_value( coords, opportunity_cost_ );

      // disappearance configuration
      coords = std::vector<size_t>(table_dim, 0);
      coords[0] = 1; // (1,0,...,0)
      table.set_value( coords, disappearance_(traxel_map[n]) );

      // move configurations
      coords = std::vector<size_t>(table_dim, 0);
      coords[0] = 1;
      // (1,0,0,0,1,0,0)
      for(size_t i = 1; i < table_dim; ++i) {
	coords[i] = 1; 
	table.set_value( coords, move_(traxel_map[n], traxel_map[g.target(arcs[i-1])]) );
	coords[i] = 0; // reset coords
      }
      
      // division configurations
      coords = std::vector<size_t>(table_dim, 0);
      coords[0] = 1;
      // (1,0,0,1,0,1,0,0) 
      for(unsigned int i = 1; i < table_dim - 1; ++i) {
	for(unsigned int j = i+1; j < table_dim; ++j) {
	  coords[i] = 1;
	  coords[j] = 1;
	  table.set_value(coords, division_(traxel_map[n],
					    traxel_map[g.target(arcs[i-1])],
					    traxel_map[g.target(arcs[j-1])]
					    ));
	  
	  // reset
  	  coords[i] = 0;
	  coords[j] = 0;
	}
      }

      table.add_to( *mrf_ );      

    }   
    LOG(logDEBUG) << "Chaingraph::add_outgoing_factor(): leaving";
}

  void Chaingraph::add_incoming_factor( const HypothesesGraph& g, const HypothesesGraph::Node& n ) {
    LOG(logDEBUG) << "Chaingraph::add_incoming_factor(): entered";
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    // collect and count incoming arcs
    vector<size_t> vi; // opengm variable indeces
    int count = 0;
    for(HypothesesGraph::InArcIt a(g, n); a != lemon::INVALID; ++a) {
      vi.push_back(arc_map_[a]);
      ++count;
    }
    vi.push_back(node_map_[n]); 
    std::reverse(vi.begin(), vi.end());
    
    //// construct factor
    // build value table
    size_t table_dim = count + 1; // detection var + n * transition var
    OpengmExplicitFactor<double> table( vi, forbidden_cost_ );
    std::vector<size_t> coords;

    // allow opportunity configuration
    // (0,0,...,0)
    coords = std::vector<size_t>(table_dim, 0);
    table.set_value( coords, 0 );

    // appearance configuration
    coords = std::vector<size_t>(table_dim, 0);
    coords[0] = 1; // (1,0,...,0)
    table.set_value( coords, appearance_(traxel_map[n]) );
    
    // allow move configurations
    coords = std::vector<size_t>(table_dim, 0);
    coords[0] = 1;
    // (1,0,0,0,1,0,0)
    for(size_t i = 1; i < table_dim; ++i) {
      coords[i] = 1; 
      table.set_value( coords, 0 );
      coords[i] = 0; // reset coords
    }

    table.add_to( *mrf_ );
    LOG(logDEBUG) << "Chaingraph::add_incoming_factor(): leaving";
  }

} /* namespace pgmlink */ 
