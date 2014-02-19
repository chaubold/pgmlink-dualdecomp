#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY
#define BOOST_PYTHON_MAX_ARITY 25

#include <vector>

#include "../include/pgmlink/tracking.h"
#include "../include/pgmlink/field_of_view.h"
#include <boost/utility.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python.hpp>
#include <Python.h>


using namespace std;
using namespace pgmlink;
using namespace boost::python;

vector<vector<Event> > pythonChaingraphTracking(ChaingraphTracking& tr, TraxelStore& ts) {
	vector<vector<Event> > result = std::vector<std::vector<Event> >(0);
	// release the GIL
	Py_BEGIN_ALLOW_THREADS
	try {
		result = tr(ts);
	} catch (std::exception& e) {
		Py_BLOCK_THREADS
		throw;
	}
	Py_END_ALLOW_THREADS
	return result;
}

void export_track() {
    class_<vector<Event> >("EventVector")
	.def(vector_indexing_suite<vector<Event> >())
    ;

    class_<vector<vector<Event> > >("NestedEventVector")
	.def(vector_indexing_suite<vector<vector<Event> > >())
    ;

    class_<map<unsigned int, bool> >("DetectionMap")
      .def(map_indexing_suite<map<unsigned int, bool> >())
    ;

    class_<vector<map<unsigned int, bool> > >("DetectionMapsVector")
      .def(vector_indexing_suite<vector<map<unsigned int, bool> > >())
    ;

    class_<ChaingraphTracking>("ChaingraphTracking", 
			       init<string,double,double,double,double,
			       	   bool,double,double,bool,
			       	   bool,double,double,double,double>(
							  args("random_forest_filename", "appearance", "disappearance", "detection", "misdetection",
									  "use_random_forest", "opportunity_cost", "forbidden_cost", "with_constraints",
									  "fixed_detections", "mean_div_dist", "min_angle", "ep_gap", "n_neighbors"
									  )))
      .def("__call__", &pythonChaingraphTracking)
      .def("detections", &ChaingraphTracking::detections)
      .def("set_with_divisions", &ChaingraphTracking::set_with_divisions)
      .def("set_cplex_timeout", &ChaingraphTracking::set_cplex_timeout)
    ;

    class_<ConsTracking>("ConsTracking",
                         init<int,double,double,string,bool,double,double,double,bool,double,double,bool,double,double, bool, int, double, double, FieldOfView, bool>(
						args("max_number_objects", "max_neighbor_distance", "division_threshold",
							"detection_rf_filename", "size_dependent_detection_prob", "forbidden_cost",
							"ep_gap", "avg_obj_size",
							"with_tracklets",
							"division_weight", "transition_weight",
							"with_divisions",
							 "disappearance_cost", "appearance_cost", "with_merger_resolution", "number_of_dimensions",
							 "transition_parameter", "border_width", "fov", "with_constraints")))
	  .def("__call__", &ConsTracking::operator())
	  .def("detections", &ConsTracking::detections)
	;

    class_<ConsTrackingDD>("ConsTrackingDD",
                         init<int,double,double,string,bool,
                         double,double,double,bool,double,double,
                         bool,double,double, bool, int, double,
                         double, FieldOfView, bool, int, int>(
                        args("max_number_objects", "max_neighbor_distance", "division_threshold",
                            "detection_rf_filename", "size_dependent_detection_prob", "forbidden_cost",
                            "ep_gap", "avg_obj_size",
                            "with_tracklets",
                            "division_weight", "transition_weight",
                            "with_divisions",
                             "disappearance_cost", "appearance_cost", "with_merger_resolution",
                             "number_of_dimensions",
                             "transition_parameter", "border_width", "fov", "with_constraints",
                             "timesteps_per_block", "num_overlapping_timesteps")))
      .def("__call__", &ConsTrackingDD::operator())
      .def("detections", &ConsTrackingDD::detections)
    ;

    enum_<Event::EventType>("EventType")
	.value("Move", Event::Move)
	.value("Division", Event::Division)
	.value("Appearance", Event::Appearance)
	.value("Disappearance", Event::Disappearance)
	.value("Merger", Event::Merger)
	.value("MultiFrameMove", Event::MultiFrameMove)
	.value("Void", Event::Void)
    ;

    class_<vector<unsigned int> >("IdVector")
	.def(vector_indexing_suite<vector<unsigned int> >())
    ;
    class_<Event>("Event")
	.def_readonly("type", &Event::type)
	.def_readonly("traxel_ids", &Event::traxel_ids)
	.add_property("energy", &Event::energy)
    ;
}
