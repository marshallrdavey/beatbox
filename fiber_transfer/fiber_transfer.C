// adapted from libMesh example, miscellaneous/ex5
#include <iostream>
#include <string>
#include "math.h"

// LibMesh include files.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/elem.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/exodusII_io_helper.h"
#include "libmesh/fe_interface.h"
#include "libmesh/getpot.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/discontinuity_measure.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/boundary_info.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/point_locator_tree.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/mesh_function.h"
#include "libmesh/exact_solution.h"
#include "libmesh/system.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// function to parse strings of integers
void parse_ID_string(std::vector<int>& IDs, const std::string& IDs_string)
{
    std::string blah;
    char foo;
    for(int ii = 0; ii < IDs_string.size(); ++ii)
    {
        foo = IDs_string.at(ii);
        if(foo != ',')
        {
            blah = blah + foo;
        }
        else
        {
            IDs.push_back(atoi(blah.c_str()));
            blah.clear();
        }
    }
    IDs.push_back(atoi(blah.c_str()));
 }

//main function
int main (int argc, char** argv)
{
	//initialize libmesh
	LibMeshInit init(argc, argv);

	//load the input file
	GetPot input_file(argv[1]);

	//read in parameters from the input file
	const unsigned int dim              	    = input_file("dimension",3);
	const unsigned int mesh_order		    = input_file("mesh_order",1);
	const std::string ref_mesh_name      	    = input_file("ref_mesh_name","");
	const std::string new_mesh_name             = input_file("new_mesh_name","");
	const libMesh::subdomain_id_type domain_ID  = input_file("subdomain_ID",1);

	//read in reference mesh
	Mesh ref_mesh(init.comm(), dim);
	ExodusII_IO ref_mesh_reader(ref_mesh);
	ref_mesh_reader.read(ref_mesh_name);
	//ref_mesh.prepare_for_use();
	//ref_mesh.print_info(); 
/*
	//read in mesh for the transfer
	Mesh new_mesh(init.comm(), dim);
	ExodusII_IO new_mesh_reader(new_mesh);
	new_mesh_reader.read(new_mesh_name);
	new_mesh.prepare_for_use();
	new_mesh.print_info();

	//print out element variable names
	const std::vector<std::string> ref_elem_vars = ref_mesh_reader.get_elem_var_names();
	const std::vector<std::string> ref_node_vars = ref_mesh_reader.get_nodal_var_names();
	const std::vector<std::string> new_elem_vars = new_mesh_reader.get_elem_var_names();
	const std::vector<std::string> new_node_vars = new_mesh_reader.get_nodal_var_names();

	for(unsigned int i = 0; i < ref_elem_vars.size(); ++i)
	{
		std::cout << ref_elem_vars[i] << "\n";
	}
	
	//build stuff for integration
	//equation system objest
	EquationSystems equation_system(mesh);
	
	//specicfy volume system
	LinearImplicitSystem& volume_system = equation_system.add_system<LinearImplicitSystem>("ComputeVolume");
	
	//specifiy nodal change systems
	System& sol_x0 = equation_system.add_system<System>("x0");
	System& sol_x1 = equation_system.add_system<System>("x1");
	System& sol_x2 = equation_system.add_system<System>("x2");
	sol_x0.add_variable("X_0", static_cast<Order>(mesh_order));
	sol_x1.add_variable("X_1", static_cast<Order>(mesh_order));
	sol_x2.add_variable("X_2", static_cast<Order>(mesh_order));
	
	// add variables to system, attach assemble function, and initialize system
	volume_system.add_variable ("nerd", static_cast<Order>(mesh_order), LAGRANGE);
	equation_system.init();
	
	//create an fe type
	FEType fe_type = volume_system.variable_type(0);
	UniquePtr<FEBase> fe_elem_face(FEBase::build(dim, fe_type));
	QGauss qface(dim-1, fe_type.default_quadrature_order());
	fe_elem_face->attach_quadrature_rule(&qface);
	
	//integral stuff
	const std::vector<Real>& JxW_face = fe_elem_face->get_JxW();
	const std::vector<libMesh::Point>& qface_normals = fe_elem_face->get_normals();
	const std::vector<Point>& qface_points = fe_elem_face->get_xyz();

	//loop to alter mesh at each time step
	for(unsigned int t = 1; t <= mesh_reader.get_num_time_steps(); ++t)
	{
	
		MeshBase::const_node_iterator no = mesh.nodes_begin();
		const MeshBase::const_node_iterator no_end = mesh.nodes_end();
		mesh_reader.copy_nodal_solution(sol_x0, "X_0", "X_0", t);
		mesh_reader.copy_nodal_solution(sol_x1, "X_1", "X_1", t);
		mesh_reader.copy_nodal_solution(sol_x2, "X_2", "X_2", t);

		for(; no != no_end; ++no)
		{
			const dof_id_type no_id = (*no)->id();
			(*(*no)) = Point(sol_x0.current_solution(no_id),sol_x1.current_solution(no_id),sol_x2.current_solution(no_id));
		}		

		//loop over elements in sideset vectors
		double integral = 0.0;
		for(unsigned int i = 0; i < sideset_faces.size(); ++i)
		{
			fe_elem_face->reinit(sideset_elems[i],sideset_faces[i]);
			for(int qp = 0; qp < qface_points.size(); ++qp)
			{
				integral += -(1.0/3.0) * qface_points[qp] * qface_normals[qp] * JxW_face[qp];
			}
		}
		double top_vol = meter_volume(sort_top);
		double bottom_vol = meter_volume(sort_bottom);
		std::cout << integral + top_vol + bottom_vol << "\n";
	}
*/
	return 0;	
}
