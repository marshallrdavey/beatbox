// adapted from libMesh example, miscellaneous/ex5
#include <iostream>
#include <string>
#include "math.h"

// LibMesh include files.
#include "libmesh/exodusII_io.h"
#include "libmesh/exodusII_io_helper.h"
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
#include "libmesh/fe_interface.h"
#include "libmesh/getpot.h"
#include "libmesh/error_vector.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/point_locator_tree.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/mesh_function.h"
#include "libmesh/exact_solution.h"
#include "libmesh/system.h"
#include "libmesh/explicit_system.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

//main function
int main (int argc, char** argv)
{
	//initialize libmesh
	LibMeshInit init(argc, argv);

	//load the input file
	GetPot input_file(argv[1]);

	//read in parameters from the input file
	const unsigned int dim              	 = input_file("dimension",3);
	const std::string ref_mesh_name      	 = input_file("ref_mesh_name","");
	const std::string new_mesh_name          = input_file("new_mesh_name","");

	//read in reference mesh
	Mesh ref_mesh(init.comm(), dim);
	ExodusII_IO ref_mesh_reader(ref_mesh);
	ref_mesh_reader.read(ref_mesh_name);
	//ref_mesh.allow_renumbering(false);
	ref_mesh.prepare_for_use();
	//ref_mesh.print_info();
/*
	//read in mesh to add fibers to 
	Mesh new_mesh(init.comm(), dim);
	ExodusII_IO new_mesh_reader(new_mesh);
	new_mesh_reader.read(new_mesh_name);
	//new_mesh.allow_renumbering(false);
	new_mesh.prepare_for_use();
	//new_mesh.print_info();
*/	
	std::cout << "elemental data\n";
	//get element variable names from reference mesh 
	const std::vector<std::string> ref_elem_vars = ref_mesh_reader.get_elem_var_names();

	for(unsigned int i = 0; i < ref_elem_vars.size(); ++i)
	{
		std::cout << ref_elem_vars[i] << "\n";
	}
/*	
	std::cout << "nodal data\n";
	//get node variable names from reference mesh 
	const std::vector<std::string> ref_node_vars = ref_mesh_reader.get_nodal_var_names();

	for(unsigned int i = 0; i < ref_node_vars.size(); ++i)
	{
		std::cout << ref_node_vars[i] << "\n";
	}
*/
    EquationSystems system(ref_mesh);
    EquationSystems ref_system(ref_mesh);

    System& fibers = system.add_system<System>("fiber");
    //System& ref_fibers = ref_system.add_system<System>("fiber_info");
    System& gogo = ref_system.add_system<System>("gogo");
    System& ref_fibers = ref_system.add_system<System>("fiber_info");

    gogo.add_variable(ref_elem_vars[3], FEType(CONSTANT, MONOMIAL));
    for(unsigned int i = 0; i < 3 /*ref_elem_vars.size()*/; ++i)
    {
        ref_fibers.add_variable(ref_elem_vars[i], FEType(CONSTANT, MONOMIAL));  
        fibers.add_variable(ref_elem_vars[i], FEType(CONSTANT, MONOMIAL));
    }
    //gogo.add_variable(ref_elem_vars[3], FEType(CONSTANT, MONOMIAL));
    
    ref_system.init();
//    system.init();
    
//    std::cout << system.get_info() << "\n"; 
    std::cout << ref_system.get_info() << "\n";

    //copy fiber solutions to ref_system
    for(unsigned int i = 0; i < 3; ++i)
    {
        ref_mesh_reader.copy_elemental_solution(ref_fibers,ref_elem_vars[i],ref_elem_vars[i],1);
    }
/*
    //initiate element iterator to transfer system data from ref_fibers to new_fibers
    MeshBase::const_element_iterator el = ref_mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = ref_mesh.active_local_elements_end();
    DofMap& dof_map = ref_fibers.get_dof_map();
    std::vector<dof_id_type> dof_indices;

    for(;el != el_end; ++el)
    {
        const Elem* elem = *el;
        dof_map.dof_indices(elem, dof_indices);

        for(unsigned int i = 0; i < 3; ++i)
        {
            const dof_id_type id = dof_indices[i];
            ref_fibers.solution->set(id,5.0);
            std::cout << (*ref_fibers.solution)(id) << "\n";
        }

    }
*/

/*
	//equation system objects
	EquationSystems ref_system(ref_mesh);
	EquationSystems new_system(new_mesh);

	//specifiy fiber systems
	System& ref_fibers = ref_system.add_system<System>("fiber_info");
	System& new_fibers = new_system.add_system<System>("fiber_info");

	for(unsigned int i = 0; i < ref_elem_vars.size(); ++i)
	{
		ref_fibers.add_variable(ref_elem_vars[i], FEType(CONSTANT, MONOMIAL));
		new_fibers.add_variable(ref_elem_vars[i], FEType(CONSTANT, MONOMIAL));
	}
	ref_system.init();
	new_system.init();
	//std::cout << ref_system.get_info() <<"\n";
	//std::cout << new_system.get_info() <<"\n";

	//copy fiber solutions to ref_system
	for(unsigned int i = 0; i < ref_elem_vars.size(); ++i)
	{
		ref_mesh_reader.copy_elemental_solution(ref_fibers,ref_elem_vars[i],ref_elem_vars[i],1);
	}

	//initiate element iterator to transfer system data from ref_fibers to new_fibers
	MeshBase::const_element_iterator el = ref_mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator el_end = ref_mesh.active_local_elements_end();
	DofMap& dof_map = ref_fibers.get_dof_map();
	std::vector<dof_id_type> dof_indices;
	for(;el != el_end; ++el)
	{
		const Elem* elem = *el;
		dof_map.dof_indices(elem, dof_indices);
		for(unsigned int i = 0; i < ref_elem_vars.size(); ++i)
		{
			const dof_id_type id = dof_indices[i];
        		new_fibers.solution->set(id,(*ref_fibers.solution)(id));

		}
	}

	ExodusII_IO file_writer(new_mesh);
  	file_writer.write("fibers_for_"+new_mesh_name);
  	file_writer.write_element_data(new_system);
	//ExodusII_IO (new_mesh).write_timestep("results_for_"+new_mesh_name,new_system,1,0);
*/
	return 0;	
}
