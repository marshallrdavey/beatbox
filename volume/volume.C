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
#include "libmesh/mesh_refinement.h"

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

//make nodeset vectors from nodeset IDs
void populate_nodesets(std::vector<std::vector<Node*>>& node_vec, std::vector<std::vector<int>>& id_vec, Mesh& mesh)
{
	const BoundaryInfo& boundary_info = *mesh.boundary_info;
	MeshBase::const_node_iterator no = mesh.nodes_begin();
	const MeshBase::const_node_iterator no_end = mesh.nodes_end();
	for( ; no != no_end; ++no)
	{
		std::vector<boundary_id_type> bdry_ids;
		boundary_info.boundary_ids(*no, bdry_ids);
		if(find_first_of(bdry_ids.begin(), bdry_ids.end(), id_vec[0].begin(), id_vec[0].end()) != bdry_ids.end())
		{
			node_vec[0].push_back(*no);
            // std::cout << (*no)->get_info() << "\n";
		}
		if(find_first_of(bdry_ids.begin(), bdry_ids.end(), id_vec[1].begin(), id_vec[1].end()) != bdry_ids.end())
		{
			node_vec[1].push_back(*no);
            // std::cout << (*no)->get_info() << "\n";
		}
	}
	// std::cout << "nodeset " << id_vec[0][0] << " contains " << node_vec[0].size() << " nodes\n";
	// std::cout << "nodeset " << id_vec[1][0] << " contains " << node_vec[1].size() << " nodes\n\n";
}

//compute the mean distance between two nodesets
void nodeset_distance(double& dist, std::vector<std::vector<Node*>>& node_vec)
{
	//double min = 100000.0; //add in to see minimum distance
	for(int i = 0; i < node_vec[0].size(); ++i)
	{
		for(int j = 0; j < node_vec[1].size(); ++j)
		{
			const double d = (*node_vec[0][i] - *node_vec[1][j]).norm();
			dist += d;
			/*if(d < min)
			{
				min = d;
			}*/
		}
	}
	dist = dist/node_vec[0].size()/node_vec[1].size();
	//std::cout << "the minimum distance is " << min << "\n";
}

// compute the max distance between two nodes on a nodeset
// this function is not exact
void diameter(std::vector<std::vector<Node*>>& node_vec)
{
	double tmax = 0.0;
	for(int i = 1; i < node_vec[0].size(); ++i)
	{
		const double d = (*node_vec[0][i] - *node_vec[0][0]).norm();
		if(d > tmax)
		{
			tmax = d;
		} 
	}

    double bmax = 0.0;
    for(int i = 1; i < node_vec[1].size(); ++i)
    {
        const double d = (*node_vec[1][i] - *node_vec[1][0]).norm();
        if(d > bmax)
        {
            bmax = d;
        }
    }
	std::cout << "\nthe top diameter is " << tmax << "\n";
    std::cout << "the bottom diameter is " << bmax << "\n";
}

//populate a vector of element pointers and vector of node pointers for the face on the subset
void populate_sidesets(std::vector<Elem*>& elem_vec, std::vector<int>& face_vec, std::vector<std::vector<Node*>>& side_vec, std::vector<int>& id_vec, Mesh& mesh)
{
	const BoundaryInfo& boundary_info = *mesh.boundary_info;
	MeshBase::const_element_iterator el = mesh.elements_begin();
	const MeshBase::const_element_iterator el_end = mesh.elements_end();
	for(; el != el_end; ++el)
	{
		for(unsigned int side = 0; side < (*el)->n_sides(); ++side)
		{
			std::vector<boundary_id_type> bdry_ids; 
			boundary_info.boundary_ids(*el, side, bdry_ids);
			if(find_first_of(bdry_ids.begin(), bdry_ids.end(), id_vec.begin(), id_vec.end()) != bdry_ids.end())
			{
				elem_vec.push_back(*el);
				face_vec.push_back(side);
				const std::vector<unsigned int> points = (*el)->nodes_on_side(side);
				std::vector<Node*> node_vec;
				for(unsigned node = 0; node < points.size(); ++node)
				{	
					node_vec.push_back((*el)->node_ptr(points[node]));
				}
				side_vec.push_back(node_vec);
			}
		}
	}
	std::cout << "there are " << elem_vec.size() << " sides in sideset " << id_vec[0] << "\n";
}

//main function compute sideset volume, v2
void sideset_volume(double& vol, std::vector<Elem*>& elem_vec, std::vector<std::vector<Node*>>& side_vec, Mesh& mesh)
{
	for(int i = 0; i < elem_vec.size(); ++i)
	{
		const Elem* el = elem_vec[i];
		const std::vector<Node*> side = side_vec[i];
		const Point center = el->centroid();
		Point com;
        	for(int j=0; j<3; ++j)
        	{
                	com += *side[j]/3.0;
        	}
        	Point n_vec = (*side[1] - *side[0]).cross(*side[2] - *side[0]).unit();
        	if(n_vec*(com - center) < 0)
        	{
                	n_vec = -1.0*n_vec;
        	}
        	double signed_volume = fabs(side[0]->cross(*side[1])*(*side[2]))/6.0;
       		if(n_vec*com > 0)
        	{
                	signed_volume = -1.0*signed_volume;
                	//std::cout << "flipped\n"; //uncomment to see when the volume is subtracted
        	}
		else if(n_vec*com == 0)
	        {
                	signed_volume = 0.0;
        	}
		vol += signed_volume;
	}
	std::cout << "the side set volume is " << vol << "\n";
}

//sort nodeset vector to make a cirle
void nodeset_sort(std::vector<Node*>& out_vec, std::vector<Node*>& node_vec)
{
	out_vec.clear();
	out_vec.push_back(node_vec[0]);
	std::vector<Node*> temp_vec;
	temp_vec.assign(node_vec.begin()+1,node_vec.end());
	std::vector<Node*> n_vec;
	while(out_vec.size() < node_vec.size())
	{
		Node* n = temp_vec[0];
		for(int i = 1; i < temp_vec.size(); ++i)
		{
			if((*out_vec.back() - *temp_vec[i]).norm() < (*out_vec.back()-*n).norm())
			{
				n_vec.push_back(n);
				n = temp_vec[i];
			}
			else
			{
				n_vec.push_back(temp_vec[i]);
			}
		}
		out_vec.push_back(n);
		temp_vec.clear();
		temp_vec.assign(n_vec.begin(), n_vec.end());
		n_vec.clear();
	}
}

//calculate volume from the sorted node vector 
double meter_volume(std::vector<Node*>& node_vec)
{
	Point com;
	const int size = node_vec.size();
	for(int i=0; i<size; ++i)
	{	
		com += *node_vec[i]/size;
	}
	double vol = fabs(((*node_vec.back()).cross(*node_vec[0]))*com)/6.0;
	for(int i = 0; i < size - 1; ++i)
	{
		vol += fabs((*node_vec[i]).cross(*node_vec[i+1])*com)/6.0;
	}
    // std::cout << "meter volume: " << vol << "\n";
	return vol;
}

//main function
int main (int argc, char** argv)
{
	//initialize libmesh
	LibMeshInit init(argc, argv);

	//load the input file
	GetPot input_file(argv[1]);

	// read in parameters from the input file
	const unsigned int dim              	    = input_file("dimension",3);
	const std::string mesh_name      	        = input_file("mesh_name","");
	const unsigned int mesh_order               = input_file("mesh_order",1);
	const std::string inside_ID                 = input_file("inside","");
	const std::string top_ID		            = input_file("top","");
	const std::string bottom_ID		            = input_file("bottom","");
	const libMesh::subdomain_id_type domain_ID  = input_file("subdomain_ID",1);

	// create a simple FE mesh.
	Mesh mesh(init.comm(), dim);
	ExodusII_IO mesh_reader(mesh);
	mesh_reader.read(mesh_name);
	mesh.prepare_for_use();
	mesh.print_info();  //print out general information about the mesh

	// load boundary
	const BoundaryInfo& boundary_info = *mesh.boundary_info;
/*
	//print out sideset IDs
	std::cout << "\n" << "sideset IDs and names are... \n";
	std::set<short int>::iterator ii;
	for(ii = boundary_info.get_side_boundary_ids().begin();  ii != boundary_info.get_side_boundary_ids().end(); ++ii)
	{
		std::cout << *ii <<" " << boundary_info.get_sideset_name(*ii) << std::endl;
	}
	std::cout << "\n \n";

	//print out nodeset IDs
	std::cout << "nodeset IDs and names are... \n";
	for(ii = boundary_info.get_node_boundary_ids().begin(); ii != boundary_info.get_node_boundary_ids().end(); ++ii)
	{
		std::cout << *ii <<" " << boundary_info.get_nodeset_name(*ii) << std::endl;
	}
    std::cout << "\n\n";	
*/
	// make node pointer vector
	std::vector<std::vector<Node*>> nodeset_nodes(2);
	std::vector<std::vector<int>> nodesets(2);
	parse_ID_string(nodesets[0], top_ID);
	parse_ID_string(nodesets[1], bottom_ID);
	populate_nodesets(nodeset_nodes, nodesets, mesh);

	/*
	//mean distance between nodesets
	double dist;
	nodeset_distance(dist, nodeset_nodes);
	std::cout << "the mean distance between nodesets is " << dist << "\n";
	diameter(nodeset_nodes);	//run the diameter function
    */	

	// make element pointer vector and corresponding side vector
	std::vector<Elem*> sideset_elems;
	std::vector<int> sideset_faces;    
	std::vector<std::vector<Node*>> sideset_nodes;
	std::vector<int> sidesets;
	parse_ID_string(sidesets, inside_ID);
	populate_sidesets(sideset_elems, sideset_faces, sideset_nodes, sidesets, mesh);

	/*
	//volume of solid without nodesets
	double vol;
	sideset_volume(vol, sideset_elems, sideset_nodes, mesh);
	*/

	// sorted node pointer list
	std::vector<Node*> sort_top;
	nodeset_sort(sort_top, nodeset_nodes[0]);
	std::vector<Node*> sort_bottom;
	nodeset_sort(sort_bottom, nodeset_nodes[1]);

	// build stuff for integration
	// equation system object
	EquationSystems equation_system(mesh);
	
	// specicfy volume system
	ExplicitSystem& sol_system = equation_system.add_system<ExplicitSystem>("position");
	const unsigned int x0_var = sol_system.add_variable("X_0", static_cast<Order>(mesh_order), LAGRANGE);
	const unsigned int x1_var = sol_system.add_variable("X_1", static_cast<Order>(mesh_order), LAGRANGE);
	const unsigned int x2_var = sol_system.add_variable("X_2", static_cast<Order>(mesh_order), LAGRANGE);
    const unsigned int J_var = sol_system.add_variable("J", static_cast<Order>(mesh_order), LAGRANGE);

	// initialize system
	equation_system.init();
	
	// create an fe type
    const DofMap& dof_map = sol_system.get_dof_map();
	FEType fe_type = dof_map.variable_type(x0_var);
	UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_type));
	QGauss qface(dim-1, fe_type.default_quadrature_order());
	fe_face->attach_quadrature_rule(&qface);
	
	// integral stuff
	const std::vector<Real>& JxW_face = fe_face->get_JxW();
	const std::vector<libMesh::Point>& qface_normals = fe_face->get_normals();
	const std::vector<Point>& qface_points = fe_face->get_xyz();

    // find number of mesh time steps 
    const int n_t_steps = mesh_reader.get_num_time_steps();
/*    
    // no time steps
    if(n_t_steps < 1)
    {
        // loop over elements in sideset vectors
		double volume = 0.0;
		for(unsigned int i = 0; i < sideset_faces.size(); ++i)
		{
		    fe_face->reinit(sideset_elems[i],sideset_faces[i]);
			for(int qp = 0; qp < qface_points.size(); ++qp)
			{
			    volume += -(1.0/3.0) * qface_points[qp] * qface_normals[qp] * JxW_face[qp];
			}
		}
        // meter volume call
		double top_vol = meter_volume(sort_top);
		double bottom_vol = meter_volume(sort_bottom);
        
        // output
		std::cout << "this mesh has no time steps\n"
                  << "volume:" << volume + top_vol + bottom_vol << "\n\n";
    } 
*/    
    // at least one time step
    unsigned int bad = 0;
    if(n_t_steps > 0)
    { 
        // vector of times
        const std::vector<Real> t_steps = mesh_reader.get_time_steps();
	    std::vector<dof_id_type> dof_indices;

        // loop to alter mesh at each time step
	    for(unsigned int t = 1; t <= n_t_steps; ++t)
	    {

		    MeshBase::const_node_iterator no = mesh.nodes_begin();
		    const MeshBase::const_node_iterator no_end = mesh.nodes_end();
		    mesh_reader.copy_nodal_solution(sol_system, "X_0", "X_0", t);
		    mesh_reader.copy_nodal_solution(sol_system, "X_1", "X_1", t);
		    mesh_reader.copy_nodal_solution(sol_system, "X_2", "X_2", t);

		    for(; no != no_end; ++no)
		    {
			    dof_map.dof_indices(*no, dof_indices);
                const Number diff =
			    ((*(*no)) - Point(sol_system.current_solution(dof_indices[x0_var]),
                                 sol_system.current_solution(dof_indices[x1_var]),
                                 sol_system.current_solution(dof_indices[x2_var]))).norm();
                if(diff > 0.1) ++bad;
		    }		

		    // loop over elements in sideset vectors
		    double volume = 0.0;
		    for(unsigned int i = 0; i < sideset_faces.size(); ++i)
		    {
			    fe_face->reinit(sideset_elems[i],sideset_faces[i]);
			    for(int qp = 0; qp < qface_points.size(); ++qp)
			    {
				    volume += -(1.0/3.0) * qface_points[qp] * qface_normals[qp] * JxW_face[qp];
			    }
		    }
            // meter volume call
		    double top_vol = meter_volume(sort_top);
		    double bottom_vol = meter_volume(sort_bottom);
       
            // output
		    std::cout << "time step:" << t << "\n"
                      << "time:" << t_steps[t-1] << "\n" 
                      << "main volume: " << volume << "\n"
                      << "top volume: " << top_vol << "\n"
                      << "bottom volume: " << bottom_vol << "\n"
                      << "total volume:" << volume + top_vol + bottom_vol << "\n";

            if(t >= 0)
            {
                double min_J = 20;
                mesh_reader.copy_nodal_solution(sol_system, "J", "J", t);
                std::cout << "J copied\n";
/*                no = mesh.nodes_begin();
                for(; no != no_end; ++no)
                {
                    if(sol_J.current_solution((*no)->id()) < min_J)
                    {
                        min_J = sol_J.current_solution((*no)->id());
                    }
                }
                std::cout << "min_J:" << min_J << "\n\n";*/
            } 
	    }
    }
    std::cout << "bad nodes: " << bad << "\n";
	/*the addition or subtraction of the nodeset volumes depends on the location of the origin
	see paper: http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
	for a better understanding of this method*/

	return 0;	
}
