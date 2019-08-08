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
		}
		if(find_first_of(bdry_ids.begin(), bdry_ids.end(), id_vec[1].begin(), id_vec[1].end()) != bdry_ids.end())
		{
			node_vec[1].push_back(*no);
		}
	}
	std::cout << "\nnodeset " << id_vec[0][0] << " contains " << node_vec[0].size() << " nodes\n";
	std::cout << "nodeset " << id_vec[1][0] << " contains " << node_vec[1].size() << " nodes\n";
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

//compute the approximate diameter of the nodesets
/*void diameter(std::vector<std::vector<Node*>>& node_vec)
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
	std::cout << "the top diameter is " << tmax << "\n";
        std::cout << "the bottom diameter is " << bmax << "\n";
}*/

/*//populate an element vecotr and face vector from a sideset ID, v1
void populate_sidesets(std::vector<Elem*>& elem_vec, std::vector<int>& side_vec, std::vector<int>& id_vec, Mesh& mesh)
{
	const BoundaryInfo& boundary_info = *mesh.boundary_info;
	MeshBase::const_element_iterator el = mesh.elements_begin();
	const MeshBase::const_element_iterator el_end = mesh.elements_end();
	for(; el != el_end; ++el)
	{
		for(int i = 0; i < (*el)->n_sides(); ++i)
		{
			std::vector<boundary_id_type> bdry_ids; 
			boundary_info.boundary_ids(*el, i, bdry_ids);
			if(find_first_of(bdry_ids.begin(), bdry_ids.end(), id_vec.begin(), id_vec.end()) != bdry_ids.end())
			{
				elem_vec.push_back(*el);
				side_vec.push_back(i);
			}
		}
	}
	std::cout << "there are " << side_vec.size() << " sides in sideset " << id_vec[0] << "\n";
}*/

//populate an element vecotr and face vector from a sideset ID, v2
void populate_sidesets(std::vector<Elem*>& elem_vec, std::vector<std::vector<Node*>>& side_vec, std::vector<int>& id_vec, Mesh& mesh)
{
	const BoundaryInfo& boundary_info = *mesh.boundary_info;
	MeshBase::const_element_iterator el = mesh.elements_begin();
	const MeshBase::const_element_iterator el_end = mesh.elements_end();
	for(; el != el_end; ++el)
	{
		for(unsigned int side = 0; side < (*el)->n_sides(); ++i)
		{
			std::vector<boundary_id_type> bdry_ids; 
			boundary_info.boundary_ids(*el, side, bdry_ids);
			if(find_first_of(bdry_ids.begin(), bdry_ids.end(), id_vec.begin(), id_vec.end()) != bdry_ids.end())
			{
				elem_vec.push_back(*el);
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

/*//compute the signed volume of a sideset face, v1
double signed_volume(std::vector<Point>& p_vec, const Point center)
{
	Point com;
	const int size = p_vec.size();
	for(int i=0; i<size; ++i)
	{	
		com += p_vec[i]/size;
	}
	Point n_vec = (p_vec[1] - p_vec[0]).cross(p_vec[2] - p_vec[0]).unit();
	if(n_vec*(com - center) < 0)
	{
		n_vec = -1.0*n_vec;
	}
	double vol = fabs(p_vec[0].cross(p_vec[1])*p_vec[2])/6.0;
	if(n_vec*com > 0)
	{
		vol = -1.0*vol;
		//std::cout << "flipped\n"; //uncomment to see when the volume is subtracted
	}
	else if(n_vec*com == 0)
	{
		vol = 0.0;
	}
	return vol;
}*/

//compute the signed volume of a sideset face, v2
double signed_volume(const std::vector<Node*> p_vec, const Point center)
{
	Point com;
	const int size = p_vec.size();
	for(int i=0; i<size; ++i)
	{	
		com += *p_vec[i]/size;
	}
	Point n_vec = (*p_vec[1] - *p_vec[0]).cross(*p_vec[2] - *p_vec[0]).unit();
	if(n_vec*(com - center) < 0)
	{
		n_vec = -1.0*n_vec;
	}
	double vol = fabs(*p_vec[0].cross(*p_vec[1])*(*p_vec[2]))/6.0;
	if(n_vec*com > 0)
	{
		vol = -1.0*vol;
		//std::cout << "flipped\n"; //uncomment to see when the volume is subtracted
	}
	else if(n_vec*com == 0)
	{
		vol = 0.0;
	}
	return vol;
}

/*//main function compute sideset volume, v1
void sideset_volume(double& vol, std::vector<Elem*>& elem_vec, std::vector<int>& side_vec, Mesh& mesh)
{
	for(int i = 0; i < elem_vec.size(); ++i)
	{
		const Elem* el = elem_vec[i];
		const int side = side_vec[i];
		std::vector<Point> nodes;
		for(int j = 0; j < el->n_nodes(); ++j)
		{
			if(el->is_node_on_side(j, side))
			{
				nodes.push_back(el->point(j));
			}
		}
		const Point centroid = el->centroid();
		vol += signed_volume(nodes, centroid);
	}
	std::cout << "the side set volume is " << vol << "\n";
}*/

//main function compute sideset volume, v2
void sideset_volume(double& vol, std::vector<Elem*>& elem_vec, std::vector<std::vector<Node*>>& side_vec, Mesh& mesh)
{
	for(int i = 0; i < elem_vec.size(); ++i)
	{
		const Elem* el = elem_vec[i];
		const std::vector<Node*> side = side_vec[i];
		const Point centroid = el->centroid();
		vol += signed_volume(side, centroid);
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
	return vol;
}

//main function
int main (int argc, char** argv)
{
	//initialize libmesh
	LibMeshInit init(argc, argv);

	//load the input file
	GetPot input_file(argv[1]);

	//read in parameters from the input file
	const unsigned int dim              	    = input_file("dimension", 3);
	const std::string mesh_name      	    = input_file("mesh_name","");
	const std::string inside_ID		    = input_file("inside","");
	const std::string top_ID		    = input_file("top","");
	const std::string bottom_ID		    = input_file("bottom","");
	const libMesh::subdomain_id_type domain_ID  = input_file("subdomain_ID",1);

	//create a simple FE mesh.
	Mesh mesh(init.comm(), dim);
	ExodusII_IO mesh_reader(mesh);
	mesh_reader.read(mesh_name);
	mesh.prepare_for_use();
	//mesh.print_info();  //print out general information about the mesh

	//load boundary
	const BoundaryInfo& boundary_info = *mesh.boundary_info;

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

	//make node pointer vector
	std::vector<std::vector<Node*>> nodeset_nodes(2);
	std::vector<std::vector<int>> nodesets(2);
	parse_ID_string(nodesets[0], top_ID);
	parse_ID_string(nodesets[1], bottom_ID);
	populate_nodesets(nodeset_nodes, nodesets, mesh);

	//mean distance between nodesets
	double dist;
	nodeset_distance(dist, nodeset_nodes);
	std::cout << "the mean distance between nodesets is " << dist << "\n";
	//diameter(nodeset_nodes);	//run the diameter function

	//make element pointer vector and corresponding side vector
	std::vector<Elem*> sideset_elems;
	std::vector<std::vector<Node*>> sideset_faces;
	std::vector<int> sidesets;
	parse_ID_string(sidesets, inside_ID);
	populate_sidesets(sideset_elems, sideset_faces, sidesets, mesh);

	//volume of solid without nodesets
	/*double vol;
	sideset_volume(vol, sideset_elems, sideset_faces, mesh);
	
	//sorted node pointer list
	std::vector<Node*> sort_top;
	nodeset_sort(sort_top, nodeset_nodes[0]);
	std::vector<Node*> sort_bottom;
	nodeset_sort(sort_bottom, nodeset_nodes[1]);

	//nodeset volume
	double top_vol = meter_volume(sort_top);
	double bottom_vol = meter_volume(sort_bottom);
	std::cout << "top volume = " << top_vol << "\nbottom volume = " << bottom_vol << std::endl;
	/*the addition or subtraction of the nodeset volumes depends on the location of the origin
	see paper: http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
	for a better understanding of this method*/
	//std::cout << "total volume = " << vol + top_vol - bottom_vol << "\n";

	return 0;
}
