// Irregular packing example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double x_min=0,x_max=50;
const double y_min=0,y_max=50;
const double z_min=0,z_max=50;

// Golden ratio constants
const double Phi=0.5*(1+sqrt(5.0));
const double phi=0.5*(1-sqrt(5.0));

// Set up the number of blocks that the container is divided
// into.
const int n_x=5,n_y=5,n_z=5;

// Set up ideal particle size
const double vdwr=1.0;

// Create a wall class that, whenever called, will replace the Voronoi cell
// with a prescribed shape, in this case a dodecahedron
class wall_initial_shape : public wall {
	public:
		wall_initial_shape() {

			// Create a dodecahedron
			v.init(-vdwr,vdwr,-vdwr,vdwr,-vdwr,vdwr);
			v.plane(0,Phi,1);v.plane(0,-Phi,1);v.plane(0,Phi,-1);
			v.plane(0,-Phi,-1);v.plane(1,0,Phi);v.plane(-1,0,Phi);
			v.plane(1,0,-Phi);v.plane(-1,0,-Phi);v.plane(Phi,1,0);
			v.plane(-Phi,1,0);v.plane(Phi,-1,0);v.plane(-Phi,-1,0);
		};
		bool point_inside(double x,double y,double z) {return true;}
		bool cut_cell(voronoicell &c,double x,double y,double z) {

			// Set the cell to be equal to the dodecahedron
			c=v;
			return true;
		}
		bool cut_cell(voronoicell_neighbor &c,double x,double y,double z) {

			// Set the cell to be equal to the dodecahedron
			c=v;
			return true;
		}
	private:
		voronoicell v;
};

int main() {

	// Create a container with the geometry given above. This is bigger
	// than the particle packing itself.
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Create the "initial shape" wall class and add it to the container
	wall_initial_shape(wis);
	con.add_wall(wis);

	// Import the irregular particle packing
	con.import("test.xyz");

	// Save the particles and Voronoi cells in POV-Ray format
        con.draw_cells_gnuplot("test.gnu");
        
	// Do a custom output routine that outputs the particle IDs, volume and neighbors
	con.print_custom("%i %v %s %f %n","test.vol");
}
