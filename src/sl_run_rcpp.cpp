#define _USE_MATH_DEFINES

#include <vector>
#include <map>
#include <string> 
#include <cmath>
#include <limits>
#include <algorithm>

#include <Rcpp.h>
using namespace Rcpp;

#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num() 0
	#define omp_get_num_threads() 0
#endif

// [[Rcpp::plugins(openmp)]]


struct relativeCoords {
	int row;
	int col;
};

struct vertex3D {
	double x;
	double y;
	double z;
};

struct interception {

	// Ray intercepting the crown part
	int idRay;

	// Target cell the ray is directed to
	int idTargetCell;

	// Tree that hold the intercepted crown part
	int idTree;

	// Crown part intercepted by the ray
	int idCrownPart;

	// Length of the path across the crown
	double length;

	// Distance between interception point (middle of full path) and target cell
	double distance;

};

class Ray {

private:
	//Id of the ray
	int id;

	// Relative coordinates of cells containing trees potentially crossed by the ray (relatively to an unknown target cell)
	std::vector<relativeCoords> potentialRelCells;

	// Energy of the ray above canopy (in MJ / m2)
	double incidentEnergy;

	// Geometry of the ray and derived variables to save computing time
	double heightAngle;
	double cosHeightAngle;
	double sinHeightAngle;

	double azimuth;
	double cosAzimuth;
	double sinAzimuth;

	// Output interceptions vector
	std::vector<interception*> interceptions;


public:
	Ray(int id, double energy, double height_angle, double azimuth) {
		this->id = id;
		this->incidentEnergy = energy;
		this->heightAngle = height_angle;
		this->cosHeightAngle = cos(height_angle);
		this->sinHeightAngle = sin(height_angle);
		this->azimuth = azimuth;
		this->cosAzimuth = cos(azimuth);
		this->sinAzimuth = sin(azimuth);
	}

	// Getters
	int getId() { return(this->id); }
	double getIncidentEnergy() { return(this->incidentEnergy); }
	double getHeightAngle() { return(this->heightAngle); }
	double getCosHeightAngle() { return(this->cosHeightAngle); }
	double getSinHeightAngle() { return(this->sinHeightAngle); }
	double getAzimuth() { return(this->azimuth); }
	double getCosAzimuth() { return(this->cosAzimuth); }
	double getSinAzimuth() { return(this->sinAzimuth); }

	void addPotentialRelCell(int row, int col) { this->potentialRelCells.push_back(relativeCoords{ row, col }); }
	relativeCoords getPotentialRelCell(int id) { return(this->potentialRelCells[id]); }
	int getNPotCells() { return(this->potentialRelCells.size()); }


private:
	//void computeInterceptionCrownParaboloid(Cell* target_cell, Tree* tree, shiftedCell& shifted_cell) {

	//	// GET POSITION OF THE PARABOLOID
	//	double x = tree->getX();
	//	double y = tree->getY();
	//	double z = tree->getZ() + tree->getCrownBaseHeight();

	//	// SHIFT POSITION OF THE PARABOLOID
	//	// 2 shifts: the first one to set the target cell as origine, and the second to acoount for outside (*tree) when torus system
	//	x = x + shifted_cell.xShift - target_cell->getX();
	//	y = y + shifted_cell.yShift - target_cell->getY();
	//	z = z + shifted_cell.zShift - target_cell->getZ();

	//	// GET PARAMETERS OF THE PARABOLOID
	//	double a = tree->getCrownRadius();
	//	double b = tree->getCrownRadius();
	//	double h = tree->getHeight() - tree->getCrownBaseHeight();


	//	// FIND SOLUTION OF THE QUADRATIC EQUATION (A * x * x + B * x + C = 0)
	//	// Equation giving distance between ray intersection with (*tree) crown and target cell center

	//	   // Compute a, b and c coefficients
	//	double coef_a = this->cosHeightAngle * this->cosHeightAngle * this->cosAzimuth * this->cosAzimuth / (a * a) +
	//		this->cosHeightAngle * this->cosHeightAngle * this->sinAzimuth * this->sinAzimuth / (b * b);

	//	double coef_b = -2 * x * this->cosHeightAngle * this->cosAzimuth / (a * a) -
	//		2 * y * this->sinAzimuth * this->cosHeightAngle / (b * b) +
	//		this->sinHeightAngle / h;

	//	double coef_c = (x * x) / (a * a) + (y * y) / (b * b) - (z + h) / h;


	//	// Find positive solutions of the above quadratic equations (distance to target cell)
	//	// if not 2 solutions, set to NaN
	//	// Because negative = no interception and null = 1 interception = do not consider tangent rays with crown
	//	double delta = coef_b * coef_b - 4 * coef_a * coef_c;

	//	double sol1 = std::numeric_limits<double>::quiet_NaN();
	//	double sol2 = std::numeric_limits<double>::quiet_NaN();

	//	if (delta > 0) {
	//		double delta_sqrt = sqrt(delta);
	//		sol1 = (-coef_b + delta_sqrt) / (2 * coef_a);
	//		sol2 = (-coef_b - delta_sqrt) / (2 * coef_a);
	//	}



	//	// COMPUTE COORDINATES OF INTERCEPTION POINTS WITH CROWN LIMITS AND TARGET POINT

	//		// interception coords of the two roots
	//	double v1_x = sol1 * this->cosHeightAngle * this->cosAzimuth;
	//	double v1_y = sol1 * this->cosHeightAngle * this->sinAzimuth;
	//	double v1_z = sol1 * this->sinHeightAngle;

	//	double v2_x = sol2 * this->cosHeightAngle * this->cosAzimuth;
	//	double v2_y = sol2 * this->cosHeightAngle * this->sinAzimuth;
	//	double v2_z = sol2 * this->sinHeightAngle;

	//	// interception with the base plane
	//	double solb = z / this->sinHeightAngle;

	//	double vb_x = solb * this->cosHeightAngle * this->cosAzimuth;
	//	double vb_y = solb * this->cosHeightAngle * this->sinAzimuth;
	//	double vb_z = z; // solb * this->sinHeightAngle



	//	// FIND SOLUTIONS

	//		// Get the limits of the box between shifted paraboloid center and crown limit
	//	double x_bbox_max = x + a;
	//	double y_bbox_max = y + b;
	//	double z_bbox_max = z + h;

	//	double x_bbox_min = x - a;
	//	double y_bbox_min = y - b;
	//	double z_bbox_min = z;

	//	// Find interception points within the crown (edge of ellipsoid or z axis if half ellipsoid = crown base plane)

	//	bool is_v1 =
	//		// interception 2 with crown is above the ground AND
	//		v1_z >= 0 &&
	//		v1_x >= x_bbox_min && v1_x <= x_bbox_max &&
	//		v1_y >= y_bbox_min && v1_y <= y_bbox_max &&
	//		// interception 2 with crown is in crown bbox AND
	//		v1_z >= z_bbox_min && v1_z <= z_bbox_max;


	//	bool is_v2 =
	//		// interception 2 with crown is above the ground AND
	//		v2_z >= 0 &&
	//		// interception 2 with crown is in crown bbox AND
	//		v2_x >= x_bbox_min && v2_x <= x_bbox_max &&
	//		v2_y >= y_bbox_min && v2_y <= y_bbox_max &&
	//		v2_z >= z_bbox_min && v2_z <= z_bbox_max;


	//	bool is_vb =
	//		// interception with crown base is above the ground AND
	//		vb_z >= 0 &&
	//		// interception with crown base is in crown bbox AND
	//		vb_x >= x_bbox_min && vb_x <= x_bbox_max &&
	//		vb_y >= y_bbox_min && vb_y <= y_bbox_max &&
	//		vb_z >= z_bbox_min && vb_z <= z_bbox_max &&
	//		// interception with crown base is in the paraboloid
	//		((vb_x - x) * (vb_x - x) / (a * a) +
	//			(vb_y - y) * (vb_y - y) / (b * b) -
	//			(vb_z - z) / (h * h)) <= 1;


	//	bool is_target =
	//		// interception with target point is in crown bbox AND
	//		0 >= x_bbox_min && 0 <= x_bbox_max &&
	//		0 >= y_bbox_min && 0 <= y_bbox_max &&
	//		0 >= z_bbox_min && 0 <= z_bbox_max &&
	//		// interception with target point is in the paraboloid
	//		(x * x / (a * a) +
	//			y * y / (b * b) -
	//			-z / (h * h)) <= 1;


	//	// FIND PATH LENGTH AND DISTANCE TO CELL CENTER
	//	// Distance between target cell and middle point between two interceptions
	//	// Length is the length of the ray path across the crown
	//	double distance = 0.0;
	//	double length = 0.0;

	//	// If 0 solution, no interception
	//	if (!(!std::isnan(sol1) && is_v1) && !(!std::isnan(sol2) && is_v2) && !is_vb && !is_target) {
	//		return;
	//	}

	//	std::stack<double> distance_tmp;
	//	std::stack<double> length_tmp;

	//	int nsols = 0;
	//	if (!std::isnan(sol1) && is_v1) {
	//		nsols++;
	//		distance_tmp.push(sol1);
	//		length_tmp.push(sol1);
	//	}

	//	if (!std::isnan(sol2) && is_v2) {
	//		nsols++;
	//		distance_tmp.push(sol2);
	//		length_tmp.push(sol2);
	//	}

	//	if (is_vb) {
	//		nsols++;
	//		distance_tmp.push(solb);
	//		length_tmp.push(solb);
	//	}

	//	if (is_target) {
	//		nsols++;
	//		distance_tmp.push(0.0);
	//		length_tmp.push(0.0);
	//	}


	//	// Compute length and distance if 2 solutions
	//	if (nsols == 2) {

	//		// Length
	//		double l1 = length_tmp.top();
	//		length_tmp.pop();
	//		double l2 = length_tmp.top();

	//		length = abs(l1 - l2);

	//		// Distance
	//		double d1 = distance_tmp.top();
	//		distance_tmp.pop();
	//		double d2 = distance_tmp.top();

	//		distance = (d1 + d2) / 2.0;
	//	}
	//	// Otherwise, there is a problem
	//	else {
	//		std::cout << "Not 0 or 2 solutions" << std::endl;
	//	}

	//	// Create and push interception to the vector of ray interceptions
	//	#ifdef _OPENMP
	//		#pragma omp critical
	//	#endif
	//	{
	//		this->interceptions.emplace_back(interception(target_cell, tree, length, distance, shifted_cell.outside));
	//	}
	//}

	//void computeInterceptionCrownEllipsoid(Cell* target_cell, Tree* tree, shiftedCell& shifted_cell) {

	//	// GET POSITION OF THE ELLIPSOID
	//	double x = tree->getX();
	//	double y = tree->getY();
	//	double z = tree->getZ() + tree->getCrownBaseHeight() + (tree->getHeight() - tree->getCrownBaseHeight()) / 2.0;

	//	// SHIFT POSITION OF THE ELLIPSOID
	//	// 2 shifts: the first one to set the target cell as origine, and the second to acoount for outside (*tree) when torus system
	//	x = x + shifted_cell.xShift - target_cell->getX();
	//	y = y + shifted_cell.yShift - target_cell->getY();
	//	z = z + shifted_cell.zShift - target_cell->getZ();

	//	// GET PARAMETERS OF THE ELLIPSOID
	//	double a = tree->getCrownRadius();
	//	double b = tree->getCrownRadius();
	//	double c = (tree->getHeight() - tree->getCrownBaseHeight()) / 2.0;


	//	// FIND SOLUTION OF THE QUADRATIC EQUATION (A * x * x + B * x + C = 0)
	//	// Equation giving distance between ray intersection with (*tree) crown and target cell center

	//	   // Compute a, b and c coefficients
	//	double coef_a = this->cosHeightAngle * this->cosHeightAngle * this->cosAzimuth * this->cosAzimuth / (a * a) +
	//		this->cosHeightAngle * this->cosHeightAngle * this->sinAzimuth * this->sinAzimuth / (b * b) +
	//		this->sinHeightAngle * this->sinHeightAngle / (c * c);

	//	double coef_b = -2 * x * this->cosHeightAngle * this->cosAzimuth / (a * a) -
	//		2 * y * this->sinAzimuth * this->cosHeightAngle / (b * b) -
	//		2 * z * this->sinHeightAngle / (c * c);

	//	double coef_c = (x * x) / (a * a) + (y * y) / (b * b) + (z * z) / (c * c) - 1;


	//	// Find positive solutions of the above quadratic equations (distance to target cell)
	//	// if not 2 solutions, set to NaN
	//	// Because negative = no interception and null = 1 interception = do not consider tangent rays with crown
	//	double delta = coef_b * coef_b - 4 * coef_a * coef_c;

	//	double sol1 = std::numeric_limits<double>::quiet_NaN();
	//	double sol2 = std::numeric_limits<double>::quiet_NaN();

	//	if (delta > 0) {
	//		double delta_sqrt = sqrt(delta);
	//		sol1 = (-coef_b + delta_sqrt) / (2 * coef_a);
	//		sol2 = (-coef_b - delta_sqrt) / (2 * coef_a);
	//	}



	//	// COMPUTE COORDINATES OF INTERCEPTION POINTS WITH CROWN LIMITS AND TARGET POINT

	//		// interception coords of the two roots
	//	double v1_x = sol1 * this->cosHeightAngle * this->cosAzimuth;
	//	double v1_y = sol1 * this->cosHeightAngle * this->sinAzimuth;
	//	double v1_z = sol1 * this->sinHeightAngle;

	//	double v2_x = sol2 * this->cosHeightAngle * this->cosAzimuth;
	//	double v2_y = sol2 * this->cosHeightAngle * this->sinAzimuth;
	//	double v2_z = sol2 * this->sinHeightAngle;

	//	// FIND SOLUTIONS

	//		// Get the limits of the box between shifted paraboloid center and crown limit
	//	double x_bbox_max = x + a;
	//	double y_bbox_max = y + b;
	//	double z_bbox_max = z + c;

	//	double x_bbox_min = x - a;
	//	double y_bbox_min = y - b;
	//	double z_bbox_min = z - c;

	//	// Find interception points within the crown (edge of ellipsoid or z axis if half ellipsoid = crown base plane)

	//	bool is_v1 =
	//		// interception 2 with crown is above the ground AND
	//		v1_z >= 0 &&
	//		v1_x >= x_bbox_min && v1_x <= x_bbox_max &&
	//		v1_y >= y_bbox_min && v1_y <= y_bbox_max &&
	//		// interception 2 with crown is in crown bbox AND
	//		v1_z >= z_bbox_min && v1_z <= z_bbox_max;


	//	bool is_v2 =
	//		// interception 2 with crown is above the ground AND
	//		v2_z >= 0 &&
	//		// interception 2 with crown is in crown bbox AND
	//		v2_x >= x_bbox_min && v2_x <= x_bbox_max &&
	//		v2_y >= y_bbox_min && v2_y <= y_bbox_max &&
	//		v2_z >= z_bbox_min && v2_z <= z_bbox_max;


	//	// FIND PATH LENGTH AND DISTANCE TO CELL CENTER
	//	// Distance between target cell and middle point between two interceptions
	//	// Length is the length of the ray path across the crown
	//	double distance = 0.0;
	//	double length = 0.0;

	//	// If 0 solution, no interception
	//	if (!(!std::isnan(sol1) && is_v1) && !(!std::isnan(sol2) && is_v2)) {
	//		return;
	//	}

	//	std::stack<double> distance_tmp;
	//	std::stack<double> length_tmp;

	//	int nsols = 0;
	//	if (!std::isnan(sol1) && is_v1) {
	//		nsols++;
	//		distance_tmp.push(sol1);
	//		length_tmp.push(sol1);
	//	}

	//	if (!std::isnan(sol2) && is_v2) {
	//		nsols++;
	//		distance_tmp.push(sol2);
	//		length_tmp.push(sol2);
	//	}

	//	// Compute length and distance if 2 solutions
	//	if (nsols == 2) {

	//		// Length
	//		double l1 = length_tmp.top();
	//		length_tmp.pop();
	//		double l2 = length_tmp.top();

	//		length = abs(l1 - l2);

	//		// Distance
	//		double d1 = distance_tmp.top();
	//		distance_tmp.pop();
	//		double d2 = distance_tmp.top();

	//		distance = (d1 + d2) / 2.0;
	//	}
	//	// Otherwise, there is a problem
	//	else {
	//		std::cout << "Not 0 or 2 solutions" << std::endl;
	//	}

	//	// Create and push interception to the vector of ray interceptions
	//	#ifdef _OPENMP
	//		#pragma omp critical
	//	#endif
	//	{
	//		this->interceptions.emplace_back(interception(target_cell, tree, length, distance, shifted_cell.outside));
	//	}
	//}

};

class RayManager {

private:
	// Vector of rays
	std::vector<Ray*> rays;

	// Energy coming from diffuse and direct rays above canopy (in MJ/m2)
	double energyAboveM2; 

public:
	RayManager(DataFrame rays, double e_above_m2) {

		// Energy per m2 above canopy
		this->energyAboveM2 = e_above_m2;

		// Number of rays
		int n_rays = rays.nrows();

		// Get vectors from the R dataframe
		IntegerVector id = rays["id_ray"];
		NumericVector energy = rays["e_incident"];
		NumericVector height_angle = rays["height_angle"];
		NumericVector azimuth = rays["azimut"];

		// Create a Tree object and store it in a vector of trees
		for (int i = 0; i < n_rays; i++) {
			Ray* ray = new Ray(id[i], energy[i], height_angle[i], azimuth[i]);
			this->rays.push_back(ray);
		}
	}

	~RayManager() {
		// Delete rays pointers
		int n_rays = this->rays.size();
		for (int i = 0; i < n_rays; i++) {
			delete this->rays[i];
		}
	}

	// Getters
	Ray* getRay(int id) { return(this->rays[id]); }
	int getNRays() { return(this->rays.size()); }
	double getEnergyAboveM2() { return(this->energyAboveM2); }
};

class CrownPart {

protected:
	// Id of the crown part and of the tree to which the crown part belong
	int idTree;
	int id;

	// Position of the center of the crown (in meters)
	double x;
	double y;
	double z;

	// Leaves charcteristics
	double crownOpeness; // between 0 and 1, percentage of ligth energy of a ray that is not absorbed by the crown (for porous envelop method)
	double crownLAD; // Leaf area density in m2/m3 (for turbid medium method)

	// Output energy (in MJ)
	double energy; // Intercepted energy
	double energyPotential; // Intercepted energy without considering neighbours


protected:
	// Get point of interception
	vertex3D getInterceptionPoint(double solution, Ray* ray) {

		vertex3D point = {
			solution * ray->getCosHeightAngle() * ray->getCosAzimuth(),
			solution * ray->getCosHeightAngle() * ray->getSinAzimuth(),
			solution * ray->getSinHeightAngle()
		};

		return(point);
	}

	// Two points are equal
	bool areEquals(vertex3D p1, vertex3D p2, double EPSILON) {
		return(abs(p1.x - p2.x) < EPSILON &&
			abs(p1.y - p2.y) < EPSILON &&
			abs(p1.z - p2.z) < EPSILON);
	}

	// Functions for finding if a point is in a volume
	bool isInBBox(vertex3D pmin, vertex3D pmax, vertex3D p, double EPSILON) {

		return (pmin.x - p.x <= EPSILON && p.x - pmax.x <= EPSILON &&
			pmin.y - p.y <= EPSILON && p.y - pmax.y <= EPSILON &&
			pmin.z - p.z <= EPSILON && p.z - pmax.z <= EPSILON);
	}

	bool isInEllipsoid(vertex3D p0, vertex3D p, double a, double b, double c) {

		double value = (p.x - p0.x) * (p.x - p0.x) / (a * a) +
			(p.y - p0.y) * (p.y - p0.y) / (b * b) +
			(p.z - p0.z) * (p.z - p0.z) / (c * c);

		return (value < 1.0);
	}


public:

	CrownPart(int id_tree, int id, 
		double x, double y, double z,
		double crown_openess, double crown_lad)
	{
		this->idTree = id_tree;
		this->id = id;
		this->x = x;
		this->y = y;
		this->z = z;
		this->crownOpeness = crown_openess;
		this->crownLAD = crown_lad;

		this->energy = 0.0;
		this->energyPotential = 0.0;
	}

	virtual ~CrownPart() {}

	// Getters
	int getId() { return(this->id); }
	int getIdTree() { return(this->idTree); }
	double getX() { return(this->x); }
	double getY() { return(this->y); }
	double getZ() { return(this->z); }
	double getCrownOpeness() { return(this->crownOpeness); }
	double getCrownLAD() { return(this->crownLAD); }
	double getEnergy() { return(this->energy); }
	double getEnergyPotential() { return(this->energyPotential); }

	// Setter
	void addEnergy(double e) { this->energy += e; }
	void addEnergyPotential(double epot) { this->energyPotential += epot; }

	// interception generic function
	virtual interception* computeInterception(int id_target_cell, Ray* ray, vertex3D& shift) { return(nullptr); }

};

class CrownPartEllipsoid8th : public CrownPart {

private:
	// Ellipsoid parameters (semi-principal axes)
	// The signs of a, b and c determine which 8th is considered.
	double a;
	double b;
	double c;

public:
	CrownPartEllipsoid8th(int id_tree, int id,
		double x, double y, double z,
		double a, double b, double c,
		double crown_openess, double crown_lad) :
		CrownPart(id_tree, id, x, y, z, crown_openess, crown_lad) {

		this->a = a;
		this->b = b;
		this->c = c;
	}

	// Getters
	double getA() { return(this->a); }
	double getB() { return(this->b); }
	double getC() { return(this->c); }


	// interception of a 8th of an ellispoid between a shifted tree and a ray coming toward a target cell
	interception* computeInterception(int id_target_cell, Ray* ray, vertex3D& shift) {

		// SHIFT POSITION OF THE ELLIPSOID
		// 2 shifts: the first one to set the target cell as origine, and the second to acoount for outside (*tree) when torus system
		vertex3D p0_shift = {
			this->x + shift.x,
			this->y + shift.y,
			this->z + shift.z
		};

		// FIND SOLUTION OF THE QUADRATIC EQUATION (A * x * x + B * x + C = 0)
		// Equation giving distance between ray intersection with (*tree) crown and target cell center
		// Compute a, b and c coefficients
		double coef_a = ray->getCosHeightAngle() * ray->getCosHeightAngle() * ray->getCosAzimuth() * ray->getCosAzimuth() / (this->a * this->a) +
			ray->getCosHeightAngle() * ray->getCosHeightAngle() * ray->getSinAzimuth() * ray->getSinAzimuth() / (this->b * this->b) +
			ray->getSinHeightAngle() * ray->getSinHeightAngle() / (this->c * this->c);

		double coef_b = -2.0 * p0_shift.x * ray->getCosHeightAngle() * ray->getCosAzimuth() / (this->a * this->a) -
			2.0 * p0_shift.y * ray->getSinAzimuth() * ray->getCosHeightAngle() / (this->b * this->b) -
			2.0 * p0_shift.z * ray->getSinHeightAngle() / (this->c * this->c);

		double coef_c = (p0_shift.x * p0_shift.x) / (this->a * this->a) +
			(p0_shift.y * p0_shift.y) / (this->b * this->b) +
			(p0_shift.z * p0_shift.z) / (this->c * this->c) - 1.0;


		// Find point of intersection between the ray and the whole ellipsoid
		// Find positive solutions of the above quadratic equations (distance to target cell)
		// if not 2 solutions, leave the function
		// Because negative = no interception and null = 1 interception = do not consider tangent rays with crown
		double delta = coef_b * coef_b - 4.0 * coef_a * coef_c;

		if (delta <= 0.0) { return nullptr; }

		double delta_sqrt = sqrt(delta);
		double solc1 = (-coef_b + delta_sqrt) / (2.0 * coef_a);
		double solc2 = (-coef_b - delta_sqrt) / (2.0 * coef_a);


		// COMPUTE COORDINATES OF INTERCEPTION POINTS WITH CROWN LIMITS AND TARGET POINT

			// Coordinates of the two interception points with the full ellipsoid
		vertex3D pc1 = getInterceptionPoint(solc1, ray);
		vertex3D pc2 = getInterceptionPoint(solc2, ray);

		// interception with respectively plane x = x0, y = y0, z = z0
		// And correct for rounding errors
		double solx0 = p0_shift.x / (ray->getCosHeightAngle() * ray->getCosAzimuth());
		double soly0 = p0_shift.y / (ray->getCosHeightAngle() * ray->getSinAzimuth());
		double solz0 = p0_shift.z / ray->getSinHeightAngle();

		vertex3D px0 = getInterceptionPoint(solx0, ray);
		vertex3D py0 = getInterceptionPoint(soly0, ray);
		vertex3D pz0 = getInterceptionPoint(solz0, ray);

		px0.x = p0_shift.x;
		py0.y = p0_shift.y;
		pz0.z = p0_shift.z;

		// interception with the target point
		vertex3D ptarget = getInterceptionPoint(0.0, ray);


		// FIND SOLUTIONS OF INTERSECTIO WITH THE EIGHTH ELLIPSOID

		// Get the limits of the bounding box of the 8th paraboloid
		double xa = p0_shift.x + a;
		double yb = p0_shift.y + b;
		double zc = p0_shift.z + c;

		vertex3D p_bbox_min = { std::min(p0_shift.x, xa), std::min(p0_shift.y, yb), std::min(p0_shift.z, zc) };
		vertex3D p_bbox_max = { std::max(p0_shift.x, xa), std::max(p0_shift.y, yb), std::max(p0_shift.z, zc) };

		double EPSILON = 1e-10; // accounts for rounding errors in boolean operations below
		std::vector<double> sols;

		// The first interception point with the full ellipsoid is in the eight ellipsoid
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc1, EPSILON) && pc1.z >= 0.0)
			sols.push_back(solc1);

		// The second interception point with the full ellipsoid is in the eight ellipsoid
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc2, EPSILON) && pc2.z >= 0.0)
			sols.push_back(solc2);

		// The ray intersects the X vertical plane of the eight ellipsoid
		if (this->isInBBox(p_bbox_min, p_bbox_max, px0, EPSILON) && this->isInEllipsoid(p0_shift, px0, this->a, this->b, this->c) && px0.z > 0.0)
			sols.push_back(solx0);

		// The ray intersects the Y vertical plane of the eight ellipsoid
		if (this->isInBBox(p_bbox_min, p_bbox_max, py0, EPSILON) && this->isInEllipsoid(p0_shift, py0, this->a, this->b, this->c) && py0.z > 0.0 && !this->areEquals(py0, px0, EPSILON))
			sols.push_back(soly0);

		// The ray intersects the Z horizontal plane of the eight ellipsoid
		if (this->isInBBox(p_bbox_min, p_bbox_max, pz0, EPSILON) && this->isInEllipsoid(p0_shift, pz0, this->a, this->b, this->c) && pz0.z > 0.0 && !this->areEquals(pz0, px0, EPSILON) && !this->areEquals(pz0, py0, EPSILON))
			sols.push_back(solz0);

		// The ray intersects the center of the target cell
		if (this->isInBBox(p_bbox_min, p_bbox_max, ptarget, EPSILON) && this->isInEllipsoid(p0_shift, ptarget, this->a, this->b, this->c))
			sols.push_back(0.0);


		// FIND PATH LENGTH AND DISTANCE TO CELL CENTER
		int nsols = sols.size();

		// intersected the ellipse outside the half quarter
		if (nsols == 0) { return(nullptr); }

		// If not 0 or 2 solutions ==> problem
		if (sols.size() != 2) {
			std::cout << "Not 0 or 2 solutions" << std::endl;
			return nullptr;
		}

		// Distance between target cell and middle point between two interceptions
		// Length is the length of the ray path across the crown
		double sol1 = sols.back();
		sols.pop_back();
		double sol2 = sols.back();

		double length = abs(sol1 - sol2);
		double distance = (sol1 + sol2) / 2.0;

		// Create interception 
		interception* interc = new interception { ray->getId(), id_target_cell, this->idTree, this->id, length, distance };

		return(interc);
	}

};

class Crown {

private:

	// Id of the tree to which the crown belong
	int idTree;

	// Radius of the crown
	double radiusMax;
	double radiusMean;

	// Crown parts in the crown
	std::map<int, CrownPart*> crownParts;


private:

	// --- FUNCTIONS TO INITIALIZE DIFFERENT CROWN FORMS ---

	// Irregular crown composed of 8 eighth of ellipsoids
	void initCrownEllipsoid8th(double x, double y, double z,
		double h, double hbase, double hmax,
		double cr_n, double cr_e, double cr_s, double cr_w,
		double crown_openess, double crown_lad) {

		// Get center of the crown 
		// x, y, z are tree coordinates
		double x0 = x;
		double y0 = y;
		double z0 = z + hmax;

		// Get depth of respectively top and bottom parts of the crown
		double depth_up = h - hmax;
		double depth_down = hmax - hbase;

		// Init top eighth parts of the crown
		if (depth_up > 0) {
			// Create pointors
			this->crownParts[1] = new CrownPartEllipsoid8th(this->idTree, 1,
				x0, y0, z0, cr_e, cr_n, depth_up, crown_openess, crown_lad);
			this->crownParts[2] = new CrownPartEllipsoid8th(this->idTree, 2,
				x0, y0, z0, cr_e, -cr_s, depth_up, crown_openess, crown_lad);
			this->crownParts[3] = new CrownPartEllipsoid8th(this->idTree, 3,
				x0, y0, z0, -cr_w, -cr_s, depth_up, crown_openess, crown_lad);
			this->crownParts[4] = new CrownPartEllipsoid8th(this->idTree, 4,
				x0, y0, z0, -cr_w, cr_n, depth_up, crown_openess, crown_lad);
		}

		// Init bottom eighth parts of the crown
		if (depth_down > 0) {
			// Create pointors
			this->crownParts[5] = new CrownPartEllipsoid8th(this->idTree, 5,
				x0, y0, z0, cr_e, cr_n, -depth_down, crown_openess, crown_lad);
			this->crownParts[6] = new CrownPartEllipsoid8th(this->idTree, 6,
				x0, y0, z0, cr_e, -cr_s, -depth_down, crown_openess, crown_lad);
			this->crownParts[7] = new CrownPartEllipsoid8th(this->idTree, 7,
				x0, y0, z0, -cr_w, -cr_s, -depth_down, crown_openess, crown_lad);
			this->crownParts[8] = new CrownPartEllipsoid8th(this->idTree, 8,
				x0, y0, z0, -cr_w, cr_n, -depth_down, crown_openess, crown_lad);
		}
	}


public:
	Crown(int id_tree, 
		std::string crown_type,
		double x, double y, double z,
		double h, double hbase, double hmax,
		double cr_n, double cr_e, double cr_s, double cr_w,
		double crown_openess, double crown_lad)
	{
		this->idTree = id_tree;

		// Maximum and average radius 
		this-> radiusMax = std::max({ cr_n, cr_e, cr_s, cr_w });
		//this->radiusMean = std::ave({ cr_n, cr_e, cr_s, cr_w });

		// Init the crown parts depending on the type
		if (crown_type == "8E") {
			this->initCrownEllipsoid8th(x, y, z,
				h, hbase, hmax, cr_n, cr_e, cr_s, cr_w,
				crown_openess, crown_lad);
		}
		else {
			std::cout << "Unrecognized crown type" << std::endl;
		}

	}

	~Crown() {
		// Delete crown parts pointers
		for (auto const& x: this->crownParts) {
			delete x.second;
		}
	}

	// Getters
	CrownPart* getCrownPart(int id) { return(this->crownParts[id]); }
	double getRadiusMax() { return(this->radiusMax); }

	double getEnergy() { 
		double energy = 0.0;
		for (auto const& x : this->crownParts) {
			energy += x.second->getEnergy();
		}
		return(energy);
	}

	double getEnergyPotential() {
		double energy_pot = 0.0;
		for (auto const& x : this->crownParts) {
			energy_pot += x.second->getEnergyPotential();
		}
		return(energy_pot);
	}


	// Find interceptions of each parts of the crown
	std::vector<interception*> computeInterception(int id_target_cell, Ray* ray, vertex3D& shift) {
		
		std::vector<interception*> interceptions_crown;

		//#ifdef _OPENMP
		//	#pragma omp parallel for
		//#endif
		for (auto const& x : this->crownParts) {
			interception* interc = x.second->computeInterception(id_target_cell, ray, shift);
			if (interc != nullptr) {
				//#ifdef _OPENMP
				//	#pragma omp critical
				//#endif
				interceptions_crown.push_back(interc);
			}
		}
		return(interceptions_crown);
	}
};

class Tree {

private:
	// Id of the tree
	int id;

	// Tree position
	double x;
	double y;
	double z;

	// Tree dimensions
	double dbh;
	double height;

	// Crown
	Crown crown;


public:
	Tree(int id, double x, double y, double z,
		double dbh, double h,
		std::string crown_type,
		double hbase, double hmax,
		double cr_n, double cr_e, double cr_s, double cr_w,
		double crown_openess, double crown_lad):

		crown(id, crown_type, x, y, z, h, hbase, hmax, cr_n, cr_e, cr_s, cr_w,
			crown_openess, crown_lad)
	{
		// Init tree variables
		this->id = id;
		this->x = x;
		this->y = y;
		this->z = z;

		// Init tree dimensions
		this->dbh = dbh;
		this->height = h;

	}


	// Getters
	int getId() { return(this->id); }
	double getX() { return(this->x); }
	double getY() { return(this->y); }
	double getZ() { return(this->z); }
	double getHeight() { return(this->height); }
	double getDbh() { return(this->dbh); }
	double getRadiusMax() { return(this->crown.getRadiusMax()); }

	CrownPart* getCrownPart(int id) { return(this->crown.getCrownPart(id)); }

	double getEnergy() { return(this->crown.getEnergy()); }
	double getEnergyPotential() { return(this->crown.getEnergyPotential()); }


	// Methods
	std::vector<interception*> computeInterception(int id_target_cell, Ray* ray, vertex3D& shift) {
		return(this->crown.computeInterception(id_target_cell, ray, shift));
	}
};

class Cell {

private:
	// Id of the cell
	int id;

	// Column and row in a grid system (i.e. Stand object)
	int row;
	int col;

	// Position of the center of the cell
	double x;
	double y;
	double z;

	// Trees in the cell
	std::vector<Tree*> trees;

	// Energy of the cell (in MJ)
	double energyAboveCanopy; // Energy coming to the cell above the canopy (energy of the cell without trees above)
	double energy; // Current energy above the cell (will be decreased by successive interception by the abov trees


public:
	Cell(int id, int row, int col, double x, double y, double z, double e_above) {
		this->id = id;
		this->row = row;
		this->col = col;
		this->x = x;
		this->y = y;
		this->z = z;

		// Energy of the cell (in MJ)
		this->energyAboveCanopy = e_above; // Energy coming to the cell above the canopy (energy of the cell without trees above)
		this->energy = e_above; // Current energy above the cell (will be decreased by successive interception by the abov trees
	}

	~Cell() {
		// Delete trees pointers
		int n_trees = this->trees.size();
		for (int i = 0; i < n_trees; i++) {
			delete this->trees[i];
		}
	}

	// Getters
	int getId() { return(this->id); }
	int getRow() { return(this->row); }
	int getCol() { return(this->col); }
	double getX() { return(this->x); }
	double getY() { return(this->y); }
	double getZ() { return(this->z); }

	bool isEmpty() { return(this->trees.empty()); }
	int getNTrees() { return(this->trees.size()); }
	Tree* getTree(int i) { return(this->trees[i]); }
	void printTreesId() {
		int n_trees = this->trees.size();
		for (int i = 0; i < n_trees; i++) {
			std::cout << (*this->trees[i]).getId() << "\t";
		}
	}

	double getEnergy() { return(this->energy); }
	double getEnergyM2(double cell_area) { return(this->energy / cell_area); }
	double getEnergyRelative() { return(this->energy / this->energyAboveCanopy); }


	// Setters
	void addTree(Tree* tree) { this->trees.push_back(tree); }
	void interceptEnergy(double e) { this->energy -= e; }

};

class Stand {

private:
	// 2D map of cell pointers (grid stand) stored by <row, column>
	std::map<std::pair<int, int>, Cell*> grid;

	// Vector of cells and trees stored by unique id
	std::map<int, Cell*> cells;
	std::map<int, Tree*> trees;

	// Number of trees
	int nTrees;

	// Tree maximum sizes
	double maxTreeHeight;
	double maxTreeCrownRadius;

	// Stand orientation (all in radians)
	double slope;
	double northToXcw;
	double aspect;
	double bottomAzimuth; // Azimuth of the vector orthogonal to the ground in the (x, y) system

	// Stand size
	double cellSize; // Length of the size of a single cell (in m)
	int nCells; // Number of cells in a row (in m)
	double cellAreaPlane; // Area	of the cell at horizontal (in m2)
	double cellArea; // Area of the cell considering slope (in m2)
	bool useTorus; // Whether we want to use a torus system for light ray tracing


public:
	Stand(DataFrame cells, DataFrame trees,
		double slope, double north_to_x_cw, double aspect,
		double cell_size, double n_cells, bool use_torus,
		double e_above_m2)
	{
		// Stand orientation (all in radians)
		this->slope = slope * M_PI / 180.0;
		this->northToXcw = north_to_x_cw * M_PI / 180.0;
		this->aspect = aspect * M_PI / 180.0;
		this->bottomAzimuth = (-aspect + north_to_x_cw) * M_PI / 180.0;

		// Stand size
		this->cellSize = cell_size;
		this->nCells = n_cells;
		this->cellAreaPlane = cell_size * cell_size;
		this->cellArea = this->cellAreaPlane / cos(this->slope);
		this->useTorus = use_torus;

		// Number of cells and trees
		int n_cells_grid = cells.nrows();
		this->nTrees = trees.nrows();

		// Get vectors from the cells R dataframe
		IntegerVector cells_id = cells["id_cell"];
		IntegerVector cells_col = cells["x_id"];
		IntegerVector cells_row = cells["y_id"];
		NumericVector cells_x = cells["x_center"];
		NumericVector cells_y = cells["y_center"];
		NumericVector cells_z = cells["z_center"];

		// Get vectors from the trees R dataframe
		IntegerVector trees_id = trees["id_tree"];
		IntegerVector trees_id_cell = trees["id_cell"];
		NumericVector trees_x = trees["x"];
		NumericVector trees_y = trees["y"];
		NumericVector trees_z = trees["z"];
		NumericVector dbh = trees["dbh_cm"];
		NumericVector height = trees["h_m"];
		NumericVector hbase = trees["hbase_m"];
		NumericVector hmax = trees["hmax_m"];
		NumericVector cr_n = trees["rn_m"];
		NumericVector cr_e = trees["re_m"];
		NumericVector cr_s = trees["rs_m"];
		NumericVector cr_w = trees["rw_m"];
		StringVector ctype = trees["crown_type"];
		NumericVector cp = trees["crown_openess"];
		NumericVector clad = trees["crown_lad"];

		// Get energy above canopy that is coming toward a cell
		double e_above_cell = e_above_m2 * this->cellSize;

		// Create a Cell pointer and store it in a vector and grid of target cells
		for (int i = 0; i < n_cells_grid; i++) {
			Cell* cell = new Cell(cells_id[i] - 1, cells_row[i], cells_col[i], cells_x[i], cells_y[i], cells_z[i], e_above_cell);
			this->cells.emplace(cell->getId(), cell);
			this->grid.emplace(std::make_pair(cell->getRow(), cell->getCol()), cell);
		}

		// Add a Tree object to the associated Cell they belong to
		for (int i = 0; i < this->nTrees; i++) {
			Tree* tree = new Tree(
				trees_id[i], trees_x[i], trees_y[i], trees_z[i],
				dbh[i], height[i], Rcpp::as< std::string >(ctype[i]), hbase[i], hmax[i],
				cr_n[i], cr_e[i], cr_s[i], cr_w[i], cp[i], clad[i]);
			this->trees.emplace(tree->getId(), tree);
			this->cells.find(trees_id_cell[i] - 1)->second->addTree(tree);
		}

		// Get maximum height and crown radius of the trees
		this->maxTreeHeight = max(height);
		this->maxTreeCrownRadius = std::max({ max(cr_n), max(cr_e), max(cr_s), max(cr_w) });
	}

	~Stand() {
		// Delete cells pointers (also destroy tree pointors within th given cell)
		int n_cells = this->cells.size();
		for (int i = 0; i < n_cells; i++) {
			delete this->cells[i];
		}
	}

	// Getters
	Cell* getCell(int row, int col) { return(this->grid[std::make_pair(row, col)]); }
	Cell* getCell(int id) { return(this->cells[id]); }
	std::map<int, Cell*> getCells() { return(this->cells); }
	int getNCells() { return(this->nCells); }

	Tree* getTree(int id) { return(this->trees[id]); }
	int getNTrees() { return(this->nTrees); }

	double getSlope() { return(this->slope); }
	double getNorthToXcw() { return(this->northToXcw); }
	double getAspect() { return(this->aspect); }
	double getBottomAzimuth() { return(this->bottomAzimuth); }
	double getCellSize() { return(this->cellSize); }
	double getCellAreaPlane() { return(this->cellAreaPlane); }
	double getCellArea() { return(this->cellArea); }
	double isUseTorus() { return(this->useTorus); }
	double getMaximumTreeCrownRadius() { return(this->maxTreeCrownRadius); }
	double getMaximumTreeHeight() { return(this->maxTreeHeight); }

};

struct shiftedCell {

	// Pointor toward the original Cell
	Cell* cell;

	// Shift applied to the cell position (in a torus system)
	vertex3D shift;
};

class InterceptionManager {

private:
	// Map of all interceptions per ray per target cell
	// First int of the pair key is the target cell, the second one is the ray
	std::map<std::pair<int, int>, std::vector<interception*>> interceptions;

	// Light interception method by the crowns : either porous envelop(use "crownOpenness") or turbid medium(use "crownLAD")
	bool turbidMedium;

	// Global variable for turbid medium
	const double EXTINCTION_COEF; // Probability of a leaf to intercept the ray (linked to leaf orientation)
	const double CLUMPING_FACTOR; // Aggregation of leaves within the crown volume (1 is homogeneous)

private:

	double applyBeerLambert(double incident_energy, double extinction_coef, double clumping_factor, double leaf_area_density, double path_length) {

		return(incident_energy * (1 - exp(-extinction_coef * clumping_factor * leaf_area_density * path_length)));

	}


public:
	InterceptionManager(bool turbid_medium):
		EXTINCTION_COEF(0.5), CLUMPING_FACTOR(1.0)
	{
		this->turbidMedium = turbid_medium;
	}

	~InterceptionManager() {
		// Delete interception pointers
		for (auto const& x : this->interceptions) {
			int n_intercs = x.second.size();
			for (int i = 0; i < n_intercs; i++) {
				delete x.second[i];
			}
		}
	}

	// Setters
	void addInterception(int cell_id, int ray_id, interception* interc) {
		this->interceptions[std::make_pair(cell_id, ray_id)].push_back(interc);
	}

	// Methods
	void orderInterceptions() {
		// Sort interceptions by distance to the target cell for each ray X target cell
		// Furthest first
		for (auto & x : this->interceptions) {
			std::sort(x.second.begin(), x.second.end(), [](interception const* i1, interception const* i2) {
				return i1->distance > i2->distance;
				});
		}
	}

	void summarizeInterceptions(RayManager& rays, Stand& stand) {

		// iterate over all rays going towards each target cell
		for (auto const& x : this->interceptions) {

			// Get ray and target cell of this interceptions vector
			Cell* target_cell = stand.getCell(x.first.first);
			Ray* ray = rays.getRay(x.first.second);

			// Compute projection of energy on plane parallel to slope
			double scalar_slope = cos(stand.getSlope()) * ray->getSinHeightAngle() +
				sin(stand.getSlope()) * ray->getCosHeightAngle() * cos(ray->getAzimuth() - stand.getBottomAzimuth());

			//  in MJ / m2 and convert it into MJ per cell
			double e_incident_slope_m2 = scalar_slope * ray->getIncidentEnergy();
			double e_incident_slope_cell = e_incident_slope_m2 * stand.getCellArea();

			// Initialize the energy of the ray coming toward the cell above the canopy
			double current_energy = e_incident_slope_cell;

			// Compute attenuation of energy across successives crown interceptions for each ray X target cell
			int n_interceptions = x.second.size();
			for (int i = 0; i < n_interceptions; i++) {

				// Get the interception
				CrownPart* crown_part = stand.getTree(x.second[i]->idTree)->getCrownPart(x.second[i]->idCrownPart);

				// Compute the potential energy to the tree (energy without attenuation by neighbours)
				// and compute energy with attenuation
				// Different considering crown as turbid medium or porous envelop
				double potential_energy = 0.0;
				double intercepted_energy = 0.0;


				// Turbid medium ==> apply beer lambert law
				if (this->turbidMedium) {

					potential_energy = this->applyBeerLambert(
						e_incident_slope_cell,
						this->EXTINCTION_COEF, this->CLUMPING_FACTOR,
						crown_part->getCrownLAD(),
						x.second[i]->length);

					intercepted_energy = this->applyBeerLambert(
						current_energy,
						this->EXTINCTION_COEF, this->CLUMPING_FACTOR,
						crown_part->getCrownLAD(),
						x.second[i]->length);
				}
				// Porous envelop ==> reduce the energy by a fixed amount
				else {
					potential_energy = e_incident_slope_cell * (1 - crown_part->getCrownOpeness());
					intercepted_energy = current_energy * (1 - crown_part->getCrownOpeness());
				}

				// Add to the potential and intercepted energy by the tree
				crown_part->addEnergyPotential(potential_energy);
				crown_part->addEnergy(intercepted_energy);

				// And remove the intercepted energy from the energy that left
				current_energy -= intercepted_energy;

				// Remove the intercepted energy from the total energy above the cell
				target_cell->interceptEnergy(intercepted_energy);
			}
		}

	}

};

class Model {

private:
	// Stand object with cells, trees and geometric variables
	Stand stand;

	// Rays that are coming toward a cell
	RayManager rays;

	// Interceptions between a crown part and a ray toward a target cell center
	InterceptionManager interceptions;


public:

	Model(DataFrame trees, DataFrame cells, DataFrame rays,
		double e_above_m2,
		double slope, double north_to_x_cw, double aspect,
		double cell_size, double n_cells,
		bool use_torus, bool turbid_medium) :

		stand(cells, trees, slope, north_to_x_cw, aspect, cell_size, n_cells, use_torus, e_above_m2),
		rays(rays, e_above_m2),
		interceptions(turbid_medium)
	{}

	// Getters
	double getMaximumTreeHeight() { return(this->stand.getMaximumTreeHeight()); }
	double getMaximumTreeCrownRadius() { return(this->stand.getMaximumTreeCrownRadius()); }


	// Methods
	void findPotentialRelCells() {
		
		// Potential relative cells for each ray coming to a cell
		int n_rays = this->rays.getNRays();
		for (int i = 0; i < n_rays; i++) {

			Ray* ray = this->rays.getRay(i);

			// Computes lateral = the boundary to add to the competition rectangle to take into account cells center instead of trees position.
			// The boundary depends on beam azimut.
			double azt = 0;
			if (ray->getAzimuth() < M_PI / 4.0) {
				azt = ray->getAzimuth();
			}
			else if (ray->getAzimuth() >= M_PI / 4.0 && ray->getAzimuth() < M_PI / 2.0) {
				azt = M_PI / 2.0 - ray->getAzimuth();
			}
			else if (ray->getAzimuth() >= M_PI / 2.0 && ray->getAzimuth() < 3.0 * M_PI / 4.0) {
				azt = ray->getAzimuth() - M_PI / 2.0;
			}
			else if (ray->getAzimuth() >= 3.0 * M_PI / 4.0 && ray->getAzimuth() < M_PI) {
				azt = M_PI - ray->getAzimuth();
			}
			else if (ray->getAzimuth() >= M_PI && ray->getAzimuth() < 5.0 * M_PI / 4.0) {
				azt = ray->getAzimuth() - M_PI;
			}
			else if (ray->getAzimuth() >= 5.0 * M_PI / 4.0 && ray->getAzimuth() < 3.0 * M_PI / 2.0) {
				azt = 3.0 * M_PI / 2.0 - ray->getAzimuth();
			}
			else if (ray->getAzimuth() >= 3.0 * M_PI / 2.0 && ray->getAzimuth() < 7.0 * M_PI / 4.0) {
				azt = ray->getAzimuth() - 3.0 * M_PI / 2.0;
			}
			else if (ray->getAzimuth() >= 7.0 * M_PI / 4.0) {
				azt = 2.0 * M_PI - ray->getAzimuth();
			}

			double lateral = this->stand.getCellSize() / sqrt(2.0) * sin(azt + M_PI / 4.0);

			// Beam width = max lateral distance from the beam to a cell center able to own a tree which can intercept the beam.
			double R = this->stand.getMaximumTreeCrownRadius() + lateral;

			// Beam reach maximum distance along the beam beyond which the cells cannot own trees which can intercept the beam (too high).
			double L = this->stand.getMaximumTreeHeight() / (tan(ray->getHeightAngle()) + cos(ray->getAzimuth() - this->stand.getBottomAzimuth()) * tan(this->stand.getSlope())) + lateral;

			// Coordinates of the four corners of the competition rectangle.
			double sinA = sin(ray->getAzimuth());
			double cosA = cos(ray->getAzimuth());

			double x1 = R * sinA + L * cosA;
			double y1 = L * sinA - R * cosA;
			double x2 = L * cosA - R * sinA;
			double y2 = L * sinA + R * cosA;
			double x3 = R * (sinA - cosA);
			double y3 = -R * (sinA + cosA);
			double x4 = -R * (sinA + cosA);
			double y4 = R * (cosA - sinA);

			double x_min = std::min({ x1, x2, x3, x4 });
			double x_max = std::max({ x1, x2, x3, x4 });
			double y_min = std::min({ y1, y2, y3, y4 });
			double y_max = std::max({ y1, y2, y3, y4 });

			// Round into relative-to-target coordinates of the cell center in which x/y_min/max are
			x_min = std::ceil(x_min / this->stand.getCellSize()) * this->stand.getCellSize();
			x_max = std::floor(x_max / this->stand.getCellSize()) * this->stand.getCellSize();
			y_min = std::ceil(y_min / this->stand.getCellSize()) * this->stand.getCellSize();
			y_max = std::floor(y_max / this->stand.getCellSize()) * this->stand.getCellSize();

			// Number of cells between min and max for both axis x and y
			int nx = std::floor((x_max - x_min) / this->stand.getCellSize() + 0.5) + 1;
			int ny = std::floor((y_max - y_min) / this->stand.getCellSize() + 0.5) + 1;

			// For each candidate relative cell of the given ray
			for (int xid = 0; xid < nx; xid++) {
				for (int yid = 0; yid < ny; yid++) {

					// Coordinates of the candidate cells (be careful : y coordinates is inverse direction of y grid)
					// id = 0 is the greater y coordinates
					double x = x_min + xid * this->stand.getCellSize();
					double y = y_max - yid * this->stand.getCellSize();

					// Add as candidate cell if it is located inside the competition rectangle
					if ((x * sinA - y * cosA < R) &&
						(x * cosA + y * sinA < L) &&
						(-x * sinA + y * cosA < R) &&
						(x * cosA + y * sinA > -R)) {

						// Relative columns and rows of the cell compared from the origin
						int col = std::floor(x / this->stand.getCellSize());
						int row = -std::floor(y / this->stand.getCellSize());

						// Add relative cell as a candidate
						ray->addPotentialRelCell(row, col);
					}
				}
			}
		}
	}

	shiftedCell getPotentialCell(Ray* ray, int rel_cell_id, Cell* target_cell) {

		relativeCoords rel_cell = ray->getPotentialRelCell(rel_cell_id);

		// Get row and column of potential cell
		int row_pot = target_cell->getRow() + rel_cell.row;
		int col_pot = target_cell->getCol() + rel_cell.col;

		// Search if cells is outside the this->stand
		bool is_outside = row_pot < 0 || row_pot >= this->stand.getNCells() || col_pot < 0 || col_pot >= this->stand.getNCells();

		// If cel is outside the main this->stand, do not keep as potential if we are not in a torus system
		if (is_outside && !this->stand.isUseTorus()) { return(shiftedCell{ nullptr, 0.0, 0.0, 0.0 }); }

		// Get row and column of the original cell in the grid cell (diffrent only within a torus system, if cell is outside)
		int row_pot_original = (this->stand.getNCells() + (row_pot % this->stand.getNCells())) % this->stand.getNCells();
		int col_pot_original = (this->stand.getNCells() + (col_pot % this->stand.getNCells())) % this->stand.getNCells();

		// Get corresponding original cell
		Cell* original_cell = this->stand.getCell(row_pot_original, col_pot_original);

		// Check if it contains trees, otherwise, do not add as potential cell
		if ((*original_cell).isEmpty()) { return(shiftedCell{ nullptr, 0.0, 0.0, 0.0 }); }

		// Compute shift in coordinates to apply on trees within the original cell (torus system)
		// Shift is a multiple of stand_size : how many this->stand size we have to shift tree coordinates times this->stand size, to apply torus system
		// If we have a potential cell with negative or gretaer than n_cells row or column
		// ==> Thus, we want to substract coordinates of original cells in order to have potential cell outside the main this->stand
		// Be careful, y - coordinates system is opposite direction of y - grid system
		double stand_size = this->stand.getCellSize() * this->stand.getNCells();

		double x_shift = (col_pot / this->stand.getNCells() - (col_pot % this->stand.getNCells() < 0 ? 1 : 0)) * stand_size; // Negative x_id ==> negative shift
		double y_shift = -(row_pot / this->stand.getNCells() - (row_pot % this->stand.getNCells() < 0 ? 1 : 0)) * stand_size; // Negative y_id ==> positive shift


		// Compute altitudinal shift (only if (x,y) shifted
		double z_shift = 0.0;
		if (x_shift != 0.0 || y_shift != 0.0) {
			double d = sqrt(x_shift * x_shift + y_shift * y_shift);
			double azimuth_xy = 0;
			if (d != 0) {
				if (y_shift >= 0) {
					azimuth_xy = acos(x_shift / d);
				}
				else {
					azimuth_xy = 2.0 * M_PI - acos(x_shift / d);
				}
			}
			z_shift = -d * cos(azimuth_xy - this->stand.getBottomAzimuth()) * tan(this->stand.getSlope());
		}

		// Return the shifted cell
		vertex3D shift = { x_shift, y_shift, z_shift };
		shiftedCell pot_cell = { original_cell, shift };

		return(pot_cell);
	}

	void computeInterceptions() {

		int n_rays = this->rays.getNRays();
		int n_cells = this->stand.getNCells();

		// For each target cell
		#ifdef _OPENMP
			#pragma omp parallel for
		#endif
		for (int c = 0; c < n_cells; c++) {

			Cell* target_cell = this->stand.getCell(c);

			// Get interceptions for all rays X all trees
			for (int r = 0; r < n_rays; r++) {

				Ray* ray = this->rays.getRay(r);

				// For each possible relative cell coordinates
				int n_relcells = ray->getNPotCells();
				for (int rc = 0; rc < n_relcells; rc++) {

					// Get a pointer to the original potential cell and its associated shift
					shiftedCell pot_cell = this->getPotentialCell(ray, rc, target_cell);

					// Cell is not potential if there is no trees or it is outside the main plot and torus system is disabled
					if (!pot_cell.cell) { continue; }

					// Compute the shift applied to the tree (2 shifts)
					// 1. Shift when torus system is enable and the tree is outside the main plot
					// 2. Set the center of the target cell as origin of tree and ray interception
					vertex3D shift_tree = {
						pot_cell.shift.x - target_cell->getX(),
						pot_cell.shift.y - target_cell->getY(),
						pot_cell.shift.z - target_cell->getZ()
					};

					// Compute interception for each tree within the potential cell (with a given shift if the potential cell is outside the plot)
					int n_pot_trees = pot_cell.cell->getNTrees();
					for (int t = 0; t < n_pot_trees; t++) {

						// Get tree pointer
						Tree* tree_p = pot_cell.cell->getTree(t);

						// Compute interception path length and distance form target cell if ray intercepted the crown parts of the tree
						std::vector<interception*> interceptions_tree = tree_p->computeInterception(target_cell->getId(), ray, shift_tree);

						// Store the interceptions
						int n_interceptions_tree = interceptions_tree.size();
						for (int i = 0; i < n_interceptions_tree; i++) {
							#ifdef _OPENMP
								#pragma omp critical
							#endif
							this->interceptions.addInterception(rc, r, interceptions_tree[i]);
						}
					}
				}
			}
		}
	}

	void orderInterceptions() {
		this->interceptions.orderInterceptions();
	}

	void summarizeInterceptions() {
		this->interceptions.summarizeInterceptions(this->rays, this->stand);
	}

	List exportResults() {

		// Get cells
		std::map<int, Cell*> cells = this->stand.getCells();

		// Get number of trees and cells
		int n_cells = cells.size();
		int n_trees = this->stand.getNTrees();

		// Init RCPP vectors for trees
		IntegerVector id_tree(n_trees);
		NumericVector e_tree(n_trees);
		NumericVector epot_tree(n_trees);
		NumericVector lci_tree(n_trees);

		// Init RCPP vectors for cells
		IntegerVector id_cell(n_cells);
		NumericVector e_cell(n_cells);
		NumericVector erel_cell(n_cells);

		// For each cell
		int icell = 0;
		int itree = 0;
		for (auto& cell : cells) {

			// Add cell to vectors
			id_cell[icell] = cell.second->getId() + 1;
			e_cell[icell] = cell.second->getEnergyM2(this->stand.getCellArea());
			erel_cell[icell] = cell.second->getEnergyRelative();

			// For each tree composing the cell
			int n_trees_cell = cell.second->getNTrees();
			for (int t = 0; t < n_trees_cell; t++) {

				Tree* tree = cell.second->getTree(t);

				id_tree[itree] = tree->getId();
				e_tree[itree] = tree->getEnergy();
				epot_tree[itree] = tree->getEnergyPotential();

				itree++;
			}

			icell++;
		}

		// Create trees and cells RCPP DataFrames
		DataFrame output_trees = DataFrame::create(
			Named("id_tree") = id_tree,
			Named("epot") = epot_tree,
			Named("e") = e_tree
		);

		DataFrame output_cells = DataFrame::create(
			Named("id_cell") = id_cell,
			Named("e") = e_cell,
			Named("erel") = erel_cell
		);

		// Return output as a List of two DataFrames
		return(List::create(
			Named("trees") = output_trees,
			Named("cells") = output_cells
		));
	}

};




// [[Rcpp::export]]
List sl_run_rcpp(
	DataFrame trees, DataFrame cells,
	DataFrame rays, double total_energy_m2,
	double slope, double north_to_x_cw, double aspect,
	double cell_size, double n_cells,
	bool use_torus, bool turbid_medium
)
{
	// TESTS:
	// - If plots is big enough for tree coordinates

	// TODO:
	// - Trunk interception
	// - Other forms E, P, C and 8P, 8C
	// - Not squared plot (rectangle or shape weird)


	// Initialize the model
	Model sl_model = Model(trees, cells,
		rays, total_energy_m2,
		slope, north_to_x_cw, aspect,
		cell_size, n_cells, use_torus, turbid_medium);

	// [OPTIMIZATION]: find the extend of possible interception of each ray
	// i.e. for each ray coming to an unknown target cell, find relative cells around the target containing tree that could possibly intercept the ray
	// We can whether or not consider a torus system
	sl_model.findPotentialRelCells();

	// For each target cell and each ray (target cells can be parallelized)
	sl_model.computeInterceptions();

	// Order the interceptions by distance to the target cell
	sl_model.orderInterceptions();

	// Get energy for each tree and each cell
	sl_model.summarizeInterceptions();

	// Convert teh output into R format
	List output = sl_model.exportResults();

	return(output);
}





