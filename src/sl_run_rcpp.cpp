#define _USE_MATH_DEFINES

constexpr double EPSILON = 10e-10; // For rounding errors

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

	// Id of the intercepted tree in the trees vector in the Stand object
	int vectIdTree;

	// Length of the path across the crown
	double length;

	// Distance between interception point (middle of full path) and target cell
	double distance;

	// Interception by crown or trunk ?
	bool withTrunk;
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
			Ray* ray = new Ray(id[i]-1, energy[i], height_angle[i], azimuth[i]);
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

class TreeVolume {

protected:
	// Id the tree to which the crown part belong
	int vectIdTree;

	// Position of the center of the crown (in meters)
	double x;
	double y;
	double z;

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
	bool areEquals(vertex3D p1, vertex3D p2) {
		return(abs(p1.x - p2.x) < EPSILON &&
			abs(p1.y - p2.y) < EPSILON &&
			abs(p1.z - p2.z) < EPSILON);
	}

	// Functions for finding if a point is in a volume
	bool isInBBox(vertex3D pmin, vertex3D pmax, vertex3D p) {

		return (pmin.x - p.x <= EPSILON && p.x - pmax.x <= EPSILON &&
			pmin.y - p.y <= EPSILON && p.y - pmax.y <= EPSILON &&
			pmin.z - p.z <= EPSILON && p.z - pmax.z <= EPSILON);
	}

public:
	TreeVolume(int vectid_tree, double x, double y, double z) {
		this->vectIdTree = vectid_tree;

		this->x = x;
		this->y = y;
		this->z = z;
	}

	int getVectIdTree() { return(this->vectIdTree); }
	double getX() { return(this->x); }
	double getY() { return(this->y); }
	double getZ() { return(this->z); }

	// interception generic function
	virtual interception* computeInterception(int id_target_cell, Ray* ray, vertex3D& shift) { return(nullptr); }

};

class CrownPart: public TreeVolume {

public:
	CrownPart(int vectid_tree,
		double x, double y, double z):
		TreeVolume(vectid_tree, x, y, z)
	{}

	virtual ~CrownPart() {}

	// interception generic function
	virtual interception* computeInterception(int id_target_cell, Ray* ray, vertex3D& shift) { return(nullptr); }
};

class CrownPartEllipsoid : public CrownPart {

private:

	// Do we consider 8th of volume, a semi-volume split on a horizontal plane or the full volume
	bool isSemi;
	bool is8th;

	// Ellipsoid parameters (semi-principal axes)
	// The signs of a, b and c determine which 8th is considered.
	double a;
	double b;
	double c;


private:
	bool isInEllipsoid(vertex3D p0, vertex3D p) {

		double value = (p.x - p0.x) * (p.x - p0.x) / (this->a * this->a) +
			(p.y - p0.y) * (p.y - p0.y) / (this->b * this->b) +
			(p.z - p0.z) * (p.z - p0.z) / (this->c * this->c);

		return (value < 1.0);
	}

public:
	CrownPartEllipsoid(int vectid_tree,
		bool is_semi, bool is_8th,
		double x, double y, double z,
		double a, double b, double c) :
		CrownPart(vectid_tree, x, y, z) {

		this->isSemi = is_8th ? true : is_semi; // If 8th thus semi ellipsoid is mandatory
		this->is8th = is_8th;

		this->a = a;
		this->b = b;
		this->c = c;
	}

	// interception of a 8th, semi or full ellispoid between a shifted tree and a ray coming toward a target cell
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

			// Interception with the horizontal base of the crown 
			// Checked only if semi or 8th ellispoid
		double solz0 = p0_shift.z / ray->getSinHeightAngle();
		vertex3D pz0 = getInterceptionPoint(solz0, ray);
		pz0.z = p0_shift.z;


			// If interception with a eighth of ellipsoid, consider interception with the X and Y vertical planes
			// Checked only in the case of a eihth of ellipsoid
		double solx0 = p0_shift.x / (ray->getCosHeightAngle() * ray->getCosAzimuth());
		double soly0 = p0_shift.y / (ray->getCosHeightAngle() * ray->getSinAzimuth());

		vertex3D px0 = getInterceptionPoint(solx0, ray);
		vertex3D py0 = getInterceptionPoint(soly0, ray);

		px0.x = p0_shift.x;
		py0.y = p0_shift.y;


			// interception with the target point
		vertex3D ptarget = getInterceptionPoint(0.0, ray);


		// FIND SOLUTIONS OF INTERSECTION WITH THE EIGHTH ELLIPSOID

		double xa_add = p0_shift.x + this->a;
		double yb_add = p0_shift.y + this->b;
		double zc_add = p0_shift.z + this->c;

		double xa_sub = p0_shift.x - this->a;
		double yb_sub = p0_shift.y - this->b;
		double zc_sub = p0_shift.z - this->c;

		// Get the limits of the bounding box of the 8th paraboloid
		vertex3D p_bbox_min = { 0.0, 0.0, 0.0 };
		vertex3D p_bbox_max = { 0.0, 0.0, 0.0 };
		if (this->is8th) {
			p_bbox_min = { std::min(p0_shift.x, xa_add), std::min(p0_shift.y, yb_add), std::min(p0_shift.z, zc_add) };
			p_bbox_max = { std::max(p0_shift.x, xa_add), std::max(p0_shift.y, yb_add), std::max(p0_shift.z, zc_add) };
		}
		// Get the limits of the bounding box of the semi-ellipsoid
		else if (this->isSemi) {
			p_bbox_min = { xa_sub, yb_sub, std::min(p0_shift.z, zc_add) };
			p_bbox_max = { xa_add, yb_add, std::max(p0_shift.z, zc_add) };
		}
		// Get the limit of the bounding box of the full ellispoid
		else {
			p_bbox_min = { xa_sub, yb_sub, zc_sub };
			p_bbox_max = { xa_add, yb_add, zc_add };
		}


		std::vector<double> sols;

		// The first interception point with the full ellipsoid is in the eight ellipsoid
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc1) && pc1.z >= 0.0)
			sols.push_back(solc1);

		// The second interception point with the full ellipsoid is in the eight ellipsoid
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc2) && pc2.z >= 0.0)
			sols.push_back(solc2);

		// The ray intersects the center of the target cell
		if (this->isInBBox(p_bbox_min, p_bbox_max, ptarget) && this->isInEllipsoid(p0_shift, ptarget))
			sols.push_back(0.0);

		// The ray intersects the Z horizontal plane of the eight ellipsoid
		if (this->isSemi) {
			if (this->isInBBox(p_bbox_min, p_bbox_max, pz0) && this->isInEllipsoid(p0_shift, pz0) && pz0.z > 0.0)
				sols.push_back(solz0);
		}

		// Add interception with the verticla planes of the 8th ellipsoid
		if (this->is8th) {
			// The ray intersects the X vertical plane of the eight ellipsoid
			if (this->isInBBox(p_bbox_min, p_bbox_max, px0) && this->isInEllipsoid(p0_shift, px0) && px0.z > 0.0 && !this->areEquals(px0, pz0))
				sols.push_back(solx0);

			// The ray intersects the Y vertical plane of the eight ellipsoid
			if (this->isInBBox(p_bbox_min, p_bbox_max, py0) && this->isInEllipsoid(p0_shift, py0) && py0.z > 0.0 && !this->areEquals(py0, pz0) && !this->areEquals(py0, px0))
				sols.push_back(soly0);
		}


		// FIND PATH LENGTH AND DISTANCE TO CELL CENTER
		int nsols = sols.size();

		// intersected the ellipse outside the half quarter
		if (nsols == 0) { return(nullptr); }

		// If not 0 or 2 solutions ==> problem
		if (sols.size() != 2) {
			std::cout << "Ellipsoid - Not 0 or 2 solutions - " << nsols << std::endl;
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
		interception* interc = new interception { this->vectIdTree, length, distance, false };

		return(interc);
	}

};

class CrownPartParaboloid : public CrownPart {

private:
	// Do we consider 4th of paraboloid or the full volume
	bool is4th;

	// Paraboloid parameters (semi-principal axes)
	// Sign of a and b determine which 4th of paraboloid we consider
	double a;
	double b;
	double h;


private:
	bool isInParaboloid(vertex3D p0, vertex3D p) {

		double value = (p.x - p0.x) * (p.x - p0.x) / (this->a * this->a) +
			(p.y - p0.y) * (p.y - p0.y) / (this->b * this->b) +  // Here, + not - because inverted paraboloid (i.e. paraboloid with apex at the top, not U-saphed paraboloid)
			(p.z - p0.z) / this->h; 

		return value < 0.0;
	}


public:
	CrownPartParaboloid(int vectid_tree,
		bool is_4th,
		double x, double y, double z,
		double a, double b, double h) :
		CrownPart(vectid_tree, x, y, z) {

		this->is4th = is_4th;

		this->a = a;
		this->b = b;
		this->h = h;
	}

	// interception of a 4th or full paraboloid between a shifted tree and a ray coming toward a target cell
	interception* computeInterception(int id_target_cell, Ray* ray, vertex3D& shift) {

		// SHIFT POSITION OF THE PARABOLOID
		// 2 shifts: the first one to set the target cell as origine, and the second to acoount for outside tree when torus system
		vertex3D p0_shift = {
			this->x + shift.x,
			this->y + shift.y,
			this->z + shift.z + this->h // Be careful, the reference point of the paraboloid is not the center of the base but the top point of the paraboloid
		};

		// Get the Z-coordinate of the base plane of the shifted paraboloid
		double zbase = p0_shift.z - this->h;
		

		// FIND SOLUTION OF THE QUADRATIC EQUATION (A * x * x + B * x + C = 0)
		// Equation giving distance between ray intersection with (*tree) crown and target cell center
		// Compute a, b and c coefficients
		double coef_a = ray->getCosHeightAngle() * ray->getCosHeightAngle() * ray->getCosAzimuth() * ray->getCosAzimuth() / (this->a * this->a) +
			ray->getCosHeightAngle() * ray->getCosHeightAngle() * ray->getSinAzimuth() * ray->getSinAzimuth() / (this->b * this->b);

		double coef_b = -2.0 * p0_shift.x * ray->getCosHeightAngle() * ray->getCosAzimuth() / (this->a * this->a) -
			2.0 * p0_shift.y * ray->getSinAzimuth() * ray->getCosHeightAngle() / (this->b * this->b) + 
			ray->getSinHeightAngle() / this->h;

		double coef_c = (p0_shift.x * p0_shift.x) / (this->a * this->a) +
			(p0_shift.y * p0_shift.y) / (this->b * this->b) - 
			p0_shift.z / this->h;


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

			// Coordinates of the two interception points with the full paraboloid
		vertex3D pc1 = getInterceptionPoint(solc1, ray);
		vertex3D pc2 = getInterceptionPoint(solc2, ray);

			// Interception with the horizontal base of the crown 
		double solbase = zbase / ray->getSinHeightAngle();
		vertex3D pbase = getInterceptionPoint(solbase, ray);
		pbase.z = zbase;


			// If interception with a fourth of paraboloid, consider interception with the X and Y vertical planes
			// Checked only in the case of a fourth of paraboloid
		double solx0 = p0_shift.x / (ray->getCosHeightAngle() * ray->getCosAzimuth());
		double soly0 = p0_shift.y / (ray->getCosHeightAngle() * ray->getSinAzimuth());

		vertex3D px0 = getInterceptionPoint(solx0, ray);
		vertex3D py0 = getInterceptionPoint(soly0, ray);

		px0.x = p0_shift.x;
		py0.y = p0_shift.y;


			// interception with the target point
		vertex3D ptarget = getInterceptionPoint(0.0, ray);


		// FIND SOLUTIONS OF INTERSECTION WITH THE EIGHTH ELLIPSOID
		double xa_sub = p0_shift.x - this->a;
		double yb_sub = p0_shift.y - this->b;

		double xa_add = p0_shift.x + this->a;
		double yb_add = p0_shift.y + this->b;

		// Get the limits of the bounding box of the 4th paraboloid
		vertex3D p_bbox_min = { 0.0, 0.0, 0.0 };
		vertex3D p_bbox_max = { 0.0, 0.0, 0.0 };
		if (this->is4th) {
			p_bbox_min = { std::min(p0_shift.x, xa_add), std::min(p0_shift.y, yb_add), zbase };
			p_bbox_max = { std::max(p0_shift.x, xa_add), std::max(p0_shift.y, yb_add), p0_shift.z };
		}
		// Get the limit of the bounding box of the full ellispoid
		else {
			p_bbox_min = { xa_sub, yb_sub, zbase };
			p_bbox_max = { xa_add, yb_add, p0_shift.z };
		}

		std::vector<double> sols;

		// The first interception point with the the infinite paraboloid is in the paraboloidal crown
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc1) && pc1.z >= 0.0)
			sols.push_back(solc1);

		// The second interception point with the infinite paraboloid is in the paraboloidal crown
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc2) && pc2.z >= 0.0)
			sols.push_back(solc2);

		// The ray intersects the Z horizontal base of the paraboloid
		if (this->isInBBox(p_bbox_min, p_bbox_max, pbase) && this->isInParaboloid(p0_shift, pbase) && pbase.z > 0.0)
			sols.push_back(solbase);

		// The ray intersects the center of the target cell
		if (this->isInBBox(p_bbox_min, p_bbox_max, ptarget) && this->isInParaboloid(p0_shift, ptarget))
			sols.push_back(0.0);

		// Add interception with the verticla planes of the 4th paraboloid
		if (this->is4th) {
			// The ray intersects the X vertical plane of the fourth paraboloid
			if (this->isInBBox(p_bbox_min, p_bbox_max, px0) && this->isInParaboloid(p0_shift, px0) && px0.z > 0.0 && !this->areEquals(px0, pbase))
				sols.push_back(solx0);

			// The ray intersects the Y vertical plane of the fourth paraboloid
			if (this->isInBBox(p_bbox_min, p_bbox_max, py0) && this->isInParaboloid(p0_shift, py0) && py0.z > 0.0 && !this->areEquals(py0, pbase) && !this->areEquals(py0, px0))
				sols.push_back(soly0);
		}


		// FIND PATH LENGTH AND DISTANCE TO CELL CENTER
		int nsols = sols.size();

		// intersected the parabol outside the half quarter
		if (nsols == 0) { return(nullptr); }

		// If not 0 or 2 solutions ==> problem
		if (sols.size() != 2) {
			std::cout << "Paraboloid crown - Not 0 or 2 solutions - " << nsols << std::endl;
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
		interception* interc = new interception{ this->vectIdTree, length, distance, false };

		return(interc);
	}

};

class Crown {

private:

	// Id of the tree to which the crown belong
	int vectIdTree;

	// Crown parts in the crown
	std::vector<CrownPart*> crownParts;

	// Leaves charcteristics
	double crownOpeness; // between 0 and 1, percentage of ligth energy of a ray that is not absorbed by the crown (for porous envelop method)
	double crownLAD; // Leaf area density in m2/m3 (for turbid medium method)
	
	// Output energy (in MJ)
	double energy; // Intercepted energy
	double energyPotential; // Intercepted energy without considering neighbours


private:

	// --- FUNCTIONS TO INITIALIZE DIFFERENT CROWN FORMS ---

	// Symetric paraboloid crown
	void initCrownParaboloid(double x, double y, double z,
		double h, double hbase,
		double cr_n, double cr_e, double cr_s, double cr_w) {

		// Radius is the same for all directions and is defined by mean crown radius of the four directions
		double crown_radius = (cr_n + cr_e + cr_s + cr_w) / 4.0;

		// Compute height of the crown
		double cdepth = (h - hbase);

		// Get center of the crown 
		// x, y, z are tree coordinates
		double x0 = x;
		double y0 = y;
		double z0 = z + hbase;

		// Create full paraboloid crown
		this->crownParts.push_back(new CrownPartParaboloid(this->vectIdTree, false,
			x0, y0, z0, crown_radius, crown_radius, cdepth));
	}

	// Irregular crown composed of 4 eighth of paraboloids
	void initCrownParaboloid4th(double x, double y, double z,
		double h, double hbase,
		double cr_n, double cr_e, double cr_s, double cr_w) {

		// Compute height of the crown
		double cdepth = h - hbase;

		// Get center of the crown 
		// x, y, z are tree coordinates
		double x0 = x;
		double y0 = y;
		double z0 = z + hbase;

		// Init fourth parts of the crown
		this->crownParts.push_back(new CrownPartParaboloid(this->vectIdTree, true,
			x0, y0, z0, cr_e, cr_n, cdepth));
		this->crownParts.push_back(new CrownPartParaboloid(this->vectIdTree, true,
			x0, y0, z0, cr_e, -cr_s, cdepth));
		this->crownParts.push_back(new CrownPartParaboloid(this->vectIdTree, true,
			x0, y0, z0, -cr_w, -cr_s, cdepth));
		this->crownParts.push_back(new CrownPartParaboloid(this->vectIdTree, true,
			x0, y0, z0, -cr_w, cr_n, cdepth));
	}

	// Symetric ellispoidal crown
	void initCrownEllipsoid(double x, double y, double z,
		double h, double hbase,
		double cr_n, double cr_e, double cr_s, double cr_w) {

		// Radius is the same for all directions and is defined by mean crown radius of the four directions
		double crown_radius = (cr_n + cr_e + cr_s + cr_w) / 4.0;

		// Compute height of maximum crown radius as the middle height of the crown depth
		double cdepth = (h - hbase);
		double hmax = hbase + cdepth / 2.0;

		// Get center of the crown 
		// x, y, z are tree coordinates
		double x0 = x;
		double y0 = y;
		double z0 = z + hmax;

		// Create full ellipsoid crown
		this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, false, false,
			x0, y0, z0, crown_radius, crown_radius, cdepth / 2.0));
	}

	// Irregular crown composed of 2 above and below semi-ellipsoids
	void initCrownEllipsoidSemi(double x, double y, double z,
		double h, double hbase, double hmax,
		double cr_n, double cr_e, double cr_s, double cr_w) {

		// Radius is the same for all directions and is defined by mean crown radius of the four directions
		double crown_radius = (cr_n + cr_e + cr_s + cr_w) / 4.0;

		// Get center of the crown 
		// x, y, z are tree coordinates
		double x0 = x;
		double y0 = y;
		double z0 = z + hmax;

		// Get depth of respectively top and bottom parts of the crown
		double depth_up = h - hmax;
		double depth_down = hmax - hbase;

		// Init top semi part of the crown
		if (depth_up > 0)
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, false,
				x0, y0, z0, crown_radius, crown_radius, depth_up));

		// Init bottom semi part of the crown
		if (depth_down > 0)
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, false,
				x0, y0, z0, crown_radius, crown_radius, -depth_down));
	}

	// Irregular crown composed of 8 eighth of ellipsoids
	void initCrownEllipsoid8th(double x, double y, double z,
		double h, double hbase, double hmax,
		double cr_n, double cr_e, double cr_s, double cr_w) {

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
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, cr_e, cr_n, depth_up));
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, cr_e, -cr_s, depth_up));
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, -cr_w, -cr_s, depth_up));
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, -cr_w, cr_n, depth_up));
		}

		// Init bottom eighth parts of the crown
		if (depth_down > 0) {
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, cr_e, cr_n, -depth_down));
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, cr_e, -cr_s, -depth_down));
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, -cr_w, -cr_s, -depth_down));
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, -cr_w, cr_n, -depth_down));
		}
	}


public:
	Crown(int vectid_tree, 
		std::string crown_type,
		double x, double y, double z,
		double h, double hbase, double hmax,
		double cr_n, double cr_e, double cr_s, double cr_w,
		double crown_openess, double crown_lad)
	{
		this->vectIdTree = vectid_tree;

		this->crownOpeness = crown_openess;
		this->crownLAD = crown_lad;

		this->energy = 0.0;
		this->energyPotential = 0.0;

		// Init the crown parts depending on the type
		if (crown_type == "8E")
			this->initCrownEllipsoid8th(x, y, z,
				h, hbase, hmax, 
				cr_n, cr_e, cr_s, cr_w);
		else if (crown_type == "2E") 
			this->initCrownEllipsoidSemi(x, y, z,
				h, hbase, hmax, 
				cr_n, cr_e, cr_s, cr_w);
		else if (crown_type == "E")
			this->initCrownEllipsoid(x, y, z,
				h, hbase, 
				cr_n, cr_e, cr_s, cr_w);
		else if (crown_type == "4P")
			this->initCrownParaboloid4th(x, y, z,
				h, hbase,
				cr_n, cr_e, cr_s, cr_w);
		else if (crown_type == "P") 
			this->initCrownParaboloid(x, y, z,
				h, hbase, 
				cr_n, cr_e, cr_s, cr_w);
		else
			std::cout << "Unrecognized crown type" << std::endl;
	}

	~Crown() {
		int n_parts = this->crownParts.size();
		for (int i = 0; i < n_parts; i++) {
			delete this->crownParts[i];
		}
	}

	// Getters
	int getVectIdTree() { return(this->vectIdTree); }
	int getNParts() { return(this->crownParts.size()); }
	double getCrownOpeness() { return(this->crownOpeness); }
	double getCrownLAD() { return(this->crownLAD); }
	double getEnergy() { return(this->energy); }
	double getEnergyPotential() { return(this->energyPotential); }

	// Setter
	void addEnergy(double e) { this->energy += e; }
	void addEnergyPotential(double epot) { this->energyPotential += epot; }


	// Find interception with a given part of the crown
	interception* computeCrownpartInterception(int id_crownpart, int id_target_cell, Ray* ray, vertex3D& shift) {
		return( this->crownParts[id_crownpart]->computeInterception(id_target_cell, ray, shift) );
	}
};

class Trunk : public TreeVolume {

private:

	// Position of the trunk (in meters)
	double x;
	double y;
	double z;

	// Dimensions of the tree
	double height;
	double radius;

private:

	bool isInCylinder(vertex3D p0, vertex3D p) {
		return (p0.x - p.x) * (p0.x - p.x) + (p0.y - p.y) * (p0.y - p.y) < this->radius * this->radius;
	}

		

public:
	Trunk(double vectid_tree, double x, double y, double z,
		double height, double radius) :
		TreeVolume(vectid_tree, x, y, z)
	{
		this->height = height;
		this->radius = radius;
	}

	// Compute interception betwene a cylinder and a ray coming toward a target cell
	interception* computeInterception(int id_target_cell, Ray* ray, vertex3D& shift) {

		// SHIFT POSITION OF THE CYLINDER
		// 2 shifts: the first one to set the target cell as origin, and the second to acoount for outside tree when torus system
		vertex3D p0_shift = {
			this->x + shift.x,
			this->y + shift.y,
			this->z + shift.z
		};
		// Correct for rounding errors after with equality between trunk base shifted center and target cell center
		if (abs(p0_shift.z) < EPSILON)
			p0_shift.z = 0.0;

		// FIND SOLUTION OF THE QUADRATIC EQUATION (A * x * x + B * x + C = 0)
		// Equation giving distance between ray intersection with infinite height cylinder and target cell center
		// Compute a, b and c coefficients
		double coef_a = ray->getCosHeightAngle() * ray->getCosHeightAngle();
		double coef_b = -2.0 * p0_shift.x * ray->getCosHeightAngle() * ray->getCosAzimuth() - 
			2.0 * p0_shift.y * ray->getCosHeightAngle() * ray->getSinAzimuth();
		double coef_c = p0_shift.x * p0_shift.x + p0_shift.y * p0_shift.y - this->radius * this->radius;

		// Find point of intersection between the ray and the whole cylinder
		// Find positive solutions of the above quadratic equations (distance to target cell)
		// if not 2 solutions, leave the function
		// Because negative = no interception and null = 1 interception = do not consider tangent rays with crown
		double delta = coef_b * coef_b - 4.0 * coef_a * coef_c;

		if (delta <= 0.0) { return nullptr; }

		double delta_sqrt = sqrt(delta);
		double solc1 = (-coef_b + delta_sqrt) / (2.0 * coef_a);
		double solc2 = (-coef_b - delta_sqrt) / (2.0 * coef_a);


		// COMPUTE COORDINATES OF INTERCEPTION POINTS WITH TRUNK LIMITS AND TARGET POINT

			// Coordinates of the two interception points with the full ellipsoid
		vertex3D pc1 = this->getInterceptionPoint(solc1, ray);
		vertex3D pc2 = this->getInterceptionPoint(solc2, ray);

		// Interception with the top horizontal base of the trunk  
		double ztop = p0_shift.z + this->height;
		double solztop = ztop / ray->getSinHeightAngle();
		vertex3D pztop = this->getInterceptionPoint(solztop, ray);
		pztop.z = ztop;

		double solzbot = p0_shift.z / ray->getSinHeightAngle();
		vertex3D pzbot = this->getInterceptionPoint(solzbot, ray);
		pzbot.z = p0_shift.z;

		// interception with the target point
		vertex3D ptarget = getInterceptionPoint(0.0, ray);


		// FIND SOLUTIONS OF INTERSECTION WITH THE EIGHTH ELLIPSOID
		// Get the limits of the bounding box of the 8th paraboloid
		vertex3D p_bbox_min = { p0_shift.x - this->radius, p0_shift.y - this->radius, p0_shift.z };
		vertex3D p_bbox_max = { p0_shift.x + this->radius, p0_shift.y + this->radius, ztop };

		std::vector<double> sols;

		// The first interception point with the full ellipsoid is in the trunk
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc1) && pc1.z >= 0.0 && pc1.z >= p0_shift.z) 
			sols.push_back(solc1);

		// The second interception point with the full ellipsoid is in the trunk
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc2) && pc2.z >= 0.0 && pc2.z >= p0_shift.z)
			sols.push_back(solc2);

		// The ray intersects the ground or top horizontal plane of the trunk
		if (this->isInCylinder(p0_shift, pztop) && ztop >= 0.0 && ztop >= p0_shift.z)
			sols.push_back(solztop);
		if (this->isInCylinder(p0_shift, pzbot) && p0_shift.z >= 0.0)
			sols.push_back(solzbot);

		// The ray intersects the center of the target cell
		if (this->isInBBox(p_bbox_min, p_bbox_max, ptarget) && this->isInCylinder(p0_shift, ptarget) && ptarget.z != pzbot.z && ptarget.z != pztop.z)
			sols.push_back(0.0);


		// FIND PATH LENGTH AND DISTANCE TO CELL CENTER
		int nsols = sols.size();

		// Did not intersect the trunk
		if (nsols == 0) { return(nullptr); }

		// If not 0 or 2 solutions ==> problem
		if (sols.size() != 2) {
			std::cout << "Trunk - Not 0 or 2 solutions - " << nsols << std::endl;
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
		interception* interc = new interception{ this->vectIdTree, length, distance, true };

		return(interc);
	}
};

class Tree {

private:
	// Id of the tree in the trees vector (in Stand object)
	// And id of the tree in the input dataset
	int vectId;
	int id;

	// Tree volumes: trunk and crown
	Crown crown;
	Trunk trunk;


public:
	Tree(int vectid, int id, double x, double y, double z,
		double dbh, double h,
		std::string crown_type,
		double hbase, double hmax,
		double cr_n, double cr_e, double cr_s, double cr_w,
		double crown_openess, double crown_lad):

		crown(vectid, crown_type, x, y, z, h, hbase, hmax, cr_n, cr_e, cr_s, cr_w,
			crown_openess, crown_lad),
		trunk(vectid, x, y, z, dbh/200.0, h)
	{
		this->vectId = vectid;
		this->id = id;
	}


	// Getters
	int getId() { return(this->id); }
	int getVectId() { return(this->vectId); }
	Crown& getCrown() { return(this->crown); }
	double getCrownEnergy() { return(this->crown.getEnergy()); }
	double getCrownEnergyPotential() { return(this->crown.getEnergyPotential()); }

	// Interception methods
	interception* computeCrownpartInterception(int id_crownpart, int id_target_cell, Ray* ray, vertex3D& shift) {
		return(this->crown.computeCrownpartInterception(id_crownpart, id_target_cell, ray, shift));
	}

	interception* computeTrunkInterception(int id_target_cell, Ray* ray, vertex3D& shift) {
		return( this->trunk.computeInterception(id_target_cell, ray, shift) );
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

	// Id of the tree in the cell (id in the sense of id in the trees vector in Stand object)
	std::vector<int> vectIdTrees;

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
		this->energy = e_above; // Current energy above the cell (will be decreased by successive interception by the abov trees)
	}

	// Getters
	int getId() { return(this->id); }
	int getRow() { return(this->row); }
	int getCol() { return(this->col); }
	double getX() { return(this->x); }
	double getY() { return(this->y); }
	double getZ() { return(this->z); }

	bool isEmpty() { return(this->vectIdTrees.empty()); }
	int getVectIdTree(int i_tree) { return(this->vectIdTrees[i_tree]); }
	int getNTrees() { return(this->vectIdTrees.size()); }

	double getEnergy() { return(this->energy); }
	double getEnergyM2(double cell_area) { return(this->energy / cell_area); }
	double getEnergyRelative() { return(this->energy / this->energyAboveCanopy); }


	// Setters
	void addTree(int vect_id) { this->vectIdTrees.push_back(vect_id); }
	void interceptEnergy(double e) { this->energy -= e; }
	void correctNullEnergy(double epsilon) {
		if (this->energy < -epsilon)
			std::cout << "Problem with energy in cell " << this->id << ": energy of " << this->energy << "MJ";
		if (abs(this->energy) < epsilon)
			this->energy = 0.0;
	}

};

class Stand {

private:
	// 2D grid of cells stored by <row, column>
	std::vector<std::vector<Cell*>> grid;

	// Vector of trees
	std::vector<Tree*> trees;

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
	int nCellsX; // Number of cells in a row
	int nCellsY; // Number of cells in a column
	int nCells; // Total number of cells
	double cellAreaPlane; // Area of the cell at horizontal (in m2)
	double cellArea; // Area of the cell considering slope (in m2)
	bool useTorus; // Whether we want to use a torus system for light ray tracing


public:
	Stand(DataFrame trees, double e_above_m2,
		double slope, double north_to_x_cw, double aspect,
		double cell_size, double n_cells_x, double n_cells_y, 
		bool use_torus
		)
	{
		// Stand orientation (all in radians)
		this->slope = slope * M_PI / 180.0;
		this->northToXcw = north_to_x_cw * M_PI / 180.0;
		this->aspect = aspect * M_PI / 180.0;
		this->bottomAzimuth = (-aspect + north_to_x_cw) * M_PI / 180.0;


		// Stand size
		this->cellSize = cell_size;

		this->nCellsX = n_cells_x;
		this->nCellsY = n_cells_y;
		this->nCells = n_cells_x * n_cells_y;

		this->cellAreaPlane = cell_size * cell_size;
		this->cellArea = this->cellAreaPlane / cos(this->slope);
		this->useTorus = use_torus;


		// Get energy above canopy that is coming toward a cell
		double e_above_cell = e_above_m2 * this->cellArea;


		// Create grid of cells
		for (int r = 0; r < this->nCellsY; r++) {
			std::vector<Cell*> cells_in_row;
			for (int c = 0; c < this->nCellsX; c++) {
				
				// Find position of cell center
				double x = this->cellSize * (c + 1.0 / 2.0);
				double y = this->cellSize * (this->nCellsX - r - 1.0 / 2.0); // Be careful, grid is first rows on the top but Y coordinate is greater for first rows
				
				// Create cell object
				cells_in_row.push_back(new Cell(
					this->nCellsY*r + c + 1, 
					r, c, x, y, this->computeZ(x, y),
					e_above_cell
				));
			}
			this->grid.push_back(cells_in_row);
		}


		// Get vectors from the trees R dataframe
		IntegerVector trees_id = trees["id_tree"];
		NumericVector trees_x = trees["x"];
		NumericVector trees_y = trees["y"];
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

		// Add a Tree object to the associated Cell they belong to
		int n_trees = trees.nrows();
		for (int i = 0; i < n_trees; i++) {

			// Create tree object
			this->trees.push_back(new Tree(
				i, trees_id[i], trees_x[i], trees_y[i], this->computeZ(trees_x[i], trees_y[i]),
				dbh[i], height[i], Rcpp::as< std::string >(ctype[i]), hbase[i], hmax[i],
				cr_n[i], cr_e[i], cr_s[i], cr_w[i], cp[i], clad[i]));

			// Find row and column of the cell and add tree id to the corresponding cell
			int cell_col = (int)(trees_x[i] / this->cellSize);
			int cell_row = this->nCellsY - (int)(trees_y[i] / this->cellSize) - 1;
			this->grid[cell_row][cell_col]->addTree(i);
		}

		// Get maximum height and crown radius of the trees
		this->maxTreeHeight = max(height);
		this->maxTreeCrownRadius = std::max({ max(cr_n), max(cr_e), max(cr_s), max(cr_w) });
	}

	~Stand() {
		// Delete cells pointers
		for (int r = 0; r < this->nCellsY; r++) {
			for (int c = 0; c < this->nCellsX; c++) {
				delete this->grid[r][c];
			}
		}
		// Delete tree pointers
		int n_trees = this->trees.size();
		for (int i = 0; i < n_trees; i++) {
			delete this->trees[i];
		}
	}

	// Getters
	Cell* getCell(int row, int col) { return(this->grid[row][col]); }
	int getNCells() { return(this->nCells); }
	int getNCellsX() { return(this->nCellsX); }
	int getNCellsY() { return(this->nCellsY); }

	Tree* getTree(int vectid) { return(this->trees[vectid]); }
	int getNTrees() { return(this->trees.size()); }

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


	// Methods
	double computeZ(double x, double y) {

		double d = sqrt(x*x + y*y);
		double azimuth_xy = 0;
		if (d != 0) {
			if (y >= 0) 
				azimuth_xy = acos(x / d);
			else
				azimuth_xy = 2.0 * M_PI - acos(x / d);
		}
		return(-d * cos(azimuth_xy - this->bottomAzimuth) * tan(this->slope));
	}

};

struct shiftedCell {

	// Pointor toward the original Cell
	Cell* cell;

	// Shift applied to the cell position (in a torus system)
	vertex3D shift;
};

class Model {

private:
	// Stand object with cells, trees and geometric variables
	Stand stand;

	// Rays that are coming toward a cell
	RayManager rays;

	// Light interception method by the crowns : either porous envelop(use "crownOpenness") or turbid medium(use "crownLAD")
	bool turbidMedium;

	// Consider interception with trunks ?
	bool trunkInterception;

	// Global variable for turbid medium
	const double EXTINCTION_COEF; // Probability of a leaf to intercept the ray (linked to leaf orientation)
	const double CLUMPING_FACTOR; // Aggregation of leaves within the crown volume (1 is homogeneous)


private:

	shiftedCell getPotentialCell(Ray* ray, int rel_cell_id, Cell* target_cell) {

		relativeCoords rel_cell = ray->getPotentialRelCell(rel_cell_id);

		// Get row and column of potential cell
		int row_pot = target_cell->getRow() + rel_cell.row;
		int col_pot = target_cell->getCol() + rel_cell.col;

		// Search if cells is outside the this->stand
		bool is_outside = row_pot < 0 || row_pot >= this->stand.getNCellsY() || col_pot < 0 || col_pot >= this->stand.getNCellsX();

		// If cel is outside the main this->stand, do not keep as potential if we are not in a torus system
		if (is_outside && !this->stand.isUseTorus()) { return(shiftedCell{ nullptr, 0.0, 0.0, 0.0 }); }

		// Get row and column of the original cell in the grid cell (diffrent only within a torus system, if cell is outside)
		int row_pot_original = (this->stand.getNCellsY() + (row_pot % this->stand.getNCellsY())) % this->stand.getNCellsY();
		int col_pot_original = (this->stand.getNCellsX() + (col_pot % this->stand.getNCellsX())) % this->stand.getNCellsX();

		// Get corresponding original cell
		Cell* original_cell = this->stand.getCell(row_pot_original, col_pot_original);

		// Check if it contains trees, otherwise, do not add as potential cell
		if ((*original_cell).isEmpty()) { return(shiftedCell{ nullptr, 0.0, 0.0, 0.0 }); }

		// Compute shift in coordinates to apply on trees within the original cell (torus system)
		// Shift is a multiple of stand_size : how many this->stand size we have to shift tree coordinates times this->stand size, to apply torus system
		// If we have a potential cell with negative or gretaer than n_cells row or column
		// ==> Thus, we want to substract coordinates of original cells in order to have potential cell outside the main this->stand
		// Be careful, y - coordinates system is opposite direction of y - grid system
		double x_shift = (col_pot / this->stand.getNCellsX() - (col_pot % this->stand.getNCellsX() < 0 ? 1 : 0)) * (this->stand.getCellSize() * this->stand.getNCellsX()); // Negative x_id ==> negative shift
		double y_shift = -(row_pot / this->stand.getNCellsY() - (row_pot % this->stand.getNCellsY() < 0 ? 1 : 0)) * (this->stand.getCellSize() * this->stand.getNCellsY()); // Negative y_id ==> positive shift

		// Compute altitudinal shift (only if (x,y) shifted
		double z_shift = 0.0;
		if (x_shift != 0.0 || y_shift != 0.0)
			z_shift = this->stand.computeZ(x_shift, y_shift);
		
		// Return the shifted cell
		vertex3D shift = { x_shift, y_shift, z_shift };
		shiftedCell pot_cell = { original_cell, shift };

		return(pot_cell);
	}

	std::vector<interception*> computeInterceptionsRayToTarget(Ray* ray, Cell* target_cell) {

		std::vector<interception*> v_interc;

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
			vertex3D shift = {
				pot_cell.shift.x - target_cell->getX(),
				pot_cell.shift.y - target_cell->getY(),
				pot_cell.shift.z - target_cell->getZ()
			};

			// Compute interception for each tree within the potential cell (with a given shift if the potential cell is outside the plot)
			int n_pot_trees = pot_cell.cell->getNTrees();
			for (int t = 0; t < n_pot_trees; t++) {

				// Get intercepted tree
				Tree* tree = this->stand.getTree(pot_cell.cell->getVectIdTree(t));
				
				// Compute interception for each crown part of the tree
				int n_parts = tree->getCrown().getNParts();
				for (int p = 0; p < n_parts; p++) {
					interception* interc_crownpart = tree->computeCrownpartInterception(p, target_cell->getId(), ray, shift);
					if (interc_crownpart != nullptr) {
						v_interc.push_back(interc_crownpart);
					}
				}

				// If specified, compute interception by trunks
				if (this->trunkInterception) {
					interception* interc_trunk = tree->computeTrunkInterception(target_cell->getId(), ray, shift);
					if (interc_trunk != nullptr) {
						v_interc.push_back(interc_trunk);
					}
				}
			}
		}

		return(v_interc);
	}

	void orderInterceptionsRayToTarget(std::vector<interception*>& v_interc) {
		std::sort(v_interc.begin(), v_interc.end(), [](interception const* i1, interception const* i2) {
			// Sort by distance to the target cell center, furthest first
			return (i1->distance > i2->distance);
			});
	}

	void summarizeInterceptionsRayToTarget(Ray* ray, Cell* target_cell, std::vector<interception*>& v_interc) {

		// Compute projection of energy on plane parallel to slope
		//  in MJ / m2 and convert it into MJ per cell
		double scalar_slope = cos(stand.getSlope()) * ray->getSinHeightAngle() +
			sin(stand.getSlope()) * ray->getCosHeightAngle() * cos(ray->getAzimuth() - stand.getBottomAzimuth());

		double e_incident_slope_m2 = scalar_slope * ray->getIncidentEnergy();
		double e_incident_slope_cell = e_incident_slope_m2 * stand.getCellArea();

		// Initialize the energy of the ray coming toward the cell above the canopy
		double current_energy = e_incident_slope_cell;
		bool trunk_intercepted = false;

		// Initialize tracking maps to check if a crown has already been intercepted (used when porous envelop)
		// and to track the incident potential energy that arrive to a crown part (used for potential energy in turbid medium when many crown parts)
		std::vector<bool> is_intercepted(this->stand.getNTrees(), false);
		std::vector<double> e_pot_incident(this->stand.getNTrees(), e_incident_slope_cell);

		//Compute attenuation of energy across successives crown interceptions for each ray X target cell
		int n_interceptions = v_interc.size();
		for (int j = 0; j < n_interceptions; j++) {

			// Intrception with a trunk
			if (v_interc[j]->withTrunk) {
				// Set the current energy to 0
				current_energy = 0.0;
				continue;
			}

			// Get the intercepted crown
			Crown& crown = stand.getTree(v_interc[j]->vectIdTree)->getCrown();

			// Compute the potential energy to the tree (energy without attenuation by neighbours)
			// and compute energy with attenuation
			// Different considering crown as turbid medium or porous envelop
			double potential_energy = 0.0;
			double intercepted_energy = 0.0;

			// Turbid medium ==> apply beer lambert law
			if (this->turbidMedium) {

				// Apply Beer Lambert law on incdent energy attenuated by previous crown parts of the tree, without considering crowns of other competitors
				potential_energy = this->applyBeerLambert(
					e_pot_incident[v_interc[j]->vectIdTree],
					this->EXTINCTION_COEF, this->CLUMPING_FACTOR,
					crown.getCrownLAD(),
					v_interc[j]->length);

				// Remove the intercepted energy from the energy that left and from the potential incident energy of the given tree crown
				e_pot_incident[v_interc[j]->vectIdTree] = e_pot_incident[v_interc[j]->vectIdTree] - potential_energy;

				// Apply beer lambert law on incident ray attenuated by the competitor crowns intercepted the crown
				intercepted_energy = trunk_intercepted ? 0.0 : this->applyBeerLambert(
					current_energy,
					this->EXTINCTION_COEF, this->CLUMPING_FACTOR,
					crown.getCrownLAD(),
					v_interc[j]->length);
				
			}
			// Porous envelop ==> reduce the energy by a fixed amount
			else {
				// Be careful, only once per crown (if many crown part, reduce only if crown has not been already intercepted by the ray)
				if (!is_intercepted[v_interc[j]->vectIdTree]) {

					// Compute potential and intercepted energy if the first criwn part is encountered by reducing the incident enrrgy by a fixed amount
					potential_energy = e_incident_slope_cell * (1 - crown.getCrownOpeness());
					intercepted_energy = trunk_intercepted ? 0.0 : (current_energy * (1 - crown.getCrownOpeness()));

					// Set the crown as being intercepted
					is_intercepted[v_interc[j]->vectIdTree] = true;
				}
				else {
					potential_energy = 0.0;
					intercepted_energy = 0.0;
				}
			}

			// Add to the potential and intercepted energy by the tree
			#ifdef _OPENMP
			#pragma omp critical
			{
			#endif
			crown.addEnergyPotential(potential_energy);
			crown.addEnergy(intercepted_energy);
			#ifdef _OPENMP
			}
			#endif

			// Remove the intercepted energy from the energy that left
			current_energy -= intercepted_energy;

			// Remove the intercepted energy by the crown from the total energy above the cell
			#ifdef _OPENMP
			#pragma omp critical
			#endif
			target_cell->interceptEnergy(intercepted_energy);

			// Delete interception pointer
			delete v_interc[j];
		}
	}


	double applyBeerLambert(double incident_energy, double extinction_coef, double clumping_factor, double leaf_area_density, double path_length) {

		return(incident_energy * (1 - exp(-extinction_coef * clumping_factor * leaf_area_density * path_length)));
	}


public:

	Model(DataFrame trees, DataFrame rays,
		double e_above_m2,
		double slope, double north_to_x_cw, double aspect,
		double cell_size, double n_cells_x, double n_cells_y,
		bool use_torus, bool turbid_medium, bool trunk_interception) :

		stand(trees, e_above_m2, slope, north_to_x_cw, aspect, cell_size, n_cells_x, n_cells_y, use_torus),
		rays(rays, e_above_m2),
		EXTINCTION_COEF(0.5), CLUMPING_FACTOR(1.0)
	{
		this->turbidMedium = turbid_medium;
		this->trunkInterception = trunk_interception;
	}

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

	void computeInterceptions() {

		int n_rays = this->rays.getNRays();
		int n_cells_x = this->stand.getNCellsX();
		int n_cells_y = this->stand.getNCellsY();

		// For each ray toward each target cell
		#ifdef _OPENMP
		#pragma omp parallel for collapse(3)
		#endif
		for (int r = 0; r < n_cells_y; r++) {
			for (int c = 0; c < n_cells_x; c++) {
				for (int i = 0; i < n_rays; i++) {

					Cell* target_cell = this->stand.getCell(r, c);
					Ray* ray = this->rays.getRay(i);

					// Get interceptions
					std::vector<interception*> v_interc = this->computeInterceptionsRayToTarget(ray, target_cell);

					// Order interceptions
					this->orderInterceptionsRayToTarget(v_interc);

					// Summarize interceptions into energy
					this->summarizeInterceptionsRayToTarget(ray, target_cell, v_interc);
				}
			}
		}
	}

	List exportResults() {

		// Get number of trees and cells
		int n_cells = this->stand.getNCells();
		int n_trees = this->stand.getNTrees();

		// Init RCPP vectors for trees
		IntegerVector id_tree(n_trees);
		NumericVector e_tree(n_trees);
		NumericVector epot_tree(n_trees);
		NumericVector lci_tree(n_trees);

		// Init RCPP vectors for cells
		IntegerVector id_cell(n_cells);
		NumericVector x_cell(n_cells);
		NumericVector y_cell(n_cells);
		NumericVector z_cell(n_cells);
		NumericVector e_cell(n_cells);
		NumericVector erel_cell(n_cells);

		// For each cell
		int icell = 0;
		int itree = 0;
		int n_rows = this->stand.getNCellsY();
		int n_cols = this->stand.getNCellsX();
		for (int r = 0; r < n_rows; r++) {
			for (int c = 0; c < n_cols; c++) {

				// Get cell
				Cell* cell = this->stand.getCell(r, c);

				// Add cell to vectors
				id_cell[icell] = cell->getId();
				x_cell[icell] = cell->getX();
				y_cell[icell] = cell->getY();
				z_cell[icell] = cell->getZ();

				// Correct for rounding errors
				cell->correctNullEnergy(1e-6);

				e_cell[icell] = cell->getEnergyM2(this->stand.getCellArea());
				erel_cell[icell] = cell->getEnergyRelative();

				// For each tree composing the cell
				int n_trees_cell = cell->getNTrees();
				for (int t = 0; t < n_trees_cell; t++) {

					Tree* tree = this->stand.getTree(cell->getVectIdTree(t));

					id_tree[itree] = tree->getId();
					e_tree[itree] = tree->getCrownEnergy();
					epot_tree[itree] = tree->getCrownEnergyPotential();

					itree++;
				}

				icell++;
			}
		}

		// Create trees and cells RCPP DataFrames
		DataFrame output_trees = DataFrame::create(
			Named("id_tree") = id_tree,
			Named("epot") = epot_tree,
			Named("e") = e_tree
		);

		DataFrame output_cells = DataFrame::create(
			Named("id_cell") = id_cell,
			Named("x_center") = x_cell,
			Named("y_center") = y_cell,
			Named("z_center") = z_cell,
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
	DataFrame trees, DataFrame rays, double e_above_m2,
	double slope, double north_to_x_cw, double aspect,
	double cell_size, double n_cells_x, double n_cells_y,
	bool use_torus, bool turbid_medium, bool trunk_interception
)
{
	// TESTS:
	// - If plots is big enough for tree coordinates
	// - If hmax is needed

	// TODO:
	// - core plot


	// Initialize the model
	Model sl_model = Model(trees,
		rays, e_above_m2,
		slope, north_to_x_cw, aspect,
		cell_size, n_cells_x, n_cells_y, 
		use_torus, turbid_medium, trunk_interception);

	// [OPTIMIZATION]: find the extend of possible interception of each ray
	// i.e. for each ray coming to an unknown target cell, find relative cells around the target containing tree that could possibly intercept the ray
	// We can whether or not consider a torus system
	sl_model.findPotentialRelCells();

	// For each target cell and each ray (target cells can be parallelized)
	sl_model.computeInterceptions();

	// Convert teh output into R format
	List output = sl_model.exportResults();

	return(output);
}


