
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

	// Is the ray a direct of diffuse one ?
	bool direct;

	// Output interceptions vector
	std::vector<interception*> interceptions;


public:
	Ray(int id, double energy, double height_angle, double azimuth, bool direct) {
		this->id = id;
		this->incidentEnergy = energy;
		this->heightAngle = height_angle;
		this->cosHeightAngle = cos(height_angle);
		this->sinHeightAngle = sin(height_angle);
		this->azimuth = azimuth;
		this->cosAzimuth = cos(azimuth);
		this->sinAzimuth = sin(azimuth);
		this->direct = direct;
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
	bool isDirect() { return(this->direct); }

	void addPotentialRelCell(int row, int col) { this->potentialRelCells.push_back(relativeCoords{ row, col }); }
	relativeCoords getPotentialRelCell(int id) { return(this->potentialRelCells[id]); }
	int getNPotCells() { return(this->potentialRelCells.size()); }

};

class RayManager {

private:
	// Vector of rays
	std::vector<Ray*> rays;

	// Energy coming from diffuse and direct rays above canopy (in MJ/m2)
	// In either a horizontal plane or the slope
	double energyDirectAboveSlopeM2;
	double energyDiffuseAboveSlopeM2;

	double energyDirectAboveHorizontalM2;
	double energyDiffuseAboveHorizontalM2;


public:
	RayManager(DataFrame rays, 
		double e_direct_above_slope_m2, double e_diffuse_above_slope_m2,
		double e_direct_above_horizontal_m2, double e_diffuse_above_horizontal_m2) {

		// Energy in MJ per m2 above canopy, direct/diffuse and on a slope or on a horizontal plane
		this->energyDirectAboveSlopeM2 = e_direct_above_slope_m2;
		this->energyDiffuseAboveSlopeM2 = e_diffuse_above_slope_m2;

		this->energyDirectAboveHorizontalM2 = e_direct_above_horizontal_m2;
		this->energyDiffuseAboveHorizontalM2 = e_diffuse_above_horizontal_m2;

		// Number of rays
		int n_rays = rays.nrows();

		// Get vectors from the R dataframe
		IntegerVector id = rays["id_ray"];
		NumericVector energy = rays["e_incident"];
		NumericVector height_angle = rays["height_angle"];
		NumericVector azimuth = rays["azimut"];
		LogicalVector direct = rays["direct"];

		// Create a Tree object and store it in a vector of trees
		for (int i = 0; i < n_rays; i++) {
			Ray* ray = new Ray(id[i]-1, energy[i], height_angle[i], azimuth[i], direct[i]);
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

	double getEnergyAboveSlopeM2() { return(this->energyDirectAboveSlopeM2 + this->energyDiffuseAboveSlopeM2); }
	double getEnergyDirectAboveSlopeM2() { return(this->energyDirectAboveSlopeM2); }
	double getEnergyDiffuseAboveSlopeM2() { return(this->energyDiffuseAboveSlopeM2); }

	double getEnergyAboveHorizontalM2() { return(this->energyDirectAboveHorizontalM2 + this->energyDiffuseAboveHorizontalM2); }
	double getEnergyDirectAboveHorizontalM2() { return(this->energyDirectAboveHorizontalM2); }
	double getEnergyDiffuseAboveHorizontalM2() { return(this->energyDiffuseAboveHorizontalM2); }
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
	virtual interception* computeInterception(Ray* ray, vertex3D& shift) { return(nullptr); }

};

class CrownPart: public TreeVolume {

public:
	CrownPart(int vectid_tree,
		double x, double y, double z):
		TreeVolume(vectid_tree, x, y, z)
	{}

	virtual ~CrownPart() {}

	// interception generic function
	virtual interception* computeInterception(Ray* ray, vertex3D& shift) { return(nullptr); }
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
	interception* computeInterception(Ray* ray, vertex3D& shift) {

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
	interception* computeInterception(Ray* ray, vertex3D& shift) {

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
	double energyDiffuse; // Same for diffuse energy
	double energyPotentialDiffuse;

	double energyDirect; // Same for direct energy
	double energyPotentialDirect;



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

		this->energyDiffuse = 0.0;
		this->energyPotentialDiffuse = 0.0;

		this->energyDirect = 0.0;
		this->energyPotentialDirect = 0.0;


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

	double getEnergy() { return(this->energyDiffuse + this->energyDirect); }
	double getEnergyPotential() { return(this->energyPotentialDiffuse + this->energyPotentialDirect); }

	double getEnergyDiffuse() { return(this->energyDiffuse); }
	double getEnergyPotentialDiffuse() { return(this->energyPotentialDiffuse); }

	double getEnergyDirect() { return(this->energyDirect); }
	double getEnergyPotentialDirect() { return(this->energyPotentialDirect); }


	// Setter
	void addEnergyDiffuse(double e) { this->energyDiffuse += e; }
	void addEnergyPotentialDiffuse(double epot) { this->energyPotentialDiffuse += epot; }

	void addEnergyDirect(double e) { this->energyDirect += e; }
	void addEnergyPotentialDirect(double epot) { this->energyPotentialDirect += epot; }


	// Find interception with a given part of the crown
	interception* computeCrownpartInterception(int id_crownpart, Ray* ray, vertex3D& shift) {
		return( this->crownParts[id_crownpart]->computeInterception(ray, shift) );
	}
};

class Trunk : public TreeVolume {

private:

	// Position of the trunk (in meters)
	double x;
	double y;
	double z;

	// Dimensions of the tree (in meters)
	double height;
	double radius;

private:

	bool isInCylinder(vertex3D p0, vertex3D p) {
		return (p0.x - p.x) * (p0.x - p.x) + (p0.y - p.y) * (p0.y - p.y) < this->radius * this->radius;
	}

		

public:
	Trunk(double vectid_tree, double x, double y, double z,
		double height_m, double dbh_cm) :
		TreeVolume(vectid_tree, x, y, z)
	{
		this->height = height_m;
		this->radius = dbh_cm/200.0; // because radius is in meters
	}

	// Compute interception betwene a cylinder and a ray coming toward a target cell
	interception* computeInterception(Ray* ray, vertex3D& shift) {

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
		trunk(vectid, x, y, z, dbh, h)
	{
		this->vectId = vectid;
		this->id = id;
	}


	// Getters
	int getId() { return(this->id); }
	int getVectId() { return(this->vectId); }

	Trunk& getTrunk() { return(this->trunk); }
	Crown& getCrown() { return(this->crown); }

	double getCrownEnergy() { return(this->crown.getEnergy()); }
	double getCrownEnergyPotential() { return(this->crown.getEnergyPotential()); }

	double getCrownEnergyDirect() { return(this->crown.getEnergyDirect()); }
	double getCrownEnergyPotentialDirect() { return(this->crown.getEnergyPotentialDirect()); }

	double getCrownEnergyDiffuse() { return(this->crown.getEnergyDiffuse()); }
	double getCrownEnergyPotentialDiffuse() { return(this->crown.getEnergyPotentialDiffuse()); }


	// Methods
	interception* computeCrownpartInterception(int id_crownpart, Ray* ray, vertex3D& shift) {
		return(this->crown.computeCrownpartInterception(id_crownpart, ray, shift));
	}

	interception* computeTrunkInterception(Ray* ray, vertex3D& shift) {
		return( this->trunk.computeInterception(ray, shift) );
	}


};

class Target {

private:

	// Point towards which to cast the rays
	// For sensor, it is its position
	// For cell, it is its center
	double x;
	double y;
	double z;

	// Column and row in a grid system (i.e. Stand object)
	// In the case of the sensor, column and row of the cell the sensor belong to 
	int row;
	int col;

	// Energy of the target (in MJ)
	double energyDirectSlope; // Current direct energy on the target on the slope
	double energyDiffuseSlope; // Current diffuse energy on the target on the slope
	double energyDirectHorizontal; // Current direct energy on the target on a horizontal plane
	double energyDiffuseHorizontal; // Current diffuse energy on the target on a horizontal plane

	// If the target is a sensor (other case is a cell)
	// In this case, do not remove energies from intercepted crown
	bool isSensor;


public:
	Target(double x, double y, double z, 
		int row, int col, 
		double e_direct_above_slope, double e_diffuse_above_slope,
		double e_direct_above_horizontal, double e_diffuse_above_horizontal,
		bool is_sensor)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->row = row;
		this->col = col;

		// will be decreased by successive interception by the above trees
		this->energyDirectSlope = e_direct_above_slope;
		this->energyDiffuseSlope = e_diffuse_above_slope;
		this->energyDirectHorizontal = e_direct_above_horizontal;
		this->energyDiffuseHorizontal = e_diffuse_above_horizontal;

		this->isSensor = is_sensor;
	}

	virtual ~Target() {}


	// Getters
	double getX() { return(this->x); }
	double getY() { return(this->y); }
	double getZ() { return(this->z); }

	int getRow() { return(this->row); }
	int getCol() { return(this->col); }

	double getEnergySlope() { return(this->energyDirectSlope + this->energyDiffuseSlope); }
	double getEnergyDirectSlope() { return(this->energyDirectSlope); }
	double getEnergyDiffuseSlope() { return(this->energyDiffuseSlope); }

	double getEnergyHorizontal() { return(this->energyDirectHorizontal + this->energyDiffuseHorizontal); }
	double getEnergyDirectHorizontal() { return(this->energyDirectHorizontal); }
	double getEnergyDiffuseHorizontal() { return(this->energyDiffuseHorizontal); }

	bool isThisSensor() { return(this->isSensor); }


	// Setters
	void interceptEnergyDirectSlope(double e) { this->energyDirectSlope -= e; }
	void interceptEnergyDiffuseSlope(double e) { this->energyDiffuseSlope -= e; }
	void interceptEnergyDirectHorizontal(double e) { this->energyDirectHorizontal -= e; }
	void interceptEnergyDiffuseHorizontal(double e) { this->energyDiffuseHorizontal -= e; }


	//void resetEnergy() { this->energy = 0.0; }

	//virtual void correctNullEnergy(double epsilon) {}
};

class Sensor : public Target {

private:
	// Id of the sensor
	int idSensor;


public:
	Sensor(double x, double y, double z, 
		int row, int col,
		double e_direct_above_slope_m2, double e_diffuse_above_slope_m2,
		double e_direct_above_horizontal_m2, double e_diffuse_above_horizontal_m2,
		int id_sensor) :
		Target(x, y, z, row, col,
			e_direct_above_slope_m2, e_diffuse_above_slope_m2, 
			e_direct_above_horizontal_m2, e_diffuse_above_horizontal_m2,
			true)
	{
		this->idSensor = id_sensor;
	}

	// Getters
	int getIdSensor() { return(this->idSensor); }

	// Methods
	/*void correctNullEnergy(double epsilon) {
		if (this->getEnergy() < -epsilon)
			std::cout << "Problem with energy in target " << this->getIdSensor() << ": energy of " << this->getEnergy() << "MJ";
		if (abs(this->getEnergy()) < epsilon)
			this->resetEnergy();
	}*/
};

class Cell : public Target {

private:
	// Id of the cell
	int idCell;

	// Id of the tree in the cell (id in the sense of id in the trees vector in Stand object)
	std::vector<int> vectIdTrees;


public:
	Cell(double x, double y, double z, 
		int row, int col, 
		double e_direct_above_slope_m2, double e_diffuse_above_slope_m2,
		double e_direct_above_horizontal_m2, double e_diffuse_above_horizontal_m2,
		int id_cell) :
		Target(x, y, z, row, col, 
			e_direct_above_slope_m2, e_diffuse_above_slope_m2, 
			e_direct_above_horizontal_m2, e_diffuse_above_horizontal_m2, 
			false)
	{
		this->idCell = id_cell;
	}

	// Getters
	int getIdCell() { return(this->idCell); }

	bool isEmpty() { return(this->vectIdTrees.empty()); }
	int getVectIdTree(int i_tree) { return(this->vectIdTrees[i_tree]); }
	int getNTrees() { return(this->vectIdTrees.size()); }


	// Setters
	void addTree(int vect_id) { this->vectIdTrees.push_back(vect_id); }

	// Methods
	//void correctNullEnergy(double epsilon) {
	//	if (this->getEnergy() < -epsilon)
	//		std::cout << "Problem with energy in cell " << this->getIdCell() << ": energy of " << this->getEnergy() << "MJ";
	//	if (abs(this->getEnergy()) < epsilon)
	//		this->resetEnergy();
	//}
};

class Stand {

private:
	// Vector of sensors
	bool areSensors;
	std::vector<Sensor*> sensors;

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
	double cellAreaHorizontal; // Area of the cell at horizontal (in m2)
	double cellAreaSlope; // Area of the cell considering slope (in m2)
	bool useTorus; // Whether we want to use a torus system for light ray tracing


public:
	Stand(DataFrame trees, DataFrame sensors,
		double e_direct_above_slope_m2, double e_diffuse_above_slope_m2,
		double e_direct_above_horizontal_m2, double e_diffuse_above_horizontal_m2,
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

		this->cellAreaHorizontal = cell_size * cell_size;
		this->cellAreaSlope = this->cellAreaHorizontal / cos(this->slope);
		this->useTorus = use_torus;


		// Get energy above canopy that is coming toward a cell
		double e_direct_above_slope_cell = e_direct_above_slope_m2 * this->cellAreaSlope;
		double e_diffuse_above_slope_cell = e_diffuse_above_slope_m2 * this->cellAreaSlope;
		double e_direct_above_horizontal_cell = e_direct_above_horizontal_m2 * this->cellAreaHorizontal;
		double e_diffuse_above_horizontal_cell = e_diffuse_above_horizontal_m2 * this->cellAreaHorizontal;


		// Create vector of sensors from the sensors R DataFrame, if sensors are specified
		this->areSensors = sensors.nrows() > 0;

		if (this->areSensors) {

			IntegerVector sensors_id = sensors["id_sensor"];
			NumericVector sensors_x = sensors["x"];
			NumericVector sensors_y = sensors["y"];
			NumericVector sensors_height = sensors["h_m"];

			int n_sensors = sensors.nrows();
			for (int i = 0; i < n_sensors; i++) {

				// Check if the sensor is within the stand limits
				if (sensors_x[i] < 0.0 ||
					sensors_x[i] > this->cellSize * this->nCellsX ||
					sensors_y[i] < 0.0 ||
					sensors_y[i] > this->cellSize * this->nCellsY) {

					// Print an error (NEED TO STOP ?)
					std::cout << "Sensor " << sensors_id[i] << " is outside the stand limits" << std::endl;

				}
				else {
					// Find the row and column the sensor belong to
					int sensor_col = findColFromPosX(sensors_x[i]);
					int sensor_row = findRowFromPosY(sensors_y[i]);

					// std::cout << "Sensor " << sensors_id[i] << " in (row, col): (" << sensor_row << "," << sensor_col << ")" << std::endl;

					// Create and push the pointor to Sensor new instance
					this->sensors.push_back(new Sensor(
						sensors_x[i], sensors_y[i], // Position of the sensor
						this->computeZ(sensors_x[i], sensors_y[i]) + sensors_height[i], // Z pos of the sensor is the height above the ground
						sensor_row, sensor_col, // Cell of the sensor
						e_direct_above_slope_cell, e_diffuse_above_slope_cell, // Direct and diffuse energy coming from above the canopy on the slope
						e_direct_above_horizontal_cell, e_diffuse_above_horizontal_cell, // Same but on a horizontal plane
						sensors_id[i] // Id of the sensor
					));
				}
			}

		}


		// Create grid of cells
		for (int r = 0; r < this->nCellsY; r++) {
			std::vector<Cell*> cells_in_row;
			for (int c = 0; c < this->nCellsX; c++) {
				
				// Find position of cell center
				double x = this->cellSize * (c + 1.0 / 2.0);
				double y = this->cellSize * (this->nCellsY - r - 1.0 / 2.0); // Be careful, grid is first rows on the top but Y coordinate is greater for first rows
				
				// Create cell object
				cells_in_row.push_back(new Cell(
					x, y, this->computeZ(x, y), // Position of the cell center
					r, c, // Row and column of the cell within the stand grid system
					e_direct_above_slope_cell, e_diffuse_above_slope_cell, // Direct and diffuse energy coming from above the canopy on the slope
					e_direct_above_horizontal_cell, e_diffuse_above_horizontal_cell, // Same but on a horizontal plane
					this->nCellsY * r + c + 1 // Id of the cell
				));
			}
			this->grid.push_back(cells_in_row);
		}


		// Create vector of trees from the trees R dataframe
		IntegerVector trees_id = trees["id_tree"];
		NumericVector trees_x = trees["x"];
		NumericVector trees_y = trees["y"];
		NumericVector dbh = trees["dbh_cm"];
		NumericVector tree_height = trees["h_m"];
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
				dbh[i], tree_height[i], Rcpp::as< std::string >(ctype[i]), hbase[i], hmax[i],
				cr_n[i], cr_e[i], cr_s[i], cr_w[i], cp[i], clad[i]));

			// Find row and column of the cell and add tree id to the corresponding cell
			int tree_col = findColFromPosX(trees_x[i]);
			int tree_row = findRowFromPosY(trees_y[i]);
			this->grid[tree_row][tree_col]->addTree(i);

			// std::cout << "Tree " << trees_id[i] << " in (row, col): (" << tree_row << "," << tree_col << ")" << std::endl;

		}

		// Get maximum height and crown radius of the trees
		this->maxTreeHeight = max(tree_height);
		this->maxTreeCrownRadius = std::max({ max(cr_n), max(cr_e), max(cr_s), max(cr_w) });
	}

	~Stand() {
		// Delete sensors pointors
		int n_sensors = this->sensors.size();
		for (int i = 0; i < n_sensors; i++) {
			delete this->sensors[i];
		}

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
	bool areThereSensors() { return(this->areSensors); }
	Sensor* getSensor(int vectid) { return(this->sensors[vectid]); }
	double getNSensors() { return(this->sensors.size()); }

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
	double getCellAreaHorizontal() { return(this->cellAreaHorizontal); }
	double getCellAreaSlope() { return(this->cellAreaSlope); }
	double isUseTorus() { return(this->useTorus); }
	double getMaximumTreeCrownRadius() { return(this->maxTreeCrownRadius); }
	double getMaximumTreeHeight() { return(this->maxTreeHeight); }


	// Methods
	double computeZ(double x, double y) {

		double d = sqrt(x*x + y*y);
		double azimuth_xy = 0.0;
		if (d != 0.0) {
			if (y >= 0.0) 
				azimuth_xy = acos(x / d);
			else
				azimuth_xy = 2.0 * M_PI - acos(x / d);
		}
		return(-d * cos(azimuth_xy - this->bottomAzimuth) * tan(this->slope));
	}

	int findColFromPosX(double x) { return((int)(x / this->cellSize)); }
	int findRowFromPosY(double y) { return(this->nCellsY - (int)(y / this->cellSize) - 1); }


};

struct shiftedCell {

	// Pointor toward the original Cell
	Cell* cell;

	// Shift applied to the cell position (in a torus system)
	vertex3D shift;
};

class Model {

private:
	// Stand object with cells, trees, sensors and geometric variables
	Stand stand;

	// Rays that are coming toward a cell
	RayManager rays;

	// Light interception method by the crowns : either porous envelop(use "crownOpenness") or turbid medium(use "crownLAD")
	bool turbidMedium;

	// Consider interception with trunks ?
	bool trunkInterception;

	// Compute interceptions only of sensors ?
	// i.e. do not cast rays towards cells and do not predict tree energy interception
	bool sensorsOnly;

	// Global variable for turbid medium
	const double EXTINCTION_COEF; // Probability of a leaf to intercept the ray (linked to leaf orientation)
	const double CLUMPING_FACTOR; // Aggregation of leaves within the crown volume (1 is homogeneous)


private:

	shiftedCell getPotentialCell(Ray* ray, int rel_cell_id, Target* target) {

		relativeCoords rel_cell = ray->getPotentialRelCell(rel_cell_id);

		// Get row and column of potential cell
		int row_pot = target->getRow() + rel_cell.row;
		int col_pot = target->getCol() + rel_cell.col;

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

	void predictEnergiesRayToTarget(Ray* ray, Target* target) {

		// Get interceptions
		std::vector<interception*> v_interc = this->computeInterceptionsRayToTarget(ray, target);

		// Order interceptions
		this->orderInterceptionsRayToTarget(v_interc);

		// Summarize interceptions into energy
		this->summarizeInterceptionsRayToTarget(ray, target, v_interc);

	}

	std::vector<interception*> computeInterceptionsRayToTarget(Ray* ray, Target* target) {

		std::vector<interception*> v_interc;

		// For each possible relative cell coordinates
		int n_relcells = ray->getNPotCells();
		for (int rc = 0; rc < n_relcells; rc++) {

			// Get a pointer to the original potential cell and its associated shift
			shiftedCell pot_cell = this->getPotentialCell(ray, rc, target);

			// Cell is not potential if there is no trees or it is outside the main plot and torus system is disabled
			if (!pot_cell.cell) { continue; }

			// Compute the shift applied to the tree (2 shifts)
			// 1. Shift when torus system is enable and the tree is outside the main plot
			// 2. Set the center of the target cell as origin of tree and ray interception
			vertex3D shift = {
				pot_cell.shift.x - target->getX(),
				pot_cell.shift.y - target->getY(),
				pot_cell.shift.z - target->getZ()
			};

			// Compute interception for each tree within the potential cell (with a given shift if the potential cell is outside the plot)
			int n_pot_trees = pot_cell.cell->getNTrees();
			for (int t = 0; t < n_pot_trees; t++) {

				// Get intercepted tree
				Tree* tree = this->stand.getTree(pot_cell.cell->getVectIdTree(t));
				
				// Compute interception for each crown part of the tree
				int n_parts = tree->getCrown().getNParts();
				for (int p = 0; p < n_parts; p++) {
					interception* interc_crownpart = tree->computeCrownpartInterception(p, ray, shift);
					if (interc_crownpart != nullptr) {
						v_interc.push_back(interc_crownpart);
					}
				}

				// If specified, compute interception by trunks
				if (this->trunkInterception) {
					interception* interc_trunk = tree->computeTrunkInterception(ray, shift);
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

	void summarizeInterceptionsRayToTarget(Ray* ray, Target* target, std::vector<interception*>& v_interc) {

		// Compute projection of energy on plane parallel to slope
		double scalar_slope = cos(stand.getSlope()) * ray->getSinHeightAngle() +
			sin(stand.getSlope()) * ray->getCosHeightAngle() * cos(ray->getAzimuth() - stand.getBottomAzimuth());

		double e_incident_slope_m2 = scalar_slope * ray->getIncidentEnergy();


		// Compute projection of energy on a horizontal plane
		double scalar_horizontal = cos(0.0) * ray->getSinHeightAngle() +
			sin(0.0) * ray->getCosHeightAngle() * cos(ray->getAzimuth() - stand.getBottomAzimuth());

		double e_incident_horizontal_m2 = scalar_horizontal * ray->getIncidentEnergy();


		// And convert incident ray energy from MJ / m2 into MJ within the whole target
		double e_incident_slope_cell = e_incident_slope_m2 * stand.getCellAreaSlope();
		double e_incident_horizontal_cell = e_incident_horizontal_m2 * stand.getCellAreaHorizontal();


		// Initialize the energy of the ray coming toward the cell above the canopy
		// Intercepted energy by the trees is from the rays with initial energy considering the slope
		// But track the attenuation of energy toward a horizontal plane
		// For tracking the enrgy arriving at at the target both at horizontal and on a slope (cell or sensor)
		double current_energy_slope = e_incident_slope_cell;
		double current_energy_horizontal = e_incident_horizontal_cell;

		// Track if the ray has been intercepted by a trunk
		bool trunk_intercepted = false;

		// Initialize tracking maps necessary when many crown parts, i.e. 4P/4E/8E
		// Check if a crown has already been intercepted (in porous envelop: as it does not consider path across the crown, attenuate the ray only at the first crown part intercepted)
		// Track the incident potential energy that arrive to a crown part (for potential energy in turbid medium)
		std::vector<bool> is_intercepted(this->stand.getNTrees(), false);
		std::vector<double> e_pot_incident(this->stand.getNTrees(), e_incident_slope_cell);

		//Compute attenuation of energy across successives crown interceptions for each ray X target cell
		int n_interceptions = v_interc.size();
		for (int j = 0; j < n_interceptions; j++) {

			// If the ray has already intercepted a trunk
			// Do not run other interception
			// But still loop on interceptions to delete the pointors
			if (trunk_intercepted) {
				delete v_interc[j];
				continue;
			}

			// First interception with a trunk
			// Next step, set current energy to 0 by setting intercepted energy to the incident energy
			if (v_interc[j]->withTrunk) {
				trunk_intercepted = true;
			}

			// Get the intercepted crown
			Crown& crown = stand.getTree(v_interc[j]->vectIdTree)->getCrown();

			// Compute the potential energy to the tree (energy without attenuation by neighbours)
			// and compute energy with attenuation
			// Different considering crown as turbid medium or porous envelop
			// Consider intercepted energy both on a slope or on a horizontal plane
			// But the intercepted energy by trees is always energy on a slope
			double potential_energy_slope = 0.0;
			double intercepted_energy_slope = 0.0;
			double intercepted_energy_horizontal = 0.0;

			// Turbid medium ==> apply beer lambert law
			if (this->turbidMedium) {

				// Apply Beer Lambert law on incdent energy attenuated by previous crown parts of the tree, without considering crowns of other competitors
				potential_energy_slope = this->applyBeerLambert(
					e_pot_incident[v_interc[j]->vectIdTree],
					this->EXTINCTION_COEF, this->CLUMPING_FACTOR,
					crown.getCrownLAD(),
					v_interc[j]->length);

				// Remove the intercepted energy from the energy that left and from the potential incident energy of the given tree crown
				e_pot_incident[v_interc[j]->vectIdTree] = e_pot_incident[v_interc[j]->vectIdTree] - potential_energy_slope;

				// Apply beer lambert law on incident ray attenuated by the competitor crowns intercepted the crown
				intercepted_energy_slope = trunk_intercepted ? 
					current_energy_slope : 
					this->applyBeerLambert(
						current_energy_slope,
						this->EXTINCTION_COEF, this->CLUMPING_FACTOR,
						crown.getCrownLAD(),
						v_interc[j]->length);

				intercepted_energy_horizontal = trunk_intercepted ? 
					current_energy_horizontal : 
					this->applyBeerLambert(
						current_energy_horizontal,
						this->EXTINCTION_COEF, this->CLUMPING_FACTOR,
						crown.getCrownLAD(),
						v_interc[j]->length);
				
			}
			// Porous envelop ==> reduce the energy by a fixed amount
			else {
				// Be careful, only once per crown (if many crown part, reduce only if crown has not been already intercepted by the ray)
				// Except the intercepted part is the trunk (to avoid missing trunk interception after a previous interception by another volume of the same tree)
				if (!is_intercepted[v_interc[j]->vectIdTree] || v_interc[j]->withTrunk) {

					// Compute potential and intercepted energy if the first crown part is encountered by reducing the incident enrrgy by a fixed amount
					potential_energy_slope = e_incident_slope_cell * (1 - crown.getCrownOpeness());

					intercepted_energy_slope = trunk_intercepted ? 
						current_energy_slope : 
						(current_energy_slope * (1 - crown.getCrownOpeness()));

					intercepted_energy_horizontal = trunk_intercepted ? 
						current_energy_horizontal : 
						(current_energy_slope * (1 - crown.getCrownOpeness()));

					// Set the crown as being intercepted
					is_intercepted[v_interc[j]->vectIdTree] = true;
				}
				else {
					// In this case, another crown part has been already intercepted
					// So avoid adding multiple attenuation of the same crown because multiple crown part have been intercepted
					potential_energy_slope = 0.0;
					intercepted_energy_slope = 0.0;
					intercepted_energy_horizontal = 0.0;
				}
			}

			// Add to the potential and intercepted energy by the tree
			// CAREFUL: ONLY IF TARGET IS A CELL 
			// If it is a sensor, do not consider energy of trees
			if (!target->isThisSensor()) {
				#ifdef _OPENMP
				#pragma omp critical
				{
				#endif

				// Add energy to wether the ray is a diffuse or direct one
				if (ray->isDirect()) {
					crown.addEnergyPotentialDirect(potential_energy_slope);
					crown.addEnergyDirect(intercepted_energy_slope);
				}
				else {
					crown.addEnergyPotentialDiffuse(potential_energy_slope);
					crown.addEnergyDiffuse(intercepted_energy_slope);
				}

				#ifdef _OPENMP
				}
				#endif
			}

			// Remove the intercepted energy from the energy of the ray
			// i.e. output energy is the transmitted energy
			current_energy_slope -= intercepted_energy_slope;
			current_energy_horizontal -= intercepted_energy_horizontal;

			// Remove the intercepted energy by the crown from the total energy above the target
			// (When considering energy  direct/diffuse,) 
			// (it is easier to decrease the energy coming from above along interceptions)
			// (than to add current energy of the ray coming to the target at the end of the interceptions)
			#ifdef _OPENMP
			#pragma omp critical
			{
			#endif

			// Diffuse or direct energy
			if (ray->isDirect()) {
				target->interceptEnergyDirectSlope(intercepted_energy_slope);
				target->interceptEnergyDirectHorizontal(intercepted_energy_horizontal);
			}
			else {
				target->interceptEnergyDiffuseSlope(intercepted_energy_slope);
				target->interceptEnergyDiffuseHorizontal(intercepted_energy_horizontal);
			}

			#ifdef _OPENMP
			}
			#endif

			// Delete interception pointer
			delete v_interc[j];
		}

	}

	double applyBeerLambert(double incident_energy, double extinction_coef, double clumping_factor, double leaf_area_density, double path_length) {

		return(incident_energy * (1 - exp(-extinction_coef * clumping_factor * leaf_area_density * path_length)));
	}


public:

	Model(DataFrame trees, 
		DataFrame sensors, bool sensors_only,
		DataFrame rays, 
		double e_direct_above_slope_m2, double e_diffuse_above_slope_m2,
		double e_direct_above_horizontal_m2, double e_diffuse_above_horizontal_m2,
		double slope, double north_to_x_cw, double aspect,
		double cell_size, double n_cells_x, double n_cells_y,
		bool use_torus, bool turbid_medium, bool trunk_interception) :

		stand(trees, sensors, 
			e_direct_above_slope_m2, e_diffuse_above_slope_m2, e_direct_above_horizontal_m2, e_diffuse_above_horizontal_m2,
			slope, north_to_x_cw, aspect, cell_size, n_cells_x, n_cells_y, use_torus),
		rays(rays, e_direct_above_slope_m2, e_diffuse_above_slope_m2, e_direct_above_horizontal_m2, e_diffuse_above_horizontal_m2),
		EXTINCTION_COEF(0.5), CLUMPING_FACTOR(1.0)
	{
		this->sensorsOnly = sensors_only;
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
		int n_sensors = this->stand.getNSensors();
		int n_cells_x = this->stand.getNCellsX();
		int n_cells_y = this->stand.getNCellsY();

		// For each ray toward each target sensor
		if (this->stand.areThereSensors()) {
			#ifdef _OPENMP
			#pragma omp parallel for collapse(2)
			#endif
			for (int s = 0; s < n_sensors; s++) {
				for (int i = 0; i < n_rays; i++) {

					Ray* ray = this->rays.getRay(i);
					Sensor* sensor = this->stand.getSensor(s);

					predictEnergiesRayToTarget(ray, sensor);
				}
			}
		}

		// For each ray toward each target cell
		// Run only in conditions
		if (!this->sensorsOnly) {
			#ifdef _OPENMP
			#pragma omp parallel for collapse(3)
			#endif
			for (int r = 0; r < n_cells_y; r++) {
				for (int c = 0; c < n_cells_x; c++) {
					for (int i = 0; i < n_rays; i++) {

						Ray* ray = this->rays.getRay(i);
						Cell* target_cell = this->stand.getCell(r, c);

						predictEnergiesRayToTarget(ray, target_cell);
					}
				}
			}
		}
	}

	List exportResults() {

		// ------------- FOR SENSORS ----------------

		// Init RCPP vectors for sensors
		int n_sensors = this->stand.getNSensors();

		IntegerVector id_sensor(n_sensors);
		NumericVector x_sensor(n_sensors);
		NumericVector y_sensor(n_sensors);
		NumericVector z_sensor(n_sensors);

		NumericVector e_slope_sensor(n_sensors);
		NumericVector pacl_slope_sensor(n_sensors);
		NumericVector e_slope_direct_sensor(n_sensors);
		NumericVector pacl_slope_direct_sensor(n_sensors);
		NumericVector e_slope_diffuse_sensor(n_sensors);
		NumericVector pacl_slope_diffuse_sensor(n_sensors);

		NumericVector e_horizontal_sensor(n_sensors);
		NumericVector pacl_horizontal_sensor(n_sensors);
		NumericVector e_horizontal_direct_sensor(n_sensors);
		NumericVector pacl_horizontal_direct_sensor(n_sensors);
		NumericVector e_horizontal_diffuse_sensor(n_sensors);
		NumericVector pacl_horizontal_diffuse_sensor(n_sensors);


		// Export the sensor into R
		for (int s = 0; s < n_sensors; s++) {

			// Get sensor
			Sensor* sensor = this->stand.getSensor(s);

			// Add sensor to vectors
			id_sensor[s] = sensor->getIdSensor();
			x_sensor[s] = sensor->getX();
			y_sensor[s] = sensor->getY();
			z_sensor[s] = sensor->getZ();

			// Correct for rounding errors
			/*sensor->correctNullEnergy(1e-6);*/

			e_slope_sensor[s] = sensor->getEnergySlope();
			e_slope_direct_sensor[s] = sensor->getEnergyDirectSlope();
			e_slope_diffuse_sensor[s] = sensor->getEnergyDiffuseSlope();

			pacl_slope_sensor[s] = sensor->getEnergySlope() / (this->rays.getEnergyAboveSlopeM2() * this->stand.getCellAreaSlope());
			pacl_slope_direct_sensor[s] = sensor->getEnergyDirectSlope() / (this->rays.getEnergyDirectAboveSlopeM2() * this->stand.getCellAreaSlope());
			pacl_slope_diffuse_sensor[s] = sensor->getEnergyDiffuseSlope() / (this->rays.getEnergyDiffuseAboveSlopeM2() * this->stand.getCellAreaSlope());
		
			e_horizontal_sensor[s] = sensor->getEnergyHorizontal();
			e_horizontal_direct_sensor[s] = sensor->getEnergyDirectHorizontal();
			e_horizontal_diffuse_sensor[s] = sensor->getEnergyDiffuseHorizontal();

			pacl_horizontal_sensor[s] = sensor->getEnergyHorizontal() / (this->rays.getEnergyAboveHorizontalM2() * this->stand.getCellAreaHorizontal());
			pacl_horizontal_direct_sensor[s] = sensor->getEnergyDirectHorizontal() / (this->rays.getEnergyDirectAboveHorizontalM2() * this->stand.getCellAreaHorizontal());
			pacl_horizontal_diffuse_sensor[s] = sensor->getEnergyDiffuseHorizontal() / (this->rays.getEnergyDiffuseAboveHorizontalM2() * this->stand.getCellAreaHorizontal());
		}

		// Create sensors RCPP DataFrames
		DataFrame output_sensors = DataFrame::create(
			Named("id_sensor") = id_sensor,
			Named("x") = x_sensor,
			Named("y") = y_sensor,
			Named("z") = z_sensor,

			Named("e_slope") = e_slope_sensor,
			Named("pacl_slope") = pacl_slope_sensor,
			Named("e_slope_direct") = e_slope_direct_sensor,
			Named("pacl_slope_direct") = pacl_slope_direct_sensor,
			Named("e_slope_diffuse") = e_slope_diffuse_sensor,
			Named("pacl_slope_diffuse") = pacl_slope_diffuse_sensor,

			Named("e_horizontal") = e_horizontal_sensor,
			Named("pacl_horizontal") = pacl_horizontal_sensor,
			Named("e_horizontal_direct") = e_horizontal_direct_sensor,
			Named("pacl_horizontal_direct") = pacl_horizontal_direct_sensor,
			Named("e_horizontal_diffuse") = e_horizontal_diffuse_sensor,
			Named("pacl_horizontal_diffuse") = pacl_horizontal_diffuse_sensor
		);


		// -------------- FOR TREES AND CELLS ------------------

		// Get number of trees and cells
		// Set to 0 if we did not compute trees and cells interception (sensor_only argument)
		int n_cells = this->stand.getNCells();
		int n_trees = this->stand.getNTrees();

		if (this->sensorsOnly) {
			n_cells = 0;
			n_trees = 0;
		}

		// Init RCPP vectors for trees
		IntegerVector id_tree(n_trees);
		NumericVector x_tree(n_trees);
		NumericVector y_tree(n_trees);
		NumericVector z_tree(n_trees);

		NumericVector e_tree(n_trees);
		NumericVector epot_tree(n_trees);
		NumericVector e_direct_tree(n_trees);
		NumericVector epot_direct_tree(n_trees);
		NumericVector e_diffuse_tree(n_trees);
		NumericVector epot_diffuse_tree(n_trees);

		// Init RCPP vectors for cells
		IntegerVector id_cell(n_cells);
		NumericVector x_cell(n_cells);
		NumericVector y_cell(n_cells);
		NumericVector z_cell(n_cells);

		NumericVector e_slope_cell(n_cells);
		NumericVector pacl_slope_cell(n_cells);
		NumericVector e_slope_direct_cell(n_cells);
		NumericVector pacl_slope_direct_cell(n_cells);
		NumericVector e_slope_diffuse_cell(n_cells);
		NumericVector pacl_slope_diffuse_cell(n_cells);

		NumericVector e_horizontal_cell(n_cells);
		NumericVector pacl_horizontal_cell(n_cells);
		NumericVector e_horizontal_direct_cell(n_cells);
		NumericVector pacl_horizontal_direct_cell(n_cells);
		NumericVector e_horizontal_diffuse_cell(n_cells);
		NumericVector pacl_horizontal_diffuse_cell(n_cells);

		if (!this->sensorsOnly) {

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
					id_cell[icell] = cell->getIdCell();
					x_cell[icell] = cell->getX();
					y_cell[icell] = cell->getY();
					z_cell[icell] = cell->getZ();

					// Correct for rounding errors
					//cell->correctNullEnergy(1e-6);

					e_slope_cell[icell] = cell->getEnergySlope();
					e_slope_direct_cell[icell] = cell->getEnergyDirectSlope();
					e_slope_diffuse_cell[icell] = cell->getEnergyDiffuseSlope();

					pacl_slope_cell[icell] = cell->getEnergySlope() / (this->rays.getEnergyAboveSlopeM2() * this->stand.getCellAreaSlope());
					pacl_slope_direct_cell[icell] = cell->getEnergyDirectSlope() / (this->rays.getEnergyDirectAboveSlopeM2() * this->stand.getCellAreaSlope());
					pacl_slope_diffuse_cell[icell] = cell->getEnergyDiffuseSlope() / (this->rays.getEnergyDiffuseAboveSlopeM2() * this->stand.getCellAreaSlope());

					e_horizontal_cell[icell] = cell->getEnergyHorizontal();
					e_horizontal_direct_cell[icell] = cell->getEnergyDirectHorizontal();
					e_horizontal_diffuse_cell[icell] = cell->getEnergyDiffuseHorizontal();

					pacl_horizontal_cell[icell] = cell->getEnergyHorizontal() / (this->rays.getEnergyAboveHorizontalM2() * this->stand.getCellAreaHorizontal());
					pacl_horizontal_direct_cell[icell] = cell->getEnergyDirectHorizontal() / (this->rays.getEnergyDirectAboveHorizontalM2() * this->stand.getCellAreaHorizontal());
					pacl_horizontal_diffuse_cell[icell] = cell->getEnergyDiffuseHorizontal() / (this->rays.getEnergyDiffuseAboveHorizontalM2() * this->stand.getCellAreaHorizontal());



					// For each tree composing the cell
					int n_trees_cell = cell->getNTrees();
					for (int t = 0; t < n_trees_cell; t++) {

						Tree* tree = this->stand.getTree(cell->getVectIdTree(t));

						id_tree[itree] = tree->getId();
						x_tree[itree] = tree->getTrunk().getX();
						y_tree[itree] = tree->getTrunk().getY();
						z_tree[itree] = tree->getTrunk().getZ();

						e_tree[itree] = tree->getCrownEnergy();
						epot_tree[itree] = tree->getCrownEnergyPotential();

						e_direct_tree[itree] = tree->getCrownEnergyDirect();
						epot_direct_tree[itree] = tree->getCrownEnergyPotentialDirect();

						e_diffuse_tree[itree] = tree->getCrownEnergyDiffuse();
						epot_diffuse_tree[itree] = tree->getCrownEnergyPotentialDiffuse();

						itree++;
					}

					icell++;
				}
			}

		}

		// Create trees and cells RCPP DataFrames
		DataFrame output_trees = DataFrame::create(
			Named("id_tree") = id_tree,
			Named("x") = x_tree,
			Named("y") = y_tree,
			Named("z") = z_tree,

			Named("epot") = epot_tree,
			Named("e") = e_tree,
			Named("epot_direct") = epot_direct_tree,
			Named("e_direct") = e_direct_tree,
			Named("epot_diffuse") = epot_diffuse_tree,
			Named("e_diffuse") = e_diffuse_tree
		);

		DataFrame output_cells = DataFrame::create(
			Named("id_cell") = id_cell,
			Named("x_center") = x_cell,
			Named("y_center") = y_cell,
			Named("z_center") = z_cell,

			Named("e_slope") = e_slope_cell,
			Named("pacl_slope") = pacl_slope_cell,
			Named("e_slope_direct") = e_slope_direct_cell,
			Named("pacl_slope_direct") = pacl_slope_direct_cell,
			Named("e_slope_diffuse") = e_slope_diffuse_cell,
			Named("pacl_slope_diffuse") = pacl_slope_diffuse_cell,

			Named("e_horizontal") = e_horizontal_cell,
			Named("pacl_horizontal") = pacl_horizontal_cell,
			Named("e_horizontal_direct") = e_horizontal_direct_cell,
			Named("pacl_horizontal_direct") = pacl_horizontal_direct_cell,
			Named("e_horizontal_diffuse") = e_horizontal_diffuse_cell,
			Named("pacl_horizontal_diffuse") = pacl_horizontal_diffuse_cell
		);

		// Return output as a List of two DataFrames
		return(List::create(
			Named("sensors") = output_sensors,
			Named("trees") = output_trees,
			Named("cells") = output_cells
		));
	}

};




// [[Rcpp::export]]
List sl_run_rcpp(
	DataFrame trees, 
	DataFrame sensors, bool sensors_only,
	DataFrame rays, 
	double e_direct_above_slope_m2, double e_diffuse_above_slope_m2,
	double e_direct_above_horizontal_m2, double e_diffuse_above_horizontal_m2,
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
	Model sl_model = Model(
		trees, 
		sensors, sensors_only,
		rays, 
		e_direct_above_slope_m2, e_diffuse_above_slope_m2,
		e_direct_above_horizontal_m2, e_diffuse_above_horizontal_m2,
		slope, north_to_x_cw, aspect,
		cell_size, n_cells_x, n_cells_y, 
		use_torus, turbid_medium, trunk_interception
	);

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


