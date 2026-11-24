
#define _USE_MATH_DEFINES

constexpr double EPSILON_TRUNK_BASE = 1e-9; // For rounding errors on the height of the trunk base relative to the target
constexpr double EPSILON_FRACTION = 1e-10; // For rounding errors in crown fractions (8E and 4P), as in Capsis (F. Andre 2019)
constexpr double EPSILON_TRUNK = 1e-6; // For rounding errors in trunk bounding box, as in Capsis (SLTrunk)

#include <vector>
#include <map>
#include <string> 
#include <cmath>
#include <limits>
#include <algorithm>
#include <atomic>

#include <Rcpp.h>
using namespace Rcpp;

#ifdef _OPENMP
#include <omp.h>
#endif

static bool g_parallel_enabled = false;
static int  g_num_threads = 1;
static bool g_threads_auto = true;

// Number of ray / crown part interceptions with an odd number of interception points that are not a known particular case
// Such interceptions are ignored (as in Capsis) and reported to the user in R (atomic because updated within OpenMP threads)
static std::atomic<long> g_n_unexpected_ellipsoid(0);
static std::atomic<long> g_n_unexpected_paraboloid(0);
static std::atomic<long> g_n_unexpected_trunk(0);


// [[Rcpp::plugins(openmp)]]


struct relativeCoords {
	int row = 0;
	int col = 0;
};

struct vertex3D {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
};

// Type of an interception point between a ray and a crown part (used to identify particular cases, as in Capsis)
enum solutionType { SOL_SURFACE, SOL_TARGET, SOL_PLANE };

struct interception {

	// Id of the intercepted tree in the trees vector in the Stand object
	int vectIdTree = -1;

	// Length of the path across the crown
	double length = 0.0;

	// Distance between interception point (middle of full path) and target cell
	double distance = 0.0;

	// Interception by crown or trunk ?
	bool withTrunk = false;
};

// Geometry of the rays toward a set of targets (sensors or cells), independent from the crown LAD and openness
// Pairs = ray x target pairs intercepted at least once (indices start at 1, as in R), in the order of the targets
// Crossings = crowns crossed by the ray of each pair (path length summed over crown parts and torus copies), linked to their pair by its index
struct geometryTargets {

	// Targets
	std::vector<int> targetsId;
	std::vector<double> eUnobsDirect; // Energy of unobstructed direct rays (in MJ/m2 on a horizontal plane)
	std::vector<double> eUnobsDiffuse; // Energy of unobstructed diffuse rays (in MJ/m2 on a horizontal plane)

	// Pairs
	std::vector<int> pairsTarget;
	std::vector<int> pairsRay;
	std::vector<int> pairsBlocked; // Is the ray intercepted by a trunk ?
	std::vector<int> pairsNTrees; // Number of crowns crossed (i.e. number of crossings of the pair)

	// Crossings
	std::vector<int> crossingsPair; // Index of the pair (row in the pairs)
	std::vector<int> crossingsTree;
	std::vector<double> crossingsLength;
};

// Rays used in the geometry
struct geometryRays {
	std::vector<double> eHorizontal; // Energy of the ray above canopy (in MJ/m2 on a horizontal plane)
	std::vector<double> eSlope; // Energy of the ray above canopy (in MJ/m2 on the slope plane)
	std::vector<int> direct; // Is the ray a direct one ?
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

	// Two points are equal (used to avoid duplicated points in crown fractions)
	bool areEquals(vertex3D p1, vertex3D p2) {
		return(std::abs(p1.x - p2.x) < EPSILON_FRACTION &&
			std::abs(p1.y - p2.y) < EPSILON_FRACTION &&
			std::abs(p1.z - p2.z) < EPSILON_FRACTION);
	}

	// Functions for finding if a point is in a volume
	// tolerance accounts for rounding errors (only used for crown fractions, as in Capsis; 0 otherwise)
	bool isInBBox(vertex3D pmin, vertex3D pmax, vertex3D p, double tolerance = 0.0) {

		return (pmin.x <= p.x + tolerance && p.x <= pmax.x + tolerance &&
			pmin.y <= p.y + tolerance && p.y <= pmax.y + tolerance &&
			pmin.z <= p.z + tolerance && p.z <= pmax.z + tolerance);
	}

	// Is the point on the periphery of the bounding box (with tolerance for rounding errors, as in Capsis)
	bool isOnBBoxPeriphery(vertex3D pmin, vertex3D pmax, vertex3D p, double tolerance) {
		return (std::abs(pmin.x - p.x) <= tolerance || std::abs(p.x - pmax.x) <= tolerance ||
			std::abs(pmin.y - p.y) <= tolerance || std::abs(p.y - pmax.y) <= tolerance ||
			std::abs(pmin.z - p.z) <= tolerance || std::abs(p.z - pmax.z) <= tolerance);
	}

	// Is an odd number of interception points with a crown part an expected particular case (as in Capsis) ?
	// Particular cases are zero-length interceptions. As any odd number of points, they are ignored,
	// but only the unexpected cases are counted and reported to the user
	//   1. the top of the crown part is at the height of the target point
	//   2. the only point is on the crown surface at the target point (z = 0), the other root is below the target
	//   3. (crown fractions) the only point is the target point, on the periphery of the bounding box
	//   4. (crown fractions) the ray crosses an axis of the crown fraction (two cut-plane points are equal)
	bool isParticularCase(int nsols, solutionType sol_type, double sol, Ray* ray, double z_top, bool is_fraction,
		vertex3D pmin, vertex3D pmax, vertex3D ptarget, vertex3D px0, vertex3D py0, vertex3D pz0, double tolerance) {

		// Case 1
		if (z_top == 0.0) { return(true); }

		// Cases 2 to 4 only concern a single interception point
		if (nsols != 1) { return(false); }

		// Case 2
		if (sol_type == SOL_SURFACE) { return(this->getInterceptionPoint(sol, ray).z == 0.0); }

		// Cases 3 and 4 only for crown fractions
		if (!is_fraction) { return(false); }

		// Case 3
		if (sol_type == SOL_TARGET && this->isOnBBoxPeriphery(pmin, pmax, ptarget, tolerance)) { return(true); }

		// Case 4
		return(this->areEquals(py0, px0) || this->areEquals(pz0, px0) || this->areEquals(pz0, py0));
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
	double a; // X-axis
	double b; // Y-axis
	double c; // Z-axis (height)


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

		// Target point on a vertical cut plane of a crown fraction (within rounding errors): set it exactly on the plane
		// (avoids ambiguous interception points within the tolerance on both sides of the plane)
		if (this->is8th) {
			if (std::abs(p0_shift.x) <= EPSILON_FRACTION) { p0_shift.x = 0.0; }
			if (std::abs(p0_shift.y) <= EPSILON_FRACTION) { p0_shift.y = 0.0; }
		}

		// Ray lying in a vertical cut plane of a crown fraction (x = x0 or y = y0):
		// the ray runs along the boundary between two fractions, which would both count the same path.
		// The path is only attributed to the fraction on the positive side of the plane, and the cut plane
		// is not considered as an interception point (departure from Capsis, which counts the path twice)
		bool ray_in_plane_x = this->is8th && p0_shift.x == 0.0 &&
			std::abs(ray->getCosHeightAngle() * ray->getCosAzimuth()) <= EPSILON_FRACTION;
		bool ray_in_plane_y = this->is8th && p0_shift.y == 0.0 &&
			std::abs(ray->getCosHeightAngle() * ray->getSinAzimuth()) <= EPSILON_FRACTION;

		if ((ray_in_plane_x && this->a < 0.0) || (ray_in_plane_y && this->b < 0.0)) { return(nullptr); }

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

		// Tolerance for rounding errors in the bounding box tests, only for crown fractions (as in Capsis)
		double tolerance = this->is8th ? EPSILON_FRACTION : 0.0;

		// Interception points (distance to the target) and their type (to identify particular cases)
		std::vector<double> sols;
		std::vector<solutionType> sols_type;
		auto add_sol = [&sols, &sols_type](double sol, solutionType type) {
			sols.push_back(sol);
			sols_type.push_back(type);
		};

		// The first interception point with the full ellipsoid is in the eight ellipsoid
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc1, tolerance) && pc1.z >= 0.0)
			add_sol(solc1, SOL_SURFACE);

		// The second interception point with the full ellipsoid is in the eight ellipsoid
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc2, tolerance) && pc2.z >= 0.0)
			add_sol(solc2, SOL_SURFACE);

		// The ray intersects the center of the target cell
		if (this->isInBBox(p_bbox_min, p_bbox_max, ptarget, tolerance) && this->isInEllipsoid(p0_shift, ptarget))
			add_sol(0.0, SOL_TARGET);

		// The ray intersects the Z horizontal plane of the eight ellipsoid
		if (this->isSemi) {
			if (this->isInBBox(p_bbox_min, p_bbox_max, pz0, tolerance) && this->isInEllipsoid(p0_shift, pz0) && pz0.z > 0.0 &&
				!(this->is8th && this->areEquals(pz0, ptarget)))
				add_sol(solz0, SOL_PLANE);
		}

		// Add interception with the verticla planes of the 8th ellipsoid
		if (this->is8th) {
			// The ray intersects the X vertical plane of the eight ellipsoid
			if (!ray_in_plane_x && this->isInBBox(p_bbox_min, p_bbox_max, px0, tolerance) && this->isInEllipsoid(p0_shift, px0) && px0.z > 0.0 && !this->areEquals(px0, pz0) &&
				!this->areEquals(px0, ptarget))
				add_sol(solx0, SOL_PLANE);

			// The ray intersects the Y vertical plane of the eight ellipsoid
			if (!ray_in_plane_y && this->isInBBox(p_bbox_min, p_bbox_max, py0, tolerance) && this->isInEllipsoid(p0_shift, py0) && py0.z > 0.0 && !this->areEquals(py0, pz0) && !this->areEquals(py0, px0) &&
				!this->areEquals(py0, ptarget))
				add_sol(soly0, SOL_PLANE);
		}


		// FIND PATH LENGTH AND DISTANCE TO CELL CENTER
		int nsols = sols.size();

		// No interception if not 2 solutions (as in Capsis)
		// 0: no interception, 1 or 3: particular case (zero-length interception) or unexpected case (counted and reported in R)
		if (nsols != 2) {
			if (nsols > 0 && !this->isParticularCase(nsols, sols_type[0], sols[0], ray, p_bbox_max.z, this->is8th,
				p_bbox_min, p_bbox_max, ptarget, px0, py0, pz0, tolerance)) {
				g_n_unexpected_ellipsoid++;
			}
			return(nullptr);
		}

		// Distance between target cell and middle point between two interceptions
		// Length is the length of the ray path across the crown
		double sol1 = sols.back();
		sols.pop_back();
		double sol2 = sols.back();

		double length = std::abs(sol1 - sol2);
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
	double a; // X-axis
	double b; // Y-axis
	double h; // Z-axis (height)


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

		// Target point on a vertical cut plane of a crown fraction (within rounding errors): set it exactly on the plane
		// (avoids ambiguous interception points within the tolerance on both sides of the plane)
		if (this->is4th) {
			if (std::abs(p0_shift.x) <= EPSILON_FRACTION) { p0_shift.x = 0.0; }
			if (std::abs(p0_shift.y) <= EPSILON_FRACTION) { p0_shift.y = 0.0; }
		}

		// Ray lying in a vertical cut plane of a crown fraction (x = x0 or y = y0):
		// the ray runs along the boundary between two fractions, which would both count the same path.
		// The path is only attributed to the fraction on the positive side of the plane, and the cut plane
		// is not considered as an interception point (departure from Capsis, which counts the path twice)
		bool ray_in_plane_x = this->is4th && p0_shift.x == 0.0 &&
			std::abs(ray->getCosHeightAngle() * ray->getCosAzimuth()) <= EPSILON_FRACTION;
		bool ray_in_plane_y = this->is4th && p0_shift.y == 0.0 &&
			std::abs(ray->getCosHeightAngle() * ray->getSinAzimuth()) <= EPSILON_FRACTION;

		if ((ray_in_plane_x && this->a < 0.0) || (ray_in_plane_y && this->b < 0.0)) { return(nullptr); }
		

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

		// Tolerance for rounding errors in the bounding box tests, only for crown fractions (as in Capsis)
		double tolerance = this->is4th ? EPSILON_FRACTION : 0.0;

		// Interception points (distance to the target) and their type (to identify particular cases)
		std::vector<double> sols;
		std::vector<solutionType> sols_type;
		auto add_sol = [&sols, &sols_type](double sol, solutionType type) {
			sols.push_back(sol);
			sols_type.push_back(type);
		};

		// The first interception point with the the infinite paraboloid is in the paraboloidal crown
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc1, tolerance) && pc1.z >= 0.0)
			add_sol(solc1, SOL_SURFACE);

		// The second interception point with the infinite paraboloid is in the paraboloidal crown
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc2, tolerance) && pc2.z >= 0.0)
			add_sol(solc2, SOL_SURFACE);

		// The ray intersects the Z horizontal base of the paraboloid
		if (this->isInBBox(p_bbox_min, p_bbox_max, pbase, tolerance) && this->isInParaboloid(p0_shift, pbase) && pbase.z > 0.0)
			add_sol(solbase, SOL_PLANE);

		// The ray intersects the center of the target cell
		// (for a 4th of paraboloid, the target must also be above the base plane, as in Capsis)
		if (this->isInBBox(p_bbox_min, p_bbox_max, ptarget, tolerance) && this->isInParaboloid(p0_shift, ptarget) &&
			(!this->is4th || ptarget.z >= pbase.z))
			add_sol(0.0, SOL_TARGET);

		// Add interception with the verticla planes of the 4th paraboloid
		if (this->is4th) {
			// The ray intersects the X vertical plane of the fourth paraboloid
			if (!ray_in_plane_x && this->isInBBox(p_bbox_min, p_bbox_max, px0, tolerance) && this->isInParaboloid(p0_shift, px0) && px0.z > 0.0 && !this->areEquals(px0, pbase) &&
				!this->areEquals(px0, ptarget))
				add_sol(solx0, SOL_PLANE);

			// The ray intersects the Y vertical plane of the fourth paraboloid
			if (!ray_in_plane_y && this->isInBBox(p_bbox_min, p_bbox_max, py0, tolerance) && this->isInParaboloid(p0_shift, py0) && py0.z > 0.0 && !this->areEquals(py0, pbase) && !this->areEquals(py0, px0) &&
				!this->areEquals(py0, ptarget))
				add_sol(soly0, SOL_PLANE);
		}


		// FIND PATH LENGTH AND DISTANCE TO CELL CENTER
		int nsols = sols.size();

		// No interception if not 2 solutions (as in Capsis)
		// 0: no interception, 1 or 3: particular case (zero-length interception) or unexpected case (counted and reported in R)
		if (nsols != 2) {
			if (nsols > 0 && !this->isParticularCase(nsols, sols_type[0], sols[0], ray, p_bbox_max.z, this->is4th,
				p_bbox_min, p_bbox_max, ptarget, px0, py0, pbase, tolerance)) {
				g_n_unexpected_paraboloid++;
			}
			return(nullptr);
		}

		// Distance between target cell and middle point between two interceptions
		// Length is the length of the ray path across the crown
		double sol1 = sols.back();
		sols.pop_back();
		double sol2 = sols.back();

		double length = std::abs(sol1 - sol2);
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
	double crownOpenness; // between 0 and 1, percentage of ligth energy of a ray that is not absorbed by the crown (for porous envelop method)
	double crownLAD; // Leaf area density in m2/m3 (for turbid medium method)
	
	// Output energies (in MJ)
	int nDirect; // Number of intercepted direct rays
	int nUnobstructedDirect; // Number of intercepted direct rays that was not intercepted by a previous other crown
	double energyDirect; // Energy intercepted by the tree from direct rays
	double energyPotentialDirect; // Potential energy intercepted by the tree (i.e. without neighbours) from direct rays
	double energyUnobstructedDirect; // Energy intercepted by the tree only from unobstructed direct rays (i.e. rays that have not been attenuated by other trees before)

	int nDiffuse; // Number of intercepted diffuse rays
	int nUnobstructedDiffuse; // Number of intercepted diffuse rays that was not intercepted by a previous other crown
	double energyDiffuse; // Energy intercepted by the tree from diffuse rays
	double energyPotentialDiffuse; // Potential energy intercepted by the tree (i.e. without neighbours) from diffuse rays
	double energyUnobstructedDiffuse; // Energy intercepted by the tree only from unobstructed diffuse rays (i.e. rays that have not been attenuated by other trees before)



private:

	// --- FUNCTIONS TO INITIALIZE DIFFERENT CROWN FORMS ---

	// Symetric paraboloid crown
	void initCrownParaboloid(double x, double y, double z,
		double h, double hbase,
		double cr_ymax, double cr_xmax, double cr_ymin, double cr_xmin) {

		// Radius is the same for all directions and is defined by mean crown radius of the four directions
		double crown_radius = (cr_ymax + cr_xmax + cr_ymin + cr_xmin) / 4.0;

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
		double cr_ymax, double cr_xmax, double cr_ymin, double cr_xmin) {

		// Compute height of the crown
		double cdepth = h - hbase;

		// Get center of the crown 
		// x, y, z are tree coordinates
		double x0 = x;
		double y0 = y;
		double z0 = z + hbase;

		// Init fourth parts of the crown
		this->crownParts.push_back(new CrownPartParaboloid(this->vectIdTree, true,
			x0, y0, z0, cr_xmax, cr_ymax, cdepth));
		this->crownParts.push_back(new CrownPartParaboloid(this->vectIdTree, true,
			x0, y0, z0, cr_xmax, -cr_ymin, cdepth));
		this->crownParts.push_back(new CrownPartParaboloid(this->vectIdTree, true,
			x0, y0, z0, -cr_xmin, -cr_ymin, cdepth));
		this->crownParts.push_back(new CrownPartParaboloid(this->vectIdTree, true,
			x0, y0, z0, -cr_xmin, cr_ymax, cdepth));
	}

	// Symetric ellispoidal crown
	void initCrownEllipsoid(double x, double y, double z,
		double h, double hbase,
		double cr_ymax, double cr_xmax, double cr_ymin, double cr_xmin) {

		// Radius is the same for all directions and is defined by mean crown radius of the four directions
		double crown_radius = (cr_ymax + cr_xmax + cr_ymin + cr_xmin) / 4.0;

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
		double cr_ymax, double cr_xmax, double cr_ymin, double cr_xmin) {

		// Radius is the same for all directions and is defined by mean crown radius of the four directions
		double crown_radius = (cr_ymax + cr_xmax + cr_ymin + cr_xmin) / 4.0;

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
		double cr_ymax, double cr_xmax, double cr_ymin, double cr_xmin) {

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
				x0, y0, z0, cr_xmax, cr_ymax, depth_up));
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, cr_xmax, -cr_ymin, depth_up));
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, -cr_xmin, -cr_ymin, depth_up));
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, -cr_xmin, cr_ymax, depth_up));
		}

		// Init bottom eighth parts of the crown
		if (depth_down > 0) {
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, cr_xmax, cr_ymax, -depth_down));
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, cr_xmax, -cr_ymin, -depth_down));
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, -cr_xmin, -cr_ymin, -depth_down));
			this->crownParts.push_back(new CrownPartEllipsoid(this->vectIdTree, true, true,
				x0, y0, z0, -cr_xmin, cr_ymax, -depth_down));
		}
	}


public:
	Crown(int vectid_tree, 
		std::string crown_type,
		double x, double y, double z,
		double h, double hbase, double hmax,
		double cr_ymax, double cr_xmax, double cr_ymin, double cr_xmin,
		double crown_openness, double crown_lad)
	{
		this->vectIdTree = vectid_tree;

		this->crownOpenness = crown_openness;
		this->crownLAD = crown_lad;

		this->nDirect = 0;
		this->nUnobstructedDirect = 0;
		this->energyDirect = 0.0;
		this->energyPotentialDirect = 0.0;
		this->energyUnobstructedDirect = 0.0;

		this->nDiffuse = 0;
		this->nUnobstructedDiffuse = 0;
		this->energyDiffuse = 0.0;
		this->energyPotentialDiffuse = 0.0;
		this->energyUnobstructedDiffuse = 0.0;


		// Init the crown parts depending on the type
		if (crown_type == "8E")
			this->initCrownEllipsoid8th(x, y, z,
				h, hbase, hmax, 
				cr_ymax, cr_xmax, cr_ymin, cr_xmin);
		else if (crown_type == "2E") 
			this->initCrownEllipsoidSemi(x, y, z,
				h, hbase, hmax, 
				cr_ymax, cr_xmax, cr_ymin, cr_xmin);
		else if (crown_type == "E")
			this->initCrownEllipsoid(x, y, z,
				h, hbase, 
				cr_ymax, cr_xmax, cr_ymin, cr_xmin);
		else if (crown_type == "4P")
			this->initCrownParaboloid4th(x, y, z,
				h, hbase,
				cr_ymax, cr_xmax, cr_ymin, cr_xmin);
		else if (crown_type == "P") 
			this->initCrownParaboloid(x, y, z,
				h, hbase, 
				cr_ymax, cr_xmax, cr_ymin, cr_xmin);
		else
			Rcpp::Rcout << "Unrecognized crown type" << std::endl;
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
	double getCrownOpenness() { return(this->crownOpenness); }
	double getCrownLAD() { return(this->crownLAD); }

	double getRci() { return( 1.0 - static_cast<double>(this->nUnobstructedDirect + this->nUnobstructedDiffuse) / static_cast<double>(this->nDirect + this->nDiffuse) ); }
	double getEnergy() { return(this->energyDiffuse + this->energyDirect); }
	double getEnergyPotential() { return(this->energyPotentialDiffuse + this->energyPotentialDirect); }
	double getEnergyUnobstructed() { return(this->energyUnobstructedDiffuse + this->energyUnobstructedDirect); }

	double getRciDirect() { return(1.0 - static_cast<double>(this->nUnobstructedDirect) / static_cast<double>(this->nDirect)); }
	double getEnergyDirect() { return(this->energyDirect); }
	double getEnergyPotentialDirect() { return(this->energyPotentialDirect); }
	double getEnergyUnobstructedDirect() { return(this->energyUnobstructedDirect); }

	double getRciDiffuse() { return(1.0 - static_cast<double>(this->nUnobstructedDiffuse) / static_cast<double>(this->nDiffuse)); }
	double getEnergyDiffuse() { return(this->energyDiffuse); }
	double getEnergyPotentialDiffuse() { return(this->energyPotentialDiffuse); }
	double getEnergyUnobstructedDiffuse() { return(this->energyUnobstructedDiffuse); }



	// Setter
	void addInterceptedRayDirect() { this->nDirect += 1; }
	void addInterceptedUnobstructedRayDirect() { this->nUnobstructedDirect += 1; }
	void addEnergyDirect(double e) { this->energyDirect += e; }
	void addEnergyPotentialDirect(double e) { this->energyPotentialDirect += e; }
	void addEnergyUnobstructedDirect(double e) { this->energyUnobstructedDirect += e; }

	void addInterceptedRayDiffuse() { this->nDiffuse += 1; }
	void addInterceptedUnobstructedRayDiffuse() { this->nUnobstructedDiffuse += 1; }
	void addEnergyDiffuse(double e) { this->energyDiffuse += e; }
	void addEnergyPotentialDiffuse(double e) { this->energyPotentialDiffuse += e; }
	void addEnergyUnobstructedDiffuse(double e) { this->energyUnobstructedDiffuse += e; }


	// Find interception with a given part of the crown
	interception* computeCrownpartInterception(int id_crownpart, Ray* ray, vertex3D& shift) {
		return( this->crownParts[id_crownpart]->computeInterception(ray, shift) );
	}
};

class Trunk : public TreeVolume {

private:

	// Position of the trunk base is stored in TreeVolume (x, y, z)

	// Dimensions of the trunk (in meters)
	// The trunk is a cylinder from the ground to the crown base (as in Capsis)
	double crownBase; // Height of the crown base above the trunk base
	double radius;

private:

	bool isInCylinder(vertex3D p0, vertex3D p) {
		return (p0.x - p.x) * (p0.x - p.x) + (p0.y - p.y) * (p0.y - p.y) <= this->radius * this->radius;
	}

		

public:
	Trunk(int vectid_tree, double x, double y, double z,
		double dbh_cm, double crown_base_m) :
		TreeVolume(vectid_tree, x, y, z)
	{
		this->crownBase = crown_base_m;
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
		if (std::abs(p0_shift.z) < EPSILON_TRUNK_BASE)
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

		// Interception with the top horizontal base of the trunk (i.e. the crown base)
		double ztop = p0_shift.z + this->crownBase;
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

		// Interception points (distance to the target) and their type (to identify particular cases)
		std::vector<double> sols;
		std::vector<solutionType> sols_type;
		auto add_sol = [&sols, &sols_type](double sol, solutionType type) {
			sols.push_back(sol);
			sols_type.push_back(type);
		};

		// The first interception point with the infinite cylinder is in the trunk
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc1, EPSILON_TRUNK) && pc1.z >= 0.0 && pc1.z >= p0_shift.z) 
			add_sol(solc1, SOL_SURFACE);

		// The second interception point with the infinite cylinder is in the trunk
		if (this->isInBBox(p_bbox_min, p_bbox_max, pc2, EPSILON_TRUNK) && pc2.z >= 0.0 && pc2.z >= p0_shift.z)
			add_sol(solc2, SOL_SURFACE);

		// The ray intersects the ground or top horizontal plane of the trunk
		if (this->isInCylinder(p0_shift, pztop) && ztop >= 0.0 && ztop >= p0_shift.z)
			add_sol(solztop, SOL_PLANE);
		if (this->isInCylinder(p0_shift, pzbot) && p0_shift.z >= 0.0)
			add_sol(solzbot, SOL_PLANE);

		// The ray intersects the center of the target cell
		if (this->isInBBox(p_bbox_min, p_bbox_max, ptarget, EPSILON_TRUNK) && this->isInCylinder(p0_shift, ptarget) && ptarget.z != pzbot.z && ptarget.z != pztop.z)
			add_sol(0.0, SOL_TARGET);


		// FIND PATH LENGTH AND DISTANCE TO CELL CENTER
		int nsols = sols.size();

		// No interception if not 2 solutions (as in Capsis)
		// 0: no interception, 1 or 3: particular case (zero-length interception) or unexpected case (counted and reported in R)
		//   Particular cases (as in Capsis SLTrunk):
		//   1. the crown base (top of the trunk) is at the height of the target point
		//   2. the only point is on the trunk surface, at the height of the target point or of the crown base
		if (nsols != 2) {
			bool particular_case = (ztop == 0.0);
			if (nsols == 1 && sols_type[0] == SOL_SURFACE) {
				double z_sol = this->getInterceptionPoint(sols[0], ray).z;
				if (z_sol == 0.0 || z_sol == ztop) { particular_case = true; }
			}
			if (nsols > 0 && !particular_case) {
				g_n_unexpected_trunk++;
			}
			return(nullptr);
		}

		// Distance between target cell and middle point between two interceptions
		// Length is the length of the ray path across the crown
		double sol1 = sols.back();
		sols.pop_back();
		double sol2 = sols.back();

		double length = std::abs(sol1 - sol2);
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
		double x_crown, double y_crown,
		double hbase, double hmax,
		double cr_ymax, double cr_xmax, double cr_ymin, double cr_xmin,
		double crown_openness, double crown_lad):

		crown(vectid, crown_type, x_crown, y_crown, z, h, hbase, hmax, cr_ymax, cr_xmax, cr_ymin, cr_xmin,
			crown_openness, crown_lad),
		trunk(vectid, x, y, z, dbh, hbase)
	{
		this->vectId = vectid;
		this->id = id;
	}


	// Getters
	int getId() { return(this->id); }
	int getVectId() { return(this->vectId); }

	Trunk& getTrunk() { return(this->trunk); }
	Crown& getCrown() { return(this->crown); }

	double getCrownRci() { return(this->crown.getRci()); }
	double getCrownEnergy() { return(this->crown.getEnergy()); }
	double getCrownEnergyPotential() { return(this->crown.getEnergyPotential()); }
	double getCrownEnergyUnobstructed() { return(this->crown.getEnergyUnobstructed()); }

	double getCrownRciDirect() { return(this->crown.getRciDirect()); }
	double getCrownEnergyDirect() { return(this->crown.getEnergyDirect()); }
	double getCrownEnergyPotentialDirect() { return(this->crown.getEnergyPotentialDirect()); }
	double getCrownEnergyUnobstructedDirect() { return(this->crown.getEnergyUnobstructedDirect()); }

	double getCrownRciDiffuse() { return(this->crown.getRciDiffuse()); }
	double getCrownEnergyDiffuse() { return(this->crown.getEnergyDiffuse()); }
	double getCrownEnergyPotentialDiffuse() { return(this->crown.getEnergyPotentialDiffuse()); }
	double getCrownEnergyUnobstructedDiffuse() { return(this->crown.getEnergyUnobstructedDiffuse()); }


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

	// Energy incident on the target (in MJ)
	// will be decreased by successive interceptions by above trees
	double energyDirectSlope; // Current direct energy on the target on the slope
	double energyDiffuseSlope; // Current diffuse energy on the target on the slope
	double energyDirectHorizontal; // Current direct energy on the target on a horizontal plane
	double energyDiffuseHorizontal; // Current diffuse energy on the target on a horizontal plane

	// Energy incident on the target coming from unobstructed rays (in MJ)
	// will be increased by the total energy of rays that thave been not been intercepted at least once by a tree
	double energyUnobstructedDirectSlope; // Unobstructed direct energy on the target on the slope
	double energyUnobstructedDiffuseSlope; // Unobstructed diffuse energy on the target on the slope
	double energyUnobstructedDirectHorizontal; // Unobstructed direct energy on the target on a horizontal plane
	double energyUnobstructedDiffuseHorizontal; // Unobstructed diffuse energy on the target on a horizontal plane

	// Vectors of interceptions between the tree and the target
	// Here, also store the total energies, whereas above, we compute it only in exportResults
	// Indeed, it is more efficient when manipulation large matrices to increment the total energy during the light interception computation
	// Rather to sum the two matrices of direct and diffuse energies at the end
	std::vector<int> interceptionsTreesNumber;
	std::vector<int> interceptionsTreesNumberDirect;
	std::vector<int> interceptionsTreesNumberDiffuse;

	std::vector<int> interceptionsTreesNumberUnobstructed;
	std::vector<int> interceptionsTreesNumberUnobstructedDirect;
	std::vector<int> interceptionsTreesNumberUnobstructedDiffuse;

	std::vector<double> interceptionsTreesEnergySlope;
	std::vector<double> interceptionsTreesEnergyDirectSlope;
	std::vector<double> interceptionsTreesEnergyDiffuseSlope;

	std::vector<double> interceptionsTreesEnergyHorizontal;
	std::vector<double> interceptionsTreesEnergyDirectHorizontal;
	std::vector<double> interceptionsTreesEnergyDiffuseHorizontal;


	// If the target is a sensor (other case is a cell)
	// In this case, do not remove energies from intercepted crown
	bool isSensor;


public:
	Target(double x, double y, double z, 
		int row, int col, 
		double e_above_direct_slope, double e_above_diffuse_slope,
		double e_above_direct_horizontal, double e_above_diffuse_horizontal,
		int n_total_trees,
		bool is_sensor) :

		// Initialize the vectors of interceptions with each tree (using the vector constructor)
		interceptionsTreesNumber(n_total_trees, 0),
		interceptionsTreesNumberDirect(n_total_trees, 0),
		interceptionsTreesNumberDiffuse(n_total_trees, 0),

		interceptionsTreesNumberUnobstructed(n_total_trees, 0),
		interceptionsTreesNumberUnobstructedDirect(n_total_trees, 0),
		interceptionsTreesNumberUnobstructedDiffuse(n_total_trees, 0),

		interceptionsTreesEnergySlope(n_total_trees, 0.0),
		interceptionsTreesEnergyDirectSlope(n_total_trees, 0.0),
		interceptionsTreesEnergyDiffuseSlope(n_total_trees, 0.0),

		interceptionsTreesEnergyHorizontal(n_total_trees, 0.0),
		interceptionsTreesEnergyDirectHorizontal(n_total_trees, 0.0),
		interceptionsTreesEnergyDiffuseHorizontal(n_total_trees, 0.0)
	{
		// Geometry of the target
		this->x = x;
		this->y = y;
		this->z = z;
		this->row = row;
		this->col = col;

		// The energy will be decreased by successive interception by the above trees
		this->energyDirectSlope = e_above_direct_slope;
		this->energyDiffuseSlope = e_above_diffuse_slope;
		this->energyDirectHorizontal = e_above_direct_horizontal;
		this->energyDiffuseHorizontal = e_above_diffuse_horizontal;

		// The energy will be increased by the total energy of rays that thave been not been intercepted at least once by a tree
		this->energyUnobstructedDirectSlope = 0.0;
		this->energyUnobstructedDiffuseSlope = 0.0;
		this->energyUnobstructedDirectHorizontal = 0.0;
		this->energyUnobstructedDiffuseHorizontal = 0.0;

		// Is it a sensor ?
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

	double getEnergyUnobstructedSlope() { return(this->energyUnobstructedDirectSlope + this->energyUnobstructedDiffuseSlope); }
	double getEnergyUnobstructedDirectSlope() { return(this->energyUnobstructedDirectSlope); }
	double getEnergyUnobstructedDiffuseSlope() { return(this->energyUnobstructedDiffuseSlope); }

	double getEnergyUnobstructedHorizontal() { return(this->energyUnobstructedDirectHorizontal + this->energyUnobstructedDiffuseHorizontal); }
	double getEnergyUnobstructedDirectHorizontal() { return(this->energyUnobstructedDirectHorizontal); }
	double getEnergyUnobstructedDiffuseHorizontal() { return(this->energyUnobstructedDiffuseHorizontal); }

	std::vector<int>& getInterceptionsTreesNumber() { return(this->interceptionsTreesNumber); }
	std::vector<int>& getInterceptionsTreesNumberDirect() { return(this->interceptionsTreesNumberDirect); }
	std::vector<int>& getInterceptionsTreesNumberDiffuse() { return(this->interceptionsTreesNumberDiffuse); }

	std::vector<int>& getInterceptionsTreesNumberUnobstructed() { return(this->interceptionsTreesNumberUnobstructed); }
	std::vector<int>& getInterceptionsTreesNumberUnobstructedDirect() { return(this->interceptionsTreesNumberUnobstructedDirect); }
	std::vector<int>& getInterceptionsTreesNumberUnobstructedDiffuse() { return(this->interceptionsTreesNumberUnobstructedDiffuse); }

	std::vector<double>& getInterceptionsTreesEnergySlope() { return(this->interceptionsTreesEnergySlope); }
	std::vector<double>& getInterceptionsTreesEnergyDirectSlope() { return(this->interceptionsTreesEnergyDirectSlope); }
	std::vector<double>& getInterceptionsTreesEnergyDiffuseSlope() { return(this->interceptionsTreesEnergyDiffuseSlope); }

	std::vector<double>& getInterceptionsTreesEnergyHorizontal() { return(this->interceptionsTreesEnergyHorizontal); }
	std::vector<double>& getInterceptionsTreesEnergyDirectHorizontal() { return(this->interceptionsTreesEnergyDirectHorizontal); }
	std::vector<double>& getInterceptionsTreesEnergyDiffuseHorizontal() { return(this->interceptionsTreesEnergyDiffuseHorizontal); }


	bool isThisSensor() { return(this->isSensor); }


	// Setters
	void interceptRayDirect(bool is_unobstructed, double e_slope, double e_horizontal, int vect_id_tree) {
		this->interceptionsTreesNumber[vect_id_tree] += 1;
		this->interceptionsTreesNumberDirect[vect_id_tree] += 1;

		if (is_unobstructed) {
			this->interceptionsTreesNumberUnobstructed[vect_id_tree] += 1;
			this->interceptionsTreesNumberUnobstructedDirect[vect_id_tree] += 1;
		}

		this->interceptEnergySlope(e_slope, vect_id_tree);
		this->interceptEnergyDirectSlope(e_slope, vect_id_tree);

		this->interceptEnergyHorizontal(e_horizontal, vect_id_tree);
		this->interceptEnergyDirectHorizontal(e_horizontal, vect_id_tree);
	}

	void interceptRayDiffuse(bool is_unobstructed, double e_slope, double e_horizontal, int vect_id_tree) {
		this->interceptionsTreesNumber[vect_id_tree] += 1;
		this->interceptionsTreesNumberDiffuse[vect_id_tree] += 1;

		if (is_unobstructed) {
			this->interceptionsTreesNumberUnobstructed[vect_id_tree] += 1;
			this->interceptionsTreesNumberUnobstructedDiffuse[vect_id_tree] += 1;
		}

		this->interceptEnergySlope(e_slope, vect_id_tree);
		this->interceptEnergyDiffuseSlope(e_slope, vect_id_tree);

		this->interceptEnergyHorizontal(e_horizontal, vect_id_tree);
		this->interceptEnergyDiffuseHorizontal(e_horizontal, vect_id_tree);
	}

	void addUnobstructedEnergiesDirect(double e_above_slope, double e_above_horizontal) {
		this->energyUnobstructedDirectSlope += e_above_slope;
		this->energyUnobstructedDirectHorizontal += e_above_horizontal;
	}

	void addUnobstructedEnergiesDiffuse(double e_above_slope, double e_above_horizontal) {
		this->energyUnobstructedDiffuseSlope += e_above_slope;
		this->energyUnobstructedDiffuseHorizontal += e_above_horizontal;
	}


	//void resetEnergy() { this->energy = 0.0; }

	//virtual void correctNullEnergy(double epsilon) {}

	private:
		void interceptEnergySlope(double e, int vect_id_tree) {
			this->interceptionsTreesEnergySlope[vect_id_tree] += e;
		}
		void interceptEnergyHorizontal(double e, int vect_id_tree) {
			this->interceptionsTreesEnergyHorizontal[vect_id_tree] += e;
		}

		void interceptEnergyDirectSlope(double e, int vect_id_tree) {
			this->energyDirectSlope -= e;
			this->interceptionsTreesEnergyDirectSlope[vect_id_tree] += e;
		}
		void interceptEnergyDirectHorizontal(double e, int vect_id_tree) {
			this->energyDirectHorizontal -= e;
			this->interceptionsTreesEnergyDirectHorizontal[vect_id_tree] += e;
		}

		void interceptEnergyDiffuseSlope(double e, int vect_id_tree) {
			this->energyDiffuseSlope -= e;
			this->interceptionsTreesEnergyDiffuseSlope[vect_id_tree] += e;
		}
		void interceptEnergyDiffuseHorizontal(double e, int vect_id_tree) {
			this->energyDiffuseHorizontal -= e;
			this->interceptionsTreesEnergyDiffuseHorizontal[vect_id_tree] += e;
		}
};

class Sensor : public Target {

private:
	// Id of the sensor
	int idSensor;


public:
	Sensor(double x, double y, double z, 
		int row, int col,
		double e_above_direct_slope, double e_above_diffuse_slope,
		double e_above_direct_horizontal, double e_above_diffuse_horizontal,
		int n_total_trees,
		int id_sensor) :
		Target(x, y, z, row, col,
			e_above_direct_slope, e_above_diffuse_slope,
			e_above_direct_horizontal, e_above_diffuse_horizontal,
			n_total_trees,
			true)
	{
		this->idSensor = id_sensor;
	}

	// Getters
	int getIdSensor() { return(this->idSensor); }

	// Methods
	/*void correctNullEnergy(double epsilon) {
		if (this->getEnergy() < -epsilon)
			Rcpp::Rcout << "Problem with energy in target " << this->getIdSensor() << ": energy of " << this->getEnergy() << "MJ";
		if (std::abs(this->getEnergy()) < epsilon)
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
		double e_above_direct_slope, double e_above_diffuse_slope,
		double e_above_direct_horizontal, double e_above_diffuse_horizontal,
		int n_total_trees,
		int id_cell) :
		Target(x, y, z, row, col, 
			e_above_direct_slope, e_above_diffuse_slope,
			e_above_direct_horizontal, e_above_diffuse_horizontal,
			n_total_trees,
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
	//		Rcpp::Rcout << "Problem with energy in cell " << this->getIdCell() << ": energy of " << this->getEnergy() << "MJ";
	//	if (std::abs(this->getEnergy()) < epsilon)
	//		this->resetEnergy();
	//}
};

class Stand {

private:
	// Vector of sensors
	int nSensors;
	bool areSensors;
	std::vector<Sensor*> sensors;

	// 2D grid of cells stored by <row, column>
	std::vector<std::vector<Cell*>> grid;

	// Vector of trees
	int nTrees;
	std::vector<Tree*> trees;
	CharacterVector rcppTreeIdChar; // used for naming the matrices of interceptions

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
	Stand(DataFrame trees, DataFrame sensors, DataFrame cells,
		double e_above_direct_slope_m2, double e_above_diffuse_slope_m2,
		double e_above_direct_horizontal_m2, double e_above_diffuse_horizontal_m2,
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


		// Number of trees
		this->nTrees = trees.nrows();

		// Number of sensors
		this->nSensors = sensors.nrows();
		this->areSensors = this->nSensors > 0;


		// Get energy above canopy that is coming toward a cell
		double e_above_direct_slope_cell = e_above_direct_slope_m2 * this->cellAreaSlope;
		double e_above_diffuse_slope_cell = e_above_diffuse_slope_m2 * this->cellAreaSlope;
		double e_above_direct_horizontal_cell = e_above_direct_horizontal_m2 * this->cellAreaHorizontal;
		double e_above_diffuse_horizontal_cell = e_above_diffuse_horizontal_m2 * this->cellAreaHorizontal;


		// Create vector of sensors from the sensors R DataFrame, if sensors are specified
		if (this->areSensors) {

			IntegerVector sensors_id = sensors["id_sensor"];
			NumericVector sensors_x = sensors["x"];
			NumericVector sensors_y = sensors["y"];
			NumericVector sensors_z = sensors["z"];

			for (int i = 0; i < this->nSensors; i++) {

				// Sensors are within the stand limits (checked in R by validate_sl_stand())
				// Find the row and column the sensor belong to
				int sensor_col = findColFromPosX(sensors_x[i]);
				int sensor_row = findRowFromPosY(sensors_y[i]);

				// Create and push the pointor to Sensor new instance
				// The sensor area is the same as the cell size
				// Area is not a problem when trating with PACL, but may be considered when analysing absolute values of energy
				this->sensors.push_back(new Sensor(
					sensors_x[i], sensors_y[i], // Position of the sensor
					sensors_z[i], // Z pos of the sensor is the height above the slope
					sensor_row, sensor_col, // Cell of the sensor
					e_above_direct_slope_cell, e_above_diffuse_slope_cell, // Direct and diffuse energy coming from above the canopy on the slope
					e_above_direct_horizontal_cell, e_above_diffuse_horizontal_cell, // Same but on a horizontal plane
					this->nTrees, // Number of trees to initialise the vector of interceptions with trees
					sensors_id[i] // Id of the sensor
				));
			}

		}


		// Create grid of cells and check consistency with cells dataframe from R
		IntegerVector cells_id = cells["id_cell"];
		NumericVector cells_x = cells["x_center"];
		NumericVector cells_y = cells["y_center"];
		NumericVector cells_z = cells["z_center"];
		IntegerVector cells_row = cells["row"];
		IntegerVector cells_col = cells["col"];
		
		int nR = cells_id.size();
		int expected = this->nCellsX * this->nCellsY;
		if (nR != expected)
			stop("R cells has %d rows, expected %d", nR, expected);

		for (int r = 0; r < this->nCellsY; r++) {
			std::vector<Cell*> cells_in_row;
			for (int c = 0; c < this->nCellsX; c++) {
				
				// Unique cell id
				int cell_id = this->nCellsX * r + c + 1;
				int idx = cell_id - 1;   // R row

				// Find position of cell center
				double x = this->cellSize * (c + 0.5f);
				double y = this->cellSize * (this->nCellsY - r - 0.5f); // Be careful, grid is first rows on the top but Y coordinate is greater for first rows
				double z = this->computeZ(x, y);

				// Consistency checks between R and CPP
				if (cells_id[idx] != cell_id)
					stop("Cell R row %d, id mismatch", idx);

				if (cells_row[idx] != r || cells_col[idx] != c)
					stop("Cell R row %d, row/col mismatch", idx);

				if (std::abs(cells_x[idx] - x) > 1e-12)
					stop("Cell R row %d, x mismatch", idx);

				if (std::abs(cells_y[idx] - y) > 1e-12)
					stop("Cell R row %d, y mismatch", idx);

				if (std::abs(cells_z[idx] - z) > 1e-10)
					stop("Cell R row %d, z mismatch (%.10f vs %.10f)",
						idx, cells_z[idx], z);


				// Create cell object
				cells_in_row.push_back(new Cell(
					x, y, z, // Position of the cell center
					r, c, // Row and column of the cell within the stand grid system
					e_above_direct_slope_cell, e_above_diffuse_slope_cell, // Direct and diffuse energy coming from above the canopy on the slope
					e_above_direct_horizontal_cell, e_above_diffuse_horizontal_cell, // Same but on a horizontal plane
					this->nTrees, // Number of trees to initialise the vector of interceptions with trees
					cell_id // Id of the cell
				));
			}
			this->grid.push_back(cells_in_row);
		}


		// Create vector of trees from the trees R dataframe
		IntegerVector trees_id = trees["id_tree"];
		NumericVector trees_x = trees["x"];
		NumericVector trees_y = trees["y"];
		NumericVector trees_z = trees["z"];
		NumericVector dbh = trees["dbh_cm"];
		NumericVector tree_height = trees["h_m"];
		NumericVector hbase = trees["hbase_m"];
		NumericVector hmax = trees["hmax_m"];
		NumericVector cr_ymax = trees["rymax_m"];
		NumericVector cr_xmax = trees["rxmax_m"];
		NumericVector cr_ymin = trees["rymin_m"];
		NumericVector cr_xmin = trees["rxmin_m"];
		StringVector ctype = trees["crown_type"];
		NumericVector cp = trees["crown_openness"];
		NumericVector clad = trees["crown_lad"];


		// Observe if column x_crown and y_crown exists, otherwise set them to x and y
		// i.e. decentered trunk
		NumericVector trees_xcrown = clone(trees_x); // Initialize crown coords as x by default
		NumericVector trees_ycrown = clone(trees_y); // Initialize crown coords as y by default

		CharacterVector cnames = trees.names(); // Get column names
		if (std::find(cnames.begin(), cnames.end(), "x_crown") != cnames.end()) { // Check if "x_crown" and "y_crown" exist
			trees_xcrown = trees["x_crown"];
		}

		if (std::find(cnames.begin(), cnames.end(), "y_crown") != cnames.end()) {
			trees_ycrown = trees["y_crown"];
		}


		// Add a Tree object to the associated Cell they belong to
		for (int i = 0; i < this->nTrees; i++) {

			// Create tree object
			this->trees.push_back(new Tree(
				i, trees_id[i], trees_x[i], trees_y[i], trees_z[i],
				dbh[i], tree_height[i], 
				Rcpp::as< std::string >(ctype[i]), 
				trees_xcrown[i], trees_ycrown[i],
				hbase[i], hmax[i],
				cr_ymax[i], cr_xmax[i], cr_ymin[i], cr_xmin[i], cp[i], clad[i]));

			// Add to trees id vectors (i.e. corrrespondance vector between trees id in vect and trees id)
			this->rcppTreeIdChar.push_back(std::to_string(trees_id[i]));

			// Find row and column of the cell and add tree id to the corresponding cell
			int tree_col = findColFromPosX(trees_x[i]);
			int tree_row = findRowFromPosY(trees_y[i]);
			this->grid[tree_row][tree_col]->addTree(i);

		}

		// Get maximum height and crown radius of the trees
		this->maxTreeHeight = max(tree_height);
		this->maxTreeCrownRadius = std::max({ max(cr_ymax), max(cr_xmax), max(cr_ymin), max(cr_xmin) });
	}

	~Stand() {
		// Delete sensors pointors
		for (int i = 0; i < this->nSensors; i++) {
			delete this->sensors[i];
		}

		// Delete cells pointers
		for (int r = 0; r < this->nCellsY; r++) {
			for (int c = 0; c < this->nCellsX; c++) {
				delete this->grid[r][c];
			}
		}
		// Delete tree pointers
		for (int i = 0; i < this->nTrees; i++) {
			delete this->trees[i];
		}
	}


	// Getters
	bool areThereSensors() { return(this->areSensors); }
	Sensor* getSensor(int vectid) { return(this->sensors[vectid]); }
	int getNSensors() { return(this->nSensors); }

	Cell* getCell(int row, int col) { return(this->grid[row][col]); }
	int getNCells() { return(this->nCells); }
	int getNCellsX() { return(this->nCellsX); }
	int getNCellsY() { return(this->nCellsY); }

	Tree* getTree(int vectid) { return(this->trees[vectid]); }
	int getNTrees() { return(this->nTrees); }
	CharacterVector& getRcppTreeIdChar() { return(this->rcppTreeIdChar); }

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
	bool almost_equal(double a, double b, double tol = 1e-6) {
		return std::abs(a - b) <= tol;
	}
	
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

	// Column and row of the cell containing a point
	// A point exactly on the upper x (or y) limit of the stand belongs to the last column (or first row)
	int findColFromPosX(double x) { return(std::min((int)(x / this->cellSize), this->nCellsX - 1)); }
	int findRowFromPosY(double y) { return(std::max(this->nCellsY - (int)(y / this->cellSize) - 1, 0)); }


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

	// Global variable for turbid medium
	double extinctionCoef; // Probability of a leaf to intercept the ray (linked to leaf orientation)
	double clumpingFactor; // Aggregation of leaves within the crown volume (1 is homogeneous)

	// Consider interception with trunks ?
	bool trunkInterception;

	// Compute interceptions only of sensors ?
	// i.e. do not cast rays towards cells and do not predict tree energy interception
	bool sensorsOnly;

	// Geometry of the rays toward the sensors and cells, independent from the crown LAD and openness (see computeGeometry())
	geometryRays raysGeometry;
	geometryTargets sensorsGeometry;
	geometryTargets cellsGeometry;



private:

	shiftedCell getPotentialCell(Ray* ray, int rel_cell_id, Target* target) {

		relativeCoords rel_cell = ray->getPotentialRelCell(rel_cell_id);

		// Get row and column of potential cell
		int row_pot = target->getRow() + rel_cell.row;
		int col_pot = target->getCol() + rel_cell.col;

		// Search if cells is outside the this->stand
		bool is_outside = row_pot < 0 || row_pot >= this->stand.getNCellsY() || col_pot < 0 || col_pot >= this->stand.getNCellsX();

		// If cel is outside the main this->stand, do not keep as potential if we are not in a torus system
		if (is_outside && !this->stand.isUseTorus()) { return shiftedCell{ nullptr, {0.0, 0.0, 0.0} }; }

		// Get row and column of the original cell in the grid cell (diffrent only within a torus system, if cell is outside)
		int row_pot_original = (this->stand.getNCellsY() + (row_pot % this->stand.getNCellsY())) % this->stand.getNCellsY();
		int col_pot_original = (this->stand.getNCellsX() + (col_pot % this->stand.getNCellsX())) % this->stand.getNCellsX();

		// Get corresponding original cell
		Cell* original_cell = this->stand.getCell(row_pot_original, col_pot_original);

		// Check if it contains trees, otherwise, do not add as potential cell
		if ((*original_cell).isEmpty()) { return shiftedCell{ nullptr, { 0.0, 0.0, 0.0 } }; }

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

		// Find interceptions
		std::vector<interception*> v_interc = this->findInterceptionsRayToTarget(ray, target);

		// Order interceptions
		this->orderInterceptionsRayToTarget(v_interc);

		// Summarize interceptions into energy on the target and in the tree crown
		this->summarizeInterceptionsRayToTarget(ray, target, v_interc);

	}

	std::vector<interception*> findInterceptionsRayToTarget(Ray* ray, Target* target) {

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
		std::stable_sort(v_interc.begin(), v_interc.end(), [](interception const* i1, interception const* i2) {
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

		// Track the first intercepted tree has been intercepted, to compute unobstructed energy by the tree
		// i.e. add energy unobstructed to the crown only if it is the first interception of the ray with a crown
		int idtree_first_intercepted = -1;

		// Necessary when many crown parts, i.e. 4P/4E/8E
		// Check if a crown has already been intercepted
		// in porous envelop: as it does not consider path across the crown, attenuate the ray only at the first crown part intercepted)
		std::vector<bool> is_intercepted(this->stand.getNTrees(), false);

		// Necessary when many crown parts, i.e. 4P/4E/8E
		// Track the incident potential energy that arrive to a crown part (for potential energy in turbid medium)
		// Indeed, potential energy of the crown consider attenuation of energy between crown parts
		std::vector<double> e_pot_incident(this->stand.getNTrees(), e_incident_slope_cell);


		//Compute attenuation of energy across successives crown interceptions for each ray X target cell
		int n_interceptions = v_interc.size();
		for (int j = 0; j < n_interceptions; j++) {

			// Is the interception with a trunk (otherwise with a crown part) ?
			bool is_trunk = v_interc[j]->withTrunk;

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

			// Trunk ==> intercepts all the remaining energy of the ray (as in Capsis)
			// The following crown parts still intercept a null energy, but their potential energy is computed
			// The trunk does not modify the potential energy of the crown
			if (is_trunk) {
				intercepted_energy_slope = current_energy_slope;
				intercepted_energy_horizontal = current_energy_horizontal;
			}
			// Turbid medium ==> apply beer lambert law
			else if (this->turbidMedium) {

				// Apply Beer Lambert law on incdent energy attenuated by previous crown parts of the tree, without considering crowns of other competitors
				potential_energy_slope = this->applyBeerLambert(
					e_pot_incident[v_interc[j]->vectIdTree],
					this->extinctionCoef, this->clumpingFactor,
					crown.getCrownLAD(),
					v_interc[j]->length);

				// When the beam penetrates a tree not yet traversed,
				// we reset the beam potential energy for this tree with the value of the initial beam energy
				// To consider energy attenuation only between crown parts for potential energy 
				e_pot_incident[v_interc[j]->vectIdTree] -= potential_energy_slope;

				// Apply beer lambert law on incident ray attenuated by the competitor crowns intercepted the crown
				intercepted_energy_slope = this->applyBeerLambert(
					current_energy_slope,
					this->extinctionCoef, this->clumpingFactor,
					crown.getCrownLAD(),
					v_interc[j]->length);

				intercepted_energy_horizontal = this->applyBeerLambert(
					current_energy_horizontal,
					this->extinctionCoef, this->clumpingFactor,
					crown.getCrownLAD(),
					v_interc[j]->length);
				
			}
			// Porous envelop ==> reduce the energy by a fixed amount
			else {
				// Be careful, only once per crown (if many crown part, reduce only if crown has not been already intercepted by the ray)
				if (!is_intercepted[v_interc[j]->vectIdTree]) {

					// Compute potential and intercepted energy if the first crown part is encountered by reducing the incident enrrgy by a fixed amount
					potential_energy_slope = e_incident_slope_cell * (1 - crown.getCrownOpenness());

					intercepted_energy_slope = current_energy_slope * (1 - crown.getCrownOpenness());

					intercepted_energy_horizontal = current_energy_horizontal * (1 - crown.getCrownOpenness());

					// Set the tree as being intercepted
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

			// Set the tree as the first intercepted one if it is
			if (idtree_first_intercepted == -1)
				idtree_first_intercepted = v_interc[j]->vectIdTree;


			// Add to the potential and intercepted energy by the tree
			// CAREFUL: ONLY IF TARGET IS A CELL 
			// If it is a sensor, do not consider energy of trees
			// AND DO NOT ADD IF IT IS AN INTERCEPTION WITH A TRUNK
			if (!target->isThisSensor() && !is_trunk) {
				#ifdef _OPENMP
				#pragma omp critical
				{
				#endif

				// Add energy intercepted real, potential and unobstructed considering the ray is a diffuse or direct one
				// We add the energy on a slope as we consider the tree
				if (ray->isDirect()) {
					crown.addInterceptedRayDirect();
					crown.addEnergyPotentialDirect(potential_energy_slope);
					crown.addEnergyDirect(intercepted_energy_slope);

					if (v_interc[j]->vectIdTree == idtree_first_intercepted) {
						crown.addInterceptedUnobstructedRayDirect();
						crown.addEnergyUnobstructedDirect(potential_energy_slope);
					}
				}
				else {
					crown.addInterceptedRayDiffuse();
					crown.addEnergyPotentialDiffuse(potential_energy_slope);
					crown.addEnergyDiffuse(intercepted_energy_slope);

					if (v_interc[j]->vectIdTree == idtree_first_intercepted) {
						crown.addInterceptedUnobstructedRayDiffuse();
						crown.addEnergyUnobstructedDiffuse(potential_energy_slope);
					}
				}

				#ifdef _OPENMP
				}
				#endif
			}


			// Remove the intercepted energy by the crown from the total energy above the target
			// (it is easier to decrease the energy coming from above along interceptions)
			// (than to add current energy of the ray coming to the target at the end of the interceptions)
			// Modifying the target is not a critical task for parallelisation as we iterate over independent cells/sensors
			// The first argument specify if the ray is unobstructed or not, thus the inverse of the ray has already been intercepted
			// Precise the id in the tree vector that have intercepted the ray, to modify the vectors of interceptions with trees
			if (ray->isDirect()) {
				target->interceptRayDirect(v_interc[j]->vectIdTree == idtree_first_intercepted, 
					intercepted_energy_slope, intercepted_energy_horizontal, v_interc[j]->vectIdTree);
			}
			else {
				target->interceptRayDiffuse(v_interc[j]->vectIdTree == idtree_first_intercepted,
					intercepted_energy_slope, intercepted_energy_horizontal, v_interc[j]->vectIdTree);
			}


			// Remove the intercepted energy from the energy of the ray
			// i.e. output energy is the transmitted energy
			current_energy_slope -= intercepted_energy_slope;
			current_energy_horizontal -= intercepted_energy_horizontal;

			// Delete interception pointer
			delete v_interc[j];
		}

		// Increase the unobstructed energy of the target if the ray has not been intercepted
		if (n_interceptions == 0) {

			if (ray->isDirect()) {
				target->addUnobstructedEnergiesDirect(e_incident_slope_cell, e_incident_horizontal_cell);
			}
			else {
				target->addUnobstructedEnergiesDiffuse(e_incident_slope_cell, e_incident_horizontal_cell);
			}
		}
	}

	double applyBeerLambert(double incident_energy, double extinction_coef, double clumping_factor, double leaf_area_density, double path_length) {

		return(incident_energy * (1 - exp(-extinction_coef * clumping_factor * leaf_area_density * path_length)));
	}

	// Geometry of the rays toward a set of targets (sensors or cells), independent from the crown LAD and openness
	// For each ray x target pair intercepted at least once: the trees whose crown is crossed, the path length across
	// each crown, and whether a trunk intercepts the ray.
	// The order of interceptions is not stored: the energy reaching a target is the product of the transmittances
	// of the crossed crowns (turbid medium or porous envelop), or null if a trunk intercepts the ray (when considered).
	geometryTargets computeGeometryTargets(std::vector<Target*>& targets, std::vector<int>& targets_id,
		std::vector<double>& rays_e) {

		int n_trees = this->stand.getNTrees();
		int n_rays = this->rays.getNRays();
		int n_targets = targets.size();

		// Buffers per target, filled in parallel and merged in the order of the targets (output independent from the number of threads)
		std::vector<std::vector<int>> pairs_ray(n_targets);
		std::vector<std::vector<int>> pairs_blocked(n_targets);
		std::vector<std::vector<int>> pairs_ntrees(n_targets);
		std::vector<std::vector<int>> crossings_pair(n_targets);
		std::vector<std::vector<int>> crossings_tree(n_targets);
		std::vector<std::vector<double>> crossings_length(n_targets);

		geometryTargets geometry;
		geometry.targetsId = targets_id;
		geometry.eUnobsDirect.assign(n_targets, 0.0);
		geometry.eUnobsDiffuse.assign(n_targets, 0.0);

		#ifdef _OPENMP
		#pragma omp parallel for schedule(static)
		#endif
		for (int k = 0; k < n_targets; k++) {

			Target* target = targets[k];

			// Path length across each crossed crown for the current ray
			std::vector<double> length_tree(n_trees, 0.0);
			std::vector<bool> is_crossed(n_trees, false);
			std::vector<int> trees_crossed;

			for (int i = 0; i < n_rays; i++) {

				Ray* ray = this->rays.getRay(i);

				std::vector<interception*> v_interc = this->findInterceptionsRayToTarget(ray, target);

				// Unobstructed ray: its whole energy reaches the target whatever the crown properties
				// (energy on a horizontal plane for sensors, on the slope plane for cells, as in exportResults())
				if (v_interc.empty()) {
					if (this->raysGeometry.direct[i]) {
						geometry.eUnobsDirect[k] += rays_e[i];
					}
					else {
						geometry.eUnobsDiffuse[k] += rays_e[i];
					}
					continue;
				}

				// Sum the path lengths per tree (several crown parts, or several torus copies of the same tree)
				// A crown is counted as crossed even with a null path length (porous envelop attenuates the ray at the first crown part)
				// Crowns are stored even if a trunk intercepts the ray: the trunk rule is applied when computing the energy
				bool blocked = false;
				int n_interceptions = v_interc.size();
				for (int j = 0; j < n_interceptions; j++) {
					if (v_interc[j]->withTrunk) {
						blocked = true;
					}
					else {
						int t = v_interc[j]->vectIdTree;
						if (!is_crossed[t]) {
							is_crossed[t] = true;
							trees_crossed.push_back(t);
						}
						length_tree[t] += v_interc[j]->length;
					}
					delete v_interc[j];
				}

				// Store the pair and its crossed crowns
				int n_crossed = trees_crossed.size();

				pairs_ray[k].push_back(i + 1);
				pairs_blocked[k].push_back(blocked);
				pairs_ntrees[k].push_back(n_crossed);

				// Index of the pair among the pairs of the target (shifted to the index among all the pairs when merging the targets)
				int pair_target = pairs_ray[k].size();

				for (int m = 0; m < n_crossed; m++) {
					int t = trees_crossed[m];
					crossings_pair[k].push_back(pair_target);
					crossings_tree[k].push_back(t + 1);
					crossings_length[k].push_back(length_tree[t]);
					length_tree[t] = 0.0;
					is_crossed[t] = false;
				}
				trees_crossed.clear();
			}
		}

		// Merge the buffers of the targets
		for (int k = 0; k < n_targets; k++) {

			// Number of pairs of the previous targets
			int n_pairs_previous = geometry.pairsTarget.size();
			for (int pair_target : crossings_pair[k]) {
				geometry.crossingsPair.push_back(n_pairs_previous + pair_target);
			}

			geometry.pairsTarget.insert(geometry.pairsTarget.end(), pairs_ray[k].size(), k + 1);
			geometry.pairsRay.insert(geometry.pairsRay.end(), pairs_ray[k].begin(), pairs_ray[k].end());
			geometry.pairsBlocked.insert(geometry.pairsBlocked.end(), pairs_blocked[k].begin(), pairs_blocked[k].end());
			geometry.pairsNTrees.insert(geometry.pairsNTrees.end(), pairs_ntrees[k].begin(), pairs_ntrees[k].end());
			geometry.crossingsTree.insert(geometry.crossingsTree.end(), crossings_tree[k].begin(), crossings_tree[k].end());
			geometry.crossingsLength.insert(geometry.crossingsLength.end(), crossings_length[k].begin(), crossings_length[k].end());
		}

		return(geometry);
	}


public:

	Model(DataFrame trees, 
		DataFrame sensors, bool sensors_only,
		DataFrame cells,
		DataFrame rays, 
		double e_direct_above_slope_m2, double e_diffuse_above_slope_m2,
		double e_direct_above_horizontal_m2, double e_diffuse_above_horizontal_m2,
		double slope, double north_to_x_cw, double aspect,
		double cell_size, double n_cells_x, double n_cells_y,
		bool use_torus, 
		bool turbid_medium, double extinction_coef, double clumping_factor,
		bool trunk_interception) :

		stand(trees, sensors, cells,
			e_direct_above_slope_m2, e_diffuse_above_slope_m2, e_direct_above_horizontal_m2, e_diffuse_above_horizontal_m2,
			slope, north_to_x_cw, aspect, cell_size, n_cells_x, n_cells_y, use_torus),
		rays(rays, e_direct_above_slope_m2, e_diffuse_above_slope_m2, e_direct_above_horizontal_m2, e_diffuse_above_horizontal_m2)
	{
		this->sensorsOnly = sensors_only;
		this->turbidMedium = turbid_medium;
		this->extinctionCoef = extinction_coef;
		this->clumpingFactor = clumping_factor;
		this->trunkInterception = trunk_interception;
	}

	// Getters
	double getMaximumTreeHeight() { return(this->stand.getMaximumTreeHeight()); }
	double getMaximumTreeCrownRadius() { return(this->stand.getMaximumTreeCrownRadius()); }


	// Methods
	void findPotentialRelCells() {
		
		// Potential relative cells for each ray coming to a cell
		// Rays are independent thus we can parallelize them
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

		// For parallelisation, there is a critical moment that is to update the targets and trees energies
		// Indeed, each ray can update the same target, anbd each ray X target can update the same tree
		// Thus, critical moment for target is managed by parallelisation over target, not rays
		// i.e. rays are computed iterativelly inside an independent target (avoiding critical moment for target energies updating)
		// But for trees critical moment, this is managed inside the function using critical flag when updating tree energies inside the interception loop

		int n_rays = this->rays.getNRays();
		int n_sensors = this->stand.getNSensors();
		int n_cells_x = this->stand.getNCellsX();
		int n_cells_y = this->stand.getNCellsY();

		// For each ray toward each target sensor
		if (this->stand.areThereSensors()) {
			#ifdef _OPENMP
			#pragma omp parallel for schedule(static)
			#endif
			for (int s = 0; s < n_sensors; s++) {
				Sensor* sensor = this->stand.getSensor(s);

				for (int i = 0; i < n_rays; i++) {
					Ray* ray = this->rays.getRay(i);
					predictEnergiesRayToTarget(ray, sensor);
				}
			}
		}

		// For each ray toward each target cell
		// run only when specified by the user
		if (!this->sensorsOnly) {
			#ifdef _OPENMP
			#pragma omp parallel for collapse(2) schedule(static)
			#endif
			for (int r = 0; r < n_cells_y; r++) {
				for (int c = 0; c < n_cells_x; c++) {
					Cell* target_cell = this->stand.getCell(r, c);

					for (int i = 0; i < n_rays; i++) {
						Ray* ray = this->rays.getRay(i);
						predictEnergiesRayToTarget(ray, target_cell);
					}
				}
			}
		}
	}

	void computeGeometry() {

		// Geometry of the rays toward the sensors and the cells, independent from the crown LAD and openness (see precompute_sl_geometry())
		// Counterpart of computeInterceptions(): cells are considered only if not sensorsOnly

		// Energy of each ray above canopy on a horizontal plane and on the slope plane (in MJ/m2),
		// as in summarizeInterceptionsRayToTarget(), and type of ray (direct or diffuse)
		int n_rays = this->rays.getNRays();
		this->raysGeometry.eHorizontal.assign(n_rays, 0.0);
		this->raysGeometry.eSlope.assign(n_rays, 0.0);
		this->raysGeometry.direct.assign(n_rays, 0);

		for (int i = 0; i < n_rays; i++) {
			Ray* ray = this->rays.getRay(i);
			double scalar_slope = cos(stand.getSlope()) * ray->getSinHeightAngle() +
				sin(stand.getSlope()) * ray->getCosHeightAngle() * cos(ray->getAzimuth() - stand.getBottomAzimuth());
			double scalar_horizontal = cos(0.0) * ray->getSinHeightAngle() +
				sin(0.0) * ray->getCosHeightAngle() * cos(ray->getAzimuth() - stand.getBottomAzimuth());
			this->raysGeometry.eSlope[i] = scalar_slope * ray->getIncidentEnergy();
			this->raysGeometry.eHorizontal[i] = scalar_horizontal * ray->getIncidentEnergy();
			this->raysGeometry.direct[i] = ray->isDirect();
		}

		// Geometry toward the sensors
		if (this->stand.areThereSensors()) {
			int n_sensors = this->stand.getNSensors();
			std::vector<Target*> sensors;
			std::vector<int> sensors_id;
			for (int s = 0; s < n_sensors; s++) {
				Sensor* sensor = this->stand.getSensor(s);
				sensors.push_back(sensor);
				sensors_id.push_back(sensor->getIdSensor());
			}
			this->sensorsGeometry = this->computeGeometryTargets(sensors, sensors_id, this->raysGeometry.eHorizontal);
		}

		// Geometry toward the cells, row by row
		if (!this->sensorsOnly) {
			std::vector<Target*> cells;
			std::vector<int> cells_id;
			for (int r = 0; r < this->stand.getNCellsY(); r++) {
				for (int c = 0; c < this->stand.getNCellsX(); c++) {
					Cell* cell = this->stand.getCell(r, c);
					cells.push_back(cell);
					cells_id.push_back(cell->getIdCell());
				}
			}
			this->cellsGeometry = this->computeGeometryTargets(cells, cells_id, this->raysGeometry.eSlope);
		}
	}

	List exportGeometry() {

		// Convert the geometry computed by computeGeometry() into R format (only function creating Rcpp objects for the geometry)

		// Convert the geometry toward a set of targets (NULL if not computed)
		auto export_targets = [](geometryTargets& geometry, bool is_computed, std::string name_id) -> SEXP {

			if (!is_computed) { return(R_NilValue); }

			return(List::create(
				Named("targets") = DataFrame::create(
					Named(name_id) = IntegerVector(geometry.targetsId.begin(), geometry.targetsId.end()),
					Named("e_unobs_direct") = NumericVector(geometry.eUnobsDirect.begin(), geometry.eUnobsDirect.end()),
					Named("e_unobs_diffuse") = NumericVector(geometry.eUnobsDiffuse.begin(), geometry.eUnobsDiffuse.end())
				),
				Named("pairs") = DataFrame::create(
					Named("target") = IntegerVector(geometry.pairsTarget.begin(), geometry.pairsTarget.end()),
					Named("ray") = IntegerVector(geometry.pairsRay.begin(), geometry.pairsRay.end()),
					Named("blocked") = LogicalVector(geometry.pairsBlocked.begin(), geometry.pairsBlocked.end()),
					Named("n_trees") = IntegerVector(geometry.pairsNTrees.begin(), geometry.pairsNTrees.end())
				),
				Named("crossings") = DataFrame::create(
					Named("pair") = IntegerVector(geometry.crossingsPair.begin(), geometry.crossingsPair.end()),
					Named("tree") = IntegerVector(geometry.crossingsTree.begin(), geometry.crossingsTree.end()),
					Named("length") = NumericVector(geometry.crossingsLength.begin(), geometry.crossingsLength.end())
				)
			));
		};

		return(List::create(
			Named("sensors") = export_targets(this->sensorsGeometry, this->stand.areThereSensors(), "id_sensor"),
			Named("cells") = export_targets(this->cellsGeometry, !this->sensorsOnly, "id_cell"),
			Named("rays") = DataFrame::create(
				Named("e_horizontal") = NumericVector(this->raysGeometry.eHorizontal.begin(), this->raysGeometry.eHorizontal.end()),
				Named("e_slope") = NumericVector(this->raysGeometry.eSlope.begin(), this->raysGeometry.eSlope.end()),
				Named("direct") = LogicalVector(this->raysGeometry.direct.begin(), this->raysGeometry.direct.end())
			)
		));
	}

	List exportResults() {

		int n_trees = this->stand.getNTrees();
		int n_sensors = this->stand.getNSensors();
		int n_cells = this->stand.getNCells();

		// ------------- FOR SENSORS ----------------

		// Init RCPP vectors/matrices for sensors

		IntegerVector id_sensors(n_sensors);

		NumericVector e_sensors(n_sensors);
		NumericVector pacl_sensors(n_sensors);
		NumericVector punobs_sensors(n_sensors);
		NumericVector e_direct_sensors(n_sensors);
		NumericVector pacl_direct_sensors(n_sensors);
		NumericVector punobs_direct_sensors(n_sensors);
		NumericVector e_diffuse_sensors(n_sensors);
		NumericVector pacl_diffuse_sensors(n_sensors);
		NumericVector punobs_diffuse_sensors(n_sensors);

		IntegerMatrix interceptions_number_sensors(n_trees, n_sensors);
		IntegerMatrix interceptions_number_direct_sensors(n_trees, n_sensors);
		IntegerMatrix interceptions_number_diffuse_sensors(n_trees, n_sensors);

		IntegerMatrix interceptions_number_unobstructed_sensors(n_trees, n_sensors);
		IntegerMatrix interceptions_number_unobstructed_direct_sensors(n_trees, n_sensors);
		IntegerMatrix interceptions_number_unobstructed_diffuse_sensors(n_trees, n_sensors);

		NumericMatrix interceptions_energy_sensors(n_trees, n_sensors);
		NumericMatrix interceptions_energy_direct_sensors(n_trees, n_sensors);
		NumericMatrix interceptions_energy_diffuse_sensors(n_trees, n_sensors);


		// Export the sensor into R
		CharacterVector rcpp_sensor_id_char(n_sensors);
		for (int s = 0; s < n_sensors; s++) {

			// Get sensor
			Sensor* sensor = this->stand.getSensor(s);

			// Update sensor character id vector
			rcpp_sensor_id_char[s] = std::to_string(sensor->getIdSensor());

			// Export sensor info to vectors
			id_sensors[s] = sensor->getIdSensor();

			// Correct for rounding errors
			/*sensor->correctNullEnergy(1e-6);*/

			// Export energies toward the sensor
			// BE CAREFUL !! Horizontal energy, in MJ/m2
			e_sensors[s] = sensor->getEnergyHorizontal() / this->stand.getCellAreaHorizontal();
			e_direct_sensors[s] = sensor->getEnergyDirectHorizontal() / this->stand.getCellAreaHorizontal();
			e_diffuse_sensors[s] = sensor->getEnergyDiffuseHorizontal() / this->stand.getCellAreaHorizontal();

			pacl_sensors[s] = e_sensors[s] / this->rays.getEnergyAboveHorizontalM2();
			pacl_direct_sensors[s] = e_direct_sensors[s] / this->rays.getEnergyDirectAboveHorizontalM2();
			pacl_diffuse_sensors[s] = e_diffuse_sensors[s] / this->rays.getEnergyDiffuseAboveHorizontalM2();

			punobs_sensors[s] = sensor->getEnergyUnobstructedHorizontal() / sensor->getEnergyHorizontal();
			punobs_direct_sensors[s] = sensor->getEnergyUnobstructedDirectHorizontal() / sensor->getEnergyDirectHorizontal();
			punobs_diffuse_sensors[s] = sensor->getEnergyUnobstructedDiffuseHorizontal() / sensor->getEnergyDiffuseHorizontal();


			// Export matrices of interceptions with trees
			std::vector<int>& interceptions_number_sensor = sensor->getInterceptionsTreesNumber();
			interceptions_number_sensors(_, s) = NumericVector(interceptions_number_sensor.begin(), interceptions_number_sensor.end());
			
			std::vector<int>& interceptions_number_direct_sensor = sensor->getInterceptionsTreesNumberDirect();
			interceptions_number_direct_sensors(_, s) = NumericVector(interceptions_number_direct_sensor.begin(), interceptions_number_direct_sensor.end());
			
			std::vector<int>& interceptions_number_diffuse_sensor = sensor->getInterceptionsTreesNumberDiffuse();
			interceptions_number_diffuse_sensors(_, s) = NumericVector(interceptions_number_diffuse_sensor.begin(), interceptions_number_diffuse_sensor.end());


			std::vector<int>& interceptions_number_unobstructed_sensor = sensor->getInterceptionsTreesNumberUnobstructed();
			interceptions_number_unobstructed_sensors(_, s) = NumericVector(interceptions_number_unobstructed_sensor.begin(), interceptions_number_unobstructed_sensor.end());

			std::vector<int>& interceptions_number_unobstructed_direct_sensor = sensor->getInterceptionsTreesNumberUnobstructedDirect();
			interceptions_number_unobstructed_direct_sensors(_, s) = NumericVector(interceptions_number_unobstructed_direct_sensor.begin(), interceptions_number_unobstructed_direct_sensor.end());

			std::vector<int>& interceptions_number_unobstructed_diffuse_sensor = sensor->getInterceptionsTreesNumberUnobstructedDiffuse();
			interceptions_number_unobstructed_diffuse_sensors(_, s) = NumericVector(interceptions_number_unobstructed_diffuse_sensor.begin(), interceptions_number_unobstructed_diffuse_sensor.end());


			std::vector<double>& interceptions_energy_sensor = sensor->getInterceptionsTreesEnergyHorizontal();
			interceptions_energy_sensors(_, s) = NumericVector(interceptions_energy_sensor.begin(), interceptions_energy_sensor.end());

			std::vector<double>& interceptions_energy_direct_sensor = sensor->getInterceptionsTreesEnergyDirectHorizontal();
			interceptions_energy_direct_sensors(_, s) = NumericVector(interceptions_energy_direct_sensor.begin(), interceptions_energy_direct_sensor.end());

			std::vector<double>& interceptions_energy_diffuse_sensor = sensor->getInterceptionsTreesEnergyDiffuseHorizontal();
			interceptions_energy_diffuse_sensors(_, s) = NumericVector(interceptions_energy_diffuse_sensor.begin(), interceptions_energy_diffuse_sensor.end());
		}

		// Set the col and row names of the matrices
		rownames(interceptions_number_sensors) = this->stand.getRcppTreeIdChar();
		rownames(interceptions_number_direct_sensors) = this->stand.getRcppTreeIdChar();
		rownames(interceptions_number_diffuse_sensors) = this->stand.getRcppTreeIdChar();
		rownames(interceptions_number_unobstructed_sensors) = this->stand.getRcppTreeIdChar();
		rownames(interceptions_number_unobstructed_direct_sensors) = this->stand.getRcppTreeIdChar();
		rownames(interceptions_number_unobstructed_diffuse_sensors) = this->stand.getRcppTreeIdChar();
		rownames(interceptions_energy_sensors) = this->stand.getRcppTreeIdChar();
		rownames(interceptions_energy_direct_sensors) = this->stand.getRcppTreeIdChar();
		rownames(interceptions_energy_diffuse_sensors) = this->stand.getRcppTreeIdChar();

		colnames(interceptions_number_sensors) = rcpp_sensor_id_char;
		colnames(interceptions_number_direct_sensors) = rcpp_sensor_id_char;
		colnames(interceptions_number_diffuse_sensors) = rcpp_sensor_id_char;
		colnames(interceptions_number_unobstructed_sensors) = rcpp_sensor_id_char;
		colnames(interceptions_number_unobstructed_direct_sensors) = rcpp_sensor_id_char;
		colnames(interceptions_number_unobstructed_diffuse_sensors) = rcpp_sensor_id_char;
		colnames(interceptions_energy_sensors) = rcpp_sensor_id_char;
		colnames(interceptions_energy_direct_sensors) = rcpp_sensor_id_char;
		colnames(interceptions_energy_diffuse_sensors) = rcpp_sensor_id_char;
		

		// Create sensors RCPP DataFrames
		DataFrame output_sensors = DataFrame::create(
			Named("id_sensor") = id_sensors,

			Named("e") = e_sensors,
			Named("e_direct") = e_direct_sensors,
			Named("e_diffuse") = e_diffuse_sensors,

			Named("pacl") = pacl_sensors,
			Named("pacl_direct") = pacl_direct_sensors,
			Named("pacl_diffuse") = pacl_diffuse_sensors,

			Named("punobs") = punobs_sensors,
			Named("punobs_direct") = punobs_direct_sensors,
			Named("punobs_diffuse") = punobs_diffuse_sensors
		);


		// -------------- FOR TREES AND CELLS ------------------

		// Get number of trees and cells
		// Set to 0 if we did not compute trees and cells interception (sensor_only argument)
		if (this->sensorsOnly) {
			n_cells = 0;
			n_trees = 0;
		}

		// Init RCPP vectors for trees
		IntegerVector id_trees(n_trees);

		NumericVector e_trees(n_trees);
		NumericVector e_direct_trees(n_trees);
		NumericVector e_diffuse_trees(n_trees);

		NumericVector epot_trees(n_trees);
		NumericVector epot_direct_trees(n_trees);
		NumericVector epot_diffuse_trees(n_trees);

		NumericVector lci_trees(n_trees);
		NumericVector lci_direct_trees(n_trees);
		NumericVector lci_diffuse_trees(n_trees);

		NumericVector eunobs_trees(n_trees);
		NumericVector eunobs_direct_trees(n_trees);
		NumericVector eunobs_diffuse_trees(n_trees);

		NumericVector rci_trees(n_trees);
		NumericVector rci_direct_trees(n_trees);
		NumericVector rci_diffuse_trees(n_trees);


		// Init RCPP vectors/matrices for cells
		IntegerVector id_cells(n_cells);

		NumericVector e_cells(n_cells);
		NumericVector e_direct_cells(n_cells);
		NumericVector e_diffuse_cells(n_cells);

		NumericVector pacl_cells(n_cells);
		NumericVector pacl_direct_cells(n_cells);
		NumericVector pacl_diffuse_cells(n_cells);

		NumericVector punobs_cells(n_cells);
		NumericVector punobs_direct_cells(n_cells);
		NumericVector punobs_diffuse_cells(n_cells);

		IntegerMatrix interceptions_number_cells(n_trees, n_cells);
		IntegerMatrix interceptions_number_direct_cells(n_trees, n_cells);
		IntegerMatrix interceptions_number_diffuse_cells(n_trees, n_cells);

		IntegerMatrix interceptions_number_unobstructed_cells(n_trees, n_cells);
		IntegerMatrix interceptions_number_unobstructed_direct_cells(n_trees, n_cells);
		IntegerMatrix interceptions_number_unobstructed_diffuse_cells(n_trees, n_cells);

		NumericMatrix interceptions_energy_cells(n_trees, n_cells);
		NumericMatrix interceptions_energy_direct_cells(n_trees, n_cells);
		NumericMatrix interceptions_energy_diffuse_cells(n_trees, n_cells);



		if (!this->sensorsOnly) {

			// For each cell
			CharacterVector rcpp_cell_id_char(n_cells);
			
			int icell = 0;
			int itree = 0;
			int n_rows = this->stand.getNCellsY();
			int n_cols = this->stand.getNCellsX();

			for (int r = 0; r < n_rows; r++) {
				for (int c = 0; c < n_cols; c++) {

					// Get cell
					Cell* cell = this->stand.getCell(r, c);

					// Update sensor character id vector
					rcpp_cell_id_char[icell] = std::to_string(cell->getIdCell());

					// Add cell to vectors
					id_cells[icell] = cell->getIdCell();

					// Correct for rounding errors
					//cell->correctNullEnergy(1e-6);

					// Export energies for cells
					// BE CAREFUL !! Slope energies, in MJ/m2
					e_cells[icell] = cell->getEnergySlope() / this->stand.getCellAreaSlope();
					e_direct_cells[icell] = cell->getEnergyDirectSlope() / this->stand.getCellAreaSlope();
					e_diffuse_cells[icell] = cell->getEnergyDiffuseSlope() / this->stand.getCellAreaSlope();

					pacl_cells[icell] = e_cells[icell] / this->rays.getEnergyAboveSlopeM2();
					pacl_direct_cells[icell] = e_direct_cells[icell] / this->rays.getEnergyDirectAboveSlopeM2();
					pacl_diffuse_cells[icell] = e_diffuse_cells[icell] / this->rays.getEnergyDiffuseAboveSlopeM2();

					punobs_cells[icell] = cell->getEnergyUnobstructedSlope() / cell->getEnergySlope();
					punobs_direct_cells[icell] = cell->getEnergyUnobstructedDirectSlope() / cell->getEnergyDirectSlope();
					punobs_diffuse_cells[icell] = cell->getEnergyUnobstructedDiffuseSlope() / cell->getEnergyDiffuseSlope();


					// Export matrices of interceptions with trees
					std::vector<int>& interceptions_number_cell = cell->getInterceptionsTreesNumber();
					interceptions_number_cells(_, icell) = NumericVector(interceptions_number_cell.begin(), interceptions_number_cell.end());

					std::vector<int>& interceptions_number_direct_cell = cell->getInterceptionsTreesNumberDirect();
					interceptions_number_direct_cells(_, icell) = NumericVector(interceptions_number_direct_cell.begin(), interceptions_number_direct_cell.end());

					std::vector<int>& interceptions_number_diffuse_cell = cell->getInterceptionsTreesNumberDiffuse();
					interceptions_number_diffuse_cells(_, icell) = NumericVector(interceptions_number_diffuse_cell.begin(), interceptions_number_diffuse_cell.end());


					std::vector<int>& interceptions_number_unobstructed_cell = cell->getInterceptionsTreesNumberUnobstructed();
					interceptions_number_unobstructed_cells(_, icell) = NumericVector(interceptions_number_unobstructed_cell.begin(), interceptions_number_unobstructed_cell.end());

					std::vector<int>& interceptions_number_unobstructed_direct_cell = cell->getInterceptionsTreesNumberUnobstructedDirect();
					interceptions_number_unobstructed_direct_cells(_, icell) = NumericVector(interceptions_number_unobstructed_direct_cell.begin(), interceptions_number_unobstructed_direct_cell.end());

					std::vector<int>& interceptions_number_unobstructed_diffuse_cell = cell->getInterceptionsTreesNumberUnobstructedDiffuse();
					interceptions_number_unobstructed_diffuse_cells(_, icell) = NumericVector(interceptions_number_unobstructed_diffuse_cell.begin(), interceptions_number_unobstructed_diffuse_cell.end());


					std::vector<double>& interceptions_energy_cell = cell->getInterceptionsTreesEnergySlope();
					interceptions_energy_cells(_, icell) = NumericVector(interceptions_energy_cell.begin(), interceptions_energy_cell.end());

					std::vector<double>& interceptions_energy_direct_cell = cell->getInterceptionsTreesEnergyDirectSlope();
					interceptions_energy_direct_cells(_, icell) = NumericVector(interceptions_energy_direct_cell.begin(), interceptions_energy_direct_cell.end());

					std::vector<double>& interceptions_energy_diffuse_cell = cell->getInterceptionsTreesEnergyDiffuseSlope();
					interceptions_energy_diffuse_cells(_, icell) = NumericVector(interceptions_energy_diffuse_cell.begin(), interceptions_energy_diffuse_cell.end());

				

					// For each tree composing the cell
					int n_trees_cell = cell->getNTrees();
					for (int t = 0; t < n_trees_cell; t++) {

						Tree* tree = this->stand.getTree(cell->getVectIdTree(t));

						id_trees[itree] = tree->getId();

						e_trees[itree] = tree->getCrownEnergy();
						e_direct_trees[itree] = tree->getCrownEnergyDirect();
						e_diffuse_trees[itree] = tree->getCrownEnergyDiffuse();

						epot_trees[itree] = tree->getCrownEnergyPotential();
						epot_direct_trees[itree] = tree->getCrownEnergyPotentialDirect();
						epot_diffuse_trees[itree] = tree->getCrownEnergyPotentialDiffuse();

						lci_trees[itree] = 1.0f - e_trees[itree] / epot_trees[itree];
						lci_direct_trees[itree] = 1.0f - e_direct_trees[itree] / epot_direct_trees[itree];
						lci_diffuse_trees[itree] = 1.0f - e_diffuse_trees[itree] / epot_diffuse_trees[itree];

						eunobs_trees[itree] = tree->getCrownEnergyUnobstructed();
						eunobs_direct_trees[itree] = tree->getCrownEnergyUnobstructedDirect();
						eunobs_diffuse_trees[itree] = tree->getCrownEnergyUnobstructedDiffuse();

						rci_trees[itree] = tree->getCrownRci();
						rci_direct_trees[itree] = tree->getCrownRciDirect();
						rci_diffuse_trees[itree] = tree->getCrownRciDiffuse();

						itree++;
					}

					icell++;
				}
			}

			// Set the col and row names of the matrices
			rownames(interceptions_number_cells) = this->stand.getRcppTreeIdChar();
			rownames(interceptions_number_direct_cells) = this->stand.getRcppTreeIdChar();
			rownames(interceptions_number_diffuse_cells) = this->stand.getRcppTreeIdChar();
			rownames(interceptions_number_unobstructed_cells) = this->stand.getRcppTreeIdChar();
			rownames(interceptions_number_unobstructed_direct_cells) = this->stand.getRcppTreeIdChar();
			rownames(interceptions_number_unobstructed_diffuse_cells) = this->stand.getRcppTreeIdChar();
			rownames(interceptions_energy_cells) = this->stand.getRcppTreeIdChar();
			rownames(interceptions_energy_direct_cells) = this->stand.getRcppTreeIdChar();
			rownames(interceptions_energy_diffuse_cells) = this->stand.getRcppTreeIdChar();

			colnames(interceptions_number_cells) = rcpp_cell_id_char;
			colnames(interceptions_number_direct_cells) = rcpp_cell_id_char;
			colnames(interceptions_number_diffuse_cells) = rcpp_cell_id_char;
			colnames(interceptions_number_unobstructed_cells) = rcpp_cell_id_char;
			colnames(interceptions_number_unobstructed_direct_cells) = rcpp_cell_id_char;
			colnames(interceptions_number_unobstructed_diffuse_cells) = rcpp_cell_id_char;
			colnames(interceptions_energy_cells) = rcpp_cell_id_char;
			colnames(interceptions_energy_direct_cells) = rcpp_cell_id_char;
			colnames(interceptions_energy_diffuse_cells) = rcpp_cell_id_char;

		}

		// Create trees and cells RCPP DataFrames
		DataFrame output_trees = DataFrame::create(
			Named("id_tree") = id_trees,

			Named("epot") = epot_trees,
			Named("epot_direct") = epot_direct_trees,
			Named("epot_diffuse") = epot_diffuse_trees,

			Named("e") = e_trees,
			Named("e_direct") = e_direct_trees,
			Named("e_diffuse") = e_diffuse_trees,

			Named("lci") = lci_trees,
			Named("lci_direct") = lci_direct_trees,
			Named("lci_diffuse") = lci_diffuse_trees,

			Named("eunobs") = eunobs_trees,
			Named("eunobs_direct") = eunobs_direct_trees,
			Named("eunobs_diffuse") = eunobs_diffuse_trees,

			Named("rci") = rci_trees,
			Named("rci_direct") = rci_direct_trees,
			Named("rci_diffuse") = rci_diffuse_trees
		);

		DataFrame output_cells = DataFrame::create(
			Named("id_cell") = id_cells,

			Named("e") = e_cells,
			Named("e_direct") = e_direct_cells,
			Named("e_diffuse") = e_diffuse_cells,

			Named("pacl") = pacl_cells,
			Named("pacl_direct") = pacl_direct_cells,
			Named("pacl_diffuse") = pacl_diffuse_cells,

			Named("punobs") = punobs_cells,
			Named("punobs_direct") = punobs_direct_cells,
			Named("punobs_diffuse") = punobs_diffuse_cells
		);

		// Return output as a List of two DataFrames
		return(List::create(
			Named("sensors") = output_sensors,
			Named("trees") = output_trees,
			Named("cells") = output_cells,
			Named("interceptions") = List::create(
				Named("sensors") = List::create(
					Named("number") = List::create(
						Named("total") = interceptions_number_sensors,
						Named("direct") = interceptions_number_direct_sensors,
						Named("diffuse") = interceptions_number_diffuse_sensors
					),
					Named("number_unobstructed") = List::create(
						Named("total") = interceptions_number_unobstructed_sensors,
						Named("direct") = interceptions_number_unobstructed_direct_sensors,
						Named("diffuse") = interceptions_number_unobstructed_diffuse_sensors
					),
					Named("energy") = List::create(
						Named("total") = interceptions_energy_sensors,
						Named("direct") = interceptions_energy_direct_sensors,
						Named("diffuse") = interceptions_energy_diffuse_sensors
					)
				),
				Named("cells") = List::create(
					Named("number") = List::create(
						Named("total") = interceptions_number_cells,
						Named("direct") = interceptions_number_direct_cells,
						Named("diffuse") = interceptions_number_diffuse_cells
					),
					Named("number_unobstructed") = List::create(
						Named("total") = interceptions_number_unobstructed_cells,
						Named("direct") = interceptions_number_unobstructed_direct_cells,
						Named("diffuse") = interceptions_number_unobstructed_diffuse_cells
					),
					Named("energy") = List::create(
						Named("total") = interceptions_energy_cells,
						Named("direct") = interceptions_energy_direct_cells,
						Named("diffuse") = interceptions_energy_diffuse_cells
					)
				)
			)
		));
	}

};


// [[Rcpp::export]]
List run_sl_rcpp(
	DataFrame trees,
	DataFrame sensors, bool sensors_only,
	DataFrame cells,
	DataFrame rays, 
	double e_direct_above_slope_m2, double e_diffuse_above_slope_m2,
	double e_direct_above_horizontal_m2, double e_diffuse_above_horizontal_m2,
	double slope, double north_to_x_cw, double aspect,
	double cell_size, double n_cells_x, double n_cells_y,
	bool use_torus,
	bool turbid_medium, double extinction_coef, double clumping_factor,
	bool trunk_interception
)
{
	// Reset the counters of unexpected interceptions
	g_n_unexpected_ellipsoid = 0;
	g_n_unexpected_paraboloid = 0;
	g_n_unexpected_trunk = 0;

	// Initialize the model
	Model sl_model = Model(
		trees,
		sensors, sensors_only,
		cells,
		rays, 
		e_direct_above_slope_m2, e_diffuse_above_slope_m2,
		e_direct_above_horizontal_m2, e_diffuse_above_horizontal_m2,
		slope, north_to_x_cw, aspect,
		cell_size, n_cells_x, n_cells_y, 
		use_torus, 
		turbid_medium, extinction_coef, clumping_factor,
		trunk_interception
	);

	// [OPTIMIZATION]: find the extend of possible interception of each ray
	// i.e. for each ray coming to an unknown target cell, find relative cells around the target containing tree that could possibly intercept the ray
	// We can whether or not consider a torus system
	sl_model.findPotentialRelCells();

	// For each target cell and each ray (target cells can be parallelized)
	sl_model.computeInterceptions();

	// Convert the output into R format
	List output = sl_model.exportResults();

	// Add the number of interceptions ignored because of an unexpected number of interception points
	output.push_back(NumericVector::create(
		_["ellipsoid"] = static_cast<double>(g_n_unexpected_ellipsoid.load()),
		_["paraboloid"] = static_cast<double>(g_n_unexpected_paraboloid.load()),
		_["trunk"] = static_cast<double>(g_n_unexpected_trunk.load())),
		"n_unexpected_interceptions");

	return(output);
}


// [[Rcpp::export]]
List precompute_sl_geometry_rcpp(
	DataFrame trees,
	DataFrame sensors, bool sensors_only,
	DataFrame cells,
	DataFrame rays,
	double e_direct_above_slope_m2, double e_diffuse_above_slope_m2,
	double e_direct_above_horizontal_m2, double e_diffuse_above_horizontal_m2,
	double slope, double north_to_x_cw, double aspect,
	double cell_size, double n_cells_x, double n_cells_y,
	bool use_torus
)
{
	// Reset the counters of unexpected interceptions
	g_n_unexpected_ellipsoid = 0;
	g_n_unexpected_paraboloid = 0;
	g_n_unexpected_trunk = 0;

	// Initialize the model
	// The interception model (turbid medium, extinction and clumping) is not used by the geometry
	// Trunks are always considered: whether they intercept the rays is chosen in estimate_pacl_from_geometry()
	Model sl_model = Model(
		trees,
		sensors, sensors_only,
		cells,
		rays,
		e_direct_above_slope_m2, e_diffuse_above_slope_m2,
		e_direct_above_horizontal_m2, e_diffuse_above_horizontal_m2,
		slope, north_to_x_cw, aspect,
		cell_size, n_cells_x, n_cells_y,
		use_torus,
		true, 0.5, 1.0,
		true
	);

	// [OPTIMIZATION]: find the extend of possible interception of each ray (as in run_sl_rcpp())
	sl_model.findPotentialRelCells();

	// Crowns crossed by each ray toward each target
	sl_model.computeGeometry();

	// Convert the geometry into R format
	List output = sl_model.exportGeometry();

	// Add the number of interceptions ignored because of an unexpected number of interception points
	output.push_back(NumericVector::create(
		_["ellipsoid"] = static_cast<double>(g_n_unexpected_ellipsoid.load()),
		_["paraboloid"] = static_cast<double>(g_n_unexpected_paraboloid.load()),
		_["trunk"] = static_cast<double>(g_n_unexpected_trunk.load())),
		"n_unexpected_interceptions");

	return(output);
}


// Energies above each target from a precomputed geometry (see precompute_sl_geometry())
// Pure C++, reading the R vectors through pointers (no copy of the geometry at each call)
// tree_coef: extinction coefficient x clumping factor x LAD (turbid medium) or crown openness (porous envelop), per tree
// Returns a column-major matrix n_targets x 4 (in MJ/m2 on a horizontal plane):
// energy intercepted from direct rays, from diffuse rays, and energy of unobstructed direct rays, diffuse rays
std::vector<double> computeTargetsEnergiesFromGeometry(
	int n_targets, const double* targets_e_unobs_direct, const double* targets_e_unobs_diffuse,
	int n_pairs, const int* pairs_target, const int* pairs_ray, const int* pairs_blocked, const int* pairs_ntrees,
	int n_crossings, const int* crossings_pair, const int* crossings_tree, const double* crossings_length,
	const double* rays_e_horizontal, const int* rays_direct,
	const double* tree_coef, bool turbid_medium, bool trunk_interception)
{
	// Attenuation of each pair by its crossed crowns
	// Turbid medium: optical depth (sum of extinction x clumping x LAD x path length), transmittance = exp(-depth)
	// Porous envelop: transmittance = product of the openness of the crossed crowns (once per crown)
	std::vector<double> attenuation(n_pairs, turbid_medium ? 0.0 : 1.0);

	for (int c = 0; c < n_crossings; c++) {
		int p = crossings_pair[c] - 1;
		double coef = tree_coef[crossings_tree[c] - 1];
		if (turbid_medium) {
			attenuation[p] += coef * crossings_length[c];
		}
		else {
			attenuation[p] *= coef;
		}
	}

	// Energies of the targets: intercepted (columns 0 and 1), unobstructed (columns 2 and 3)
	std::vector<double> energies(n_targets * 4, 0.0);
	for (int k = 0; k < n_targets; k++) {
		energies[2 * n_targets + k] = targets_e_unobs_direct[k];
		energies[3 * n_targets + k] = targets_e_unobs_diffuse[k];
	}

	for (int p = 0; p < n_pairs; p++) {

		int k = pairs_target[p] - 1;
		int r = pairs_ray[p] - 1;
		int col = rays_direct[r] ? 0 : 1;

		// Ray intercepted by a trunk (if considered): no energy transmitted
		if (trunk_interception && pairs_blocked[p]) {
			energies[col * n_targets + k] += rays_e_horizontal[r];
		}
		// Ray intercepted only by trunks, not considered: unobstructed ray (as in run_sl_rcpp() without trunks)
		else if (pairs_ntrees[p] == 0) {
			energies[(col + 2) * n_targets + k] += rays_e_horizontal[r];
		}
		// Ray attenuated by the crossed crowns
		else {
			double transmittance = turbid_medium ? exp(-attenuation[p]) : attenuation[p];
			energies[col * n_targets + k] += rays_e_horizontal[r] * (1.0 - transmittance);
		}
	}

	return(energies);
}


// [[Rcpp::export]]
NumericMatrix estimate_pacl_from_geometry_rcpp(
	NumericVector targets_e_unobs_direct, NumericVector targets_e_unobs_diffuse,
	IntegerVector pairs_target, IntegerVector pairs_ray, LogicalVector pairs_blocked, IntegerVector pairs_ntrees,
	IntegerVector crossings_pair, IntegerVector crossings_tree, NumericVector crossings_length,
	NumericVector rays_e_horizontal, LogicalVector rays_direct,
	NumericVector tree_coef, bool turbid_medium, bool trunk_interception
)
{
	int n_targets = targets_e_unobs_direct.size();

	// Compute the energies without copying the geometry (R logical vectors are stored as integers)
	std::vector<double> energies = computeTargetsEnergiesFromGeometry(
		n_targets, targets_e_unobs_direct.begin(), targets_e_unobs_diffuse.begin(),
		pairs_target.size(), pairs_target.begin(), pairs_ray.begin(), pairs_blocked.begin(), pairs_ntrees.begin(),
		crossings_pair.size(), crossings_pair.begin(), crossings_tree.begin(), crossings_length.begin(),
		rays_e_horizontal.begin(), rays_direct.begin(),
		tree_coef.begin(), turbid_medium, trunk_interception);

	// Convert into R format
	NumericMatrix output(n_targets, 4, energies.begin());
	colnames(output) = CharacterVector::create("e_intercepted_direct", "e_intercepted_diffuse", "e_unobs_direct", "e_unobs_diffuse");

	return(output);
}



// [[Rcpp::export]]
void sl_set_openmp(bool parallel_mode, int num_threads = -1, bool verbose = true) {

	#ifdef _OPENMP
	omp_set_dynamic(0);   // absolutely critical

	if (!parallel_mode) {
		omp_set_num_threads(1);
		g_parallel_enabled = false;
		g_num_threads = 1;
		g_threads_auto = false;
	}
	else {
		int nt;
		if (num_threads < 1) {
			nt = omp_get_num_procs();   // physical cores
			g_threads_auto = true;
		}
		else {
			nt = num_threads;
			g_threads_auto = false;
		}

		omp_set_num_threads(nt);
		g_parallel_enabled = true;
		g_num_threads = nt;
	}

	#else
	if (parallel_mode && verbose) {
		Rcpp::warning("OpenMP not available: running sequentially.");
	}
	g_parallel_enabled = false;
	g_num_threads = 1;
	g_threads_auto = false;
	#endif
}


// [[Rcpp::export]]
void sl_print_openmp_status() {

	#ifdef _OPENMP
	// Number of threads used by the next parallel region (dynamic adjustment disabled in sl_set_openmp())
	int runtime_threads = omp_get_max_threads();

	Rcpp::Rcout <<
		"parallel mode: "
		<< (g_parallel_enabled ? "ENABLED" : "DISABLED")
		<< " | Threads = " << runtime_threads
		<< (g_threads_auto ? " (auto)" : " (user)")
		<< " | Cores available = " << omp_get_num_procs()
		<< "\n";

	#else
	Rcpp::Rcout << "parallel mode disabled because OpenMP was not available\n";
	#endif
}

