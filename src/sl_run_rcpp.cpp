#define _USE_MATH_DEFINES

constexpr auto PARABOLOID = 1;
constexpr auto ELLIPSOID = 2;

#include <stack>
#include <vector>
#include <map>
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

class Tree {

private:
	// Id of the tree
	int id;

	// Tree position
	double x;
	double y;
	double z;

	// Tree crown characteristics
	double height; // in meters
	double crownRadius; // in meters
	double crownBaseHeight; // in meters
	double crownType; // either P for paraboloid or E for ellipsoid

	// Leaves characteristics
	double crownOpenness; // between 0 and 1, percentage of ligth energy of a ray that is not absorbed by the crown (for porous envelop method)
	double crownLAD; // Leaf area density in m2/m3 (for turbid medium method)
	double extinctionCoef; // Probability of a leaf to intercept the ray (linked to leaf orientation)
	double clumpingFactor; // Aggregation of leaves within the crown volume (1 is homogeneous)

	// Output energy (in MJ)
	double energy; // Intercepted energy
	double energyPotential; // Intercepted energy without considering neighbours


public:
	Tree(int id, double x, double y, double z,
		double h, double cr, double cbh, double c_type, double c_p, double c_lad,
		double extinction_coef, double clumping_factor)
	{
		this->id = id;
		this->x = x;
		this->y = y;
		this->z = z;
		this->height = h;
		this->crownRadius = cr;
		this->crownBaseHeight = cbh;
		this->crownType = c_type;
		this->crownOpenness = c_p;
		this->crownLAD = c_lad;
		this->extinctionCoef = extinction_coef;
		this->clumpingFactor = clumping_factor;

		this->energy = 0.0;
		this->energyPotential = 0.0;
	}

	// Getters
	int getId() { return(this->id); }
	double getX() { return(this->x); }
	double getY() { return(this->y); }
	double getZ() { return(this->z); }
	double getHeight() { return(this->height); }
	double getCrownRadius() { return(this->crownRadius); }
	double getCrownBaseHeight() { return(this->crownBaseHeight); }
	double getCrownType() { return(this->crownType); }
	double getCrownOpenness() { return(this->crownOpenness); }
	double getCrownLAD() { return(this->crownLAD); }
	double getExtinctionCoef() { return(this->extinctionCoef); }
	double getClumpingFactor() { return(this->clumpingFactor); }
	double getEnergy() { return(this->energy); }
	double getEnergyPotential() { return(this->energyPotential); }
	double getLci() { return(this->energy > 0.0 ? 1 - this->energy / this->energyPotential : 1); }

	// Setter
	void addEnergy(double e) { this->energy += e; }
	void addEnergyPotential(double epot) { this->energyPotential += epot; }
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
	std::map<int, Cell*> cells;

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
		double e_above_m2,
		double slope, double north_to_x_cw, double aspect,
		double cell_size, double n_cells, bool use_torus)
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

		// Compute energy above the canopy arriving to a target cell
		double e_above_cell = e_above_m2 * this->cellArea;

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
		NumericVector height = trees["height_m"];
		NumericVector cradius = trees["cradius_m"];
		NumericVector cbh = trees["cbh_m"];
		IntegerVector ctype = trees["crown_type"];
		NumericVector cp = trees["crown_p"];
		NumericVector clad = trees["crown_lad"];

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
				height[i], cradius[i], cbh[i], ctype[i], cp[i], clad[i],
				0.5, 1.0);
			this->cells.find(trees_id_cell[i] - 1)->second->addTree(tree);
		}

		// Get maximum height and crown radius of the trees
		this->maxTreeHeight = max(height);
		this->maxTreeCrownRadius = max(cradius);
	}

	~Stand() {
		// Delete cells pointers
		int n_cells = this->cells.size();
		for (int i = 0; i < n_cells; i++) {
			delete this->cells[i];
		}
	}

	// Getters
	Cell* getCell(int row, int col) { return(this->grid[std::make_pair(row, col)]); }
	std::map<int, Cell*> getCells() { return(this->cells); }
	int getNTrees() { return(this->nTrees); }
	double getSlope() { return(this->slope); }
	double getNorthToXcw() { return(this->northToXcw); }
	double getAspect() { return(this->aspect); }
	double getBottomAzimuth() { return(this->bottomAzimuth); }
	double getCellSize() { return(this->cellSize); }
	int getNCells() { return(this->nCells); }
	double getCellAreaPlane() { return(this->cellAreaPlane); }
	double getCellArea() { return(this->cellArea); }
	double isUseTorus() { return(this->useTorus); }
	double getMaximumTreeCrownRadius() { return(this->maxTreeCrownRadius); }
	double getMaximumTreeHeight() { return(this->maxTreeHeight); }
	void printCellKeys() {
		for (auto const& pair : this->cells) {
			std::cout << pair.first << "\t";
		}
	}
};

struct ShiftedCell {

	// Pointor toward the original Cell
	Cell* cell;

	// Shift applied to the cell position (in a torus system)
	double xShift;
	double yShift;
	double zShift;

	// Is the cell outside the main plot ? (torus system)
	bool outside;

};

struct RelativeCoords {
	int row;
	int col;
};

class Interception {

private:
	// Target cell the ray is directed to
	Cell* targetCell;

	// Pointors on the tree intercepted by the ray
	Tree* tree;

	// Length of the path across the crown
	double length;

	// Distance between interception point (middle of full path) and target cell
	double distance;

	// Does the intercepted tree is in the main plot (torus system)
	bool outside;



public:
	Interception(Cell* target_cell, Tree* tree, double length, double distance, bool outside) {
		this->targetCell = target_cell;
		this->tree = tree;
		this->length = length;
		this->distance = distance;
		this->outside = outside;
	}

	Cell* getTargetCell() const { return(this->targetCell); }
	Tree* getTree() const { return(this->tree); }
	double getLength() const { return(this->length); }
	double getDistance() const { return(this->distance); }
	bool isOutside() const { return(this->outside); }


	void print() const {
		std::cout << "Ray toward " << this->targetCell->getId() << " intercepting tree " << this->tree->getId() << "==> l=" << this->length << "; d=" << this->distance << std::endl;
	}

	bool operator > (const Interception& interc) const
	{
		// First, sort by if of the target cellZ
		// Second, sort by distance to the target cell
		return this->targetCell->getId() < interc.getTargetCell()->getId() || (
			!(interc.getTargetCell()->getId() < this->targetCell->getId()) && this->distance > interc.getDistance()
			);
	}
};

class Ray {

private:
	//Id of the ray
	int id;

	// Relative coordinates of cells containing trees potentially crossed by the ray (relatively to an unknown target cell)
	std::vector<RelativeCoords> potentialRelCells;

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
	std::vector<Interception> interceptions;


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
	double getAzimuth() { return(this->azimuth); }
	int getNPotCells() { return(this->potentialRelCells.size()); }


	// Methods
	void findPotentialCells(Stand& stand) {

		// Computes lateral = the boundary to add to the competition rectangle to take into account cells center instead of trees position.
		// The boundary depends on beam azimut.
		double azt = 0;
		if (this->azimuth < M_PI / 4.0) {
			azt = this->azimuth;
		}
		else if (this->azimuth >= M_PI / 4.0 && this->azimuth < M_PI / 2.0) {
			azt = M_PI / 2.0 - this->azimuth;
		}
		else if (this->azimuth >= M_PI / 2.0 && this->azimuth < 3.0 * M_PI / 4.0) {
			azt = this->azimuth - M_PI / 2.0;
		}
		else if (this->azimuth >= 3.0 * M_PI / 4.0 && this->azimuth < M_PI) {
			azt = M_PI - this->azimuth;
		}
		else if (this->azimuth >= M_PI && this->azimuth < 5.0 * M_PI / 4.0) {
			azt = this->azimuth - M_PI;
		}
		else if (this->azimuth >= 5.0 * M_PI / 4.0 && this->azimuth < 3.0 * M_PI / 2.0) {
			azt = 3.0 * M_PI / 2.0 - this->azimuth;
		}
		else if (this->azimuth >= 3.0 * M_PI / 2.0 && this->azimuth < 7.0 * M_PI / 4.0) {
			azt = this->azimuth - 3.0 * M_PI / 2.0;
		}
		else if (this->azimuth >= 7.0 * M_PI / 4.0) {
			azt = 2.0 * M_PI - this->azimuth;
		}

		double lateral = stand.getCellSize() / sqrt(2.0) * sin(azt + M_PI / 4.0);

		// Beam width = max lateral distance from the beam to a cell center able to own a tree which can intercept the beam.
		double R = stand.getMaximumTreeCrownRadius() + lateral;

		// Beam reach maximum distance along the beam beyond which the cells cannot own trees which can intercept the beam (too high).
		double L = stand.getMaximumTreeHeight() / (tan(this->heightAngle) + cos(this->azimuth - stand.getBottomAzimuth()) * tan(stand.getSlope())) + lateral;

		// Coordinates of the four corners of the competition rectangle.
		double sinA = sin(this->azimuth);
		double cosA = cos(this->azimuth);

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
		x_min = std::ceil(x_min / stand.getCellSize()) * stand.getCellSize();
		x_max = std::floor(x_max / stand.getCellSize()) * stand.getCellSize();
		y_min = std::ceil(y_min / stand.getCellSize()) * stand.getCellSize();
		y_max = std::floor(y_max / stand.getCellSize()) * stand.getCellSize();

		// Number of cells between min and max for both axis x and y
		int nx = std::floor((x_max - x_min) / stand.getCellSize() + 0.5) + 1;
		int ny = std::floor((y_max - y_min) / stand.getCellSize() + 0.5) + 1;

		// For each candidate relative cell of the given ray
		for (int xid = 0; xid < nx; xid++) {
			for (int yid = 0; yid < ny; yid++) {

				// Coordinates of the candidate cells (be careful : y coordinates is inverse direction of y grid)
				// id = 0 is the greater y coordinates
				double x = x_min + xid * stand.getCellSize();
				double y = y_max - yid * stand.getCellSize();

				// Add as candidate cell if it is located inside the competition rectangle
				if ((x * sinA - y * cosA < R) &&
					(x * cosA + y * sinA < L) &&
					(-x * sinA + y * cosA < R) &&
					(x * cosA + y * sinA > -R)) {

					// Relative columns and rows of the cell compared from the origin
					int col = std::floor(x / stand.getCellSize());
					int row = -std::floor(y / stand.getCellSize());

					// Add as a candidate
					this->potentialRelCells.push_back(RelativeCoords{ row, col });
				}
			}
		}
	}

	void printPotentialCells() {

		std::cout << this->id << "\t";

		int n_potcells = this->potentialRelCells.size();
		for (int i = 0; i < n_potcells; i++) {
			std::cout << "(" << this->potentialRelCells[i].row << "," << this->potentialRelCells[i].col << ")\t";
		}

	}

	ShiftedCell getPotentialCell(int rel_cell_id, Cell* target_cell, Stand& stand) {

		// Get row and column of potential cell
		int row_pot = target_cell->getRow() + this->potentialRelCells[rel_cell_id].row;
		int col_pot = target_cell->getCol() + this->potentialRelCells[rel_cell_id].col;

		// Search if cells is outside the stand
		bool is_outside = row_pot < 0 || row_pot >= stand.getNCells() || col_pot < 0 || col_pot >= stand.getNCells();

		// If cel is outside the main stand, do not keep as potential if we are not in a torus system
		if (is_outside && !stand.isUseTorus()) { return(ShiftedCell{ nullptr, 0.0, 0.0, 0.0 }); }

		// Get row and column of the original cell in the grid cell (diffrent only within a torus system, if cell is outside)
		int row_pot_original = (stand.getNCells() + (row_pot % stand.getNCells())) % stand.getNCells();
		int col_pot_original = (stand.getNCells() + (col_pot % stand.getNCells())) % stand.getNCells();

		// Get corresponding original cell
		Cell* original_cell = stand.getCell(row_pot_original, col_pot_original);

		// Check if it contains trees, otherwise, do not add as potential cell
		if ((*original_cell).isEmpty()) { return(ShiftedCell{ nullptr, 0.0, 0.0, 0.0 }); }

		// Compute shift in coordinates to apply on trees within the original cell (torus system)
		// Shift is a multiple of stand_size : how many stand size we have to shift tree coordinates times stand size, to apply torus system
		// If we have a potential cell with negative or gretaer than n_cells row or column
		// ==> Thus, we want to substract coordinates of original cells in order to have potential cell outside the main stand
		// Be careful, y - coordinates system is opposite direction of y - grid system
		double stand_size = stand.getCellSize() * stand.getNCells();

		double x_shift = (col_pot / stand.getNCells() - (col_pot % stand.getNCells() < 0 ? 1 : 0)) * stand_size; // Negative x_id ==> negative shift
		double y_shift = -(row_pot / stand.getNCells() - (row_pot % stand.getNCells() < 0 ? 1 : 0)) * stand_size; // Negative y_id ==> positive shift


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
			z_shift = -d * cos(azimuth_xy - stand.getBottomAzimuth()) * tan(stand.getSlope());
		}


		// Return the shifted cell
		return(ShiftedCell{ original_cell, x_shift, y_shift, z_shift, is_outside });
	}

	void computeInterceptions(Cell* target_cell, Stand& stand) {

		// For each possible relative cell coordinates
		int n_relcells = this->potentialRelCells.size();
		for (int i = 0; i < n_relcells; i++) {

			// Get a pointer to the original potential cell and its associated shift
			ShiftedCell pot_cell = this->getPotentialCell(i, target_cell, stand);

			// Cell is not potential if there is no trees or it is outside the main plot and torus system is disabled
			if (!pot_cell.cell) { continue; }

			// Compute interception for each tree within the potential cell (with a given shift if the potential cell is outside the plot)
			int n_pot_trees = pot_cell.cell->getNTrees();
			for (int t = 0; t < n_pot_trees; t++) {

				// Get tree pointer
				Tree* tree_p = pot_cell.cell->getTree(t);

				// Compute interception path length and distance form target cell
				// Different depending on the crown form
				if (tree_p->getCrownType() == ELLIPSOID) {
					this->computeInterceptionCrownEllipsoid(target_cell, tree_p, pot_cell);
				}
				else if (tree_p->getCrownType() == PARABOLOID) {
					this->computeInterceptionCrownParaboloid(target_cell, tree_p, pot_cell);
				}
			}
		}
	}

	void orderInterceptions() {

		std::sort(this->interceptions.begin(), this->interceptions.end(), std::greater<Interception>());

	}

	double applyBeerLambert(double incident_energy, double extinction_coef, double clumping_factor, double leaf_area_density, double path_length) {

		return(incident_energy * (1 - exp(-extinction_coef * clumping_factor * leaf_area_density * path_length)));

	}

	void summarizeInterceptions(Stand &stand, bool turbid_medium) {

		// Remove uselss memory on the interception vector
		this->interceptions.shrink_to_fit();

		// Compute projection of energy on plane parallel to slope
		double scalar_slope = cos(stand.getSlope()) * sin(this->heightAngle) +
			sin(stand.getSlope()) * cos(this->heightAngle) * cos(this->azimuth - stand.getBottomAzimuth());

		//  in MJ / m2 and convert it into MJ per cell
		double e_incident_slope_m2 = scalar_slope * this->incidentEnergy;
		double e_incident_slope_cell = e_incident_slope_m2 * stand.getCellArea();

		// Compute attenuation of energy across successives crown interceptions for each ray X target cell

		// Initialize intermediate variables
		int id_target = 0;
		double current_energy = e_incident_slope_cell;

		int n_interceptions = this->interceptions.size();
		for (int i = 0; i < n_interceptions; i++) {

			// If another target cell, reinitialized the current energy of the ray as the absolute energy coming above the canopy
			if (id_target != this->interceptions[i].getTargetCell()->getId()) {

				current_energy = e_incident_slope_cell;
				id_target = this->interceptions[i].getTargetCell()->getId();

			}

			// Compute the potential energy to the tree (energy without attenuation by neighbours)
			// and compute energy with attenuation
			// Different considering crown as turbid medium or porous envelop
			double potential_energy = 0;
			double intercepted_energy = 0;

			// Turbid medium ==> apply beer lambert law
			if (turbid_medium) {

				potential_energy = this->applyBeerLambert(
					e_incident_slope_cell,
					this->interceptions[i].getTree()->getExtinctionCoef(), 
					this->interceptions[i].getTree()->getClumpingFactor(), 
					this->interceptions[i].getTree()->getCrownLAD(),
					this->interceptions[i].getLength());

				intercepted_energy = this->applyBeerLambert(
					current_energy,
					this->interceptions[i].getTree()->getExtinctionCoef(), 
					this->interceptions[i].getTree()->getClumpingFactor(), 
					this->interceptions[i].getTree()->getCrownLAD(),
					this->interceptions[i].getLength());
			}
			// Porous envelop ==> reduce the energy by a fixed amount
			else {
				potential_energy = e_incident_slope_cell * (1 - this->interceptions[i].getTree()->getCrownOpenness());
				intercepted_energy = current_energy * (1 - this->interceptions[i].getTree()->getCrownOpenness());
			}

			// Add to the potential and intercepted energy by the tree
			this->interceptions[i].getTree()->addEnergyPotential(potential_energy);
			this->interceptions[i].getTree()->addEnergy(intercepted_energy);

			// And remove the intercepted energy from the energy that left
			current_energy -= intercepted_energy;
			
			// Remove the intercepted energy from the total energy above the cell
			this->interceptions[i].getTargetCell()->interceptEnergy(intercepted_energy);
		}

	}


private:
	void computeInterceptionCrownParaboloid(Cell* target_cell, Tree* tree, ShiftedCell& shifted_cell) {

		// GET POSITION OF THE PARABOLOID
		double x = tree->getX();
		double y = tree->getY();
		double z = tree->getZ() + tree->getCrownBaseHeight();

		// SHIFT POSITION OF THE PARABOLOID
		// 2 shifts: the first one to set the target cell as origine, and the second to acoount for outside (*tree) when torus system
		x = x + shifted_cell.xShift - target_cell->getX();
		y = y + shifted_cell.yShift - target_cell->getY();
		z = z + shifted_cell.zShift - target_cell->getZ();

		// GET PARAMETERS OF THE PARABOLOID
		double a = tree->getCrownRadius();
		double b = tree->getCrownRadius();
		double h = tree->getHeight() - tree->getCrownBaseHeight();


		// FIND SOLUTION OF THE QUADRATIC EQUATION (A * x * x + B * x + C = 0)
		// Equation giving distance between ray intersection with (*tree) crown and target cell center

		   // Compute a, b and c coefficients
		double coef_a = this->cosHeightAngle * this->cosHeightAngle * this->cosAzimuth * this->cosAzimuth / (a * a) +
			this->cosHeightAngle * this->cosHeightAngle * this->sinAzimuth * this->sinAzimuth / (b * b);

		double coef_b = -2 * x * this->cosHeightAngle * this->cosAzimuth / (a * a) -
			2 * y * this->sinAzimuth * this->cosHeightAngle / (b * b) +
			this->sinHeightAngle / h;

		double coef_c = (x * x) / (a * a) + (y * y) / (b * b) - (z + h) / h;


		// Find positive solutions of the above quadratic equations (distance to target cell)
		// if not 2 solutions, set to NaN
		// Because negative = no interception and null = 1 interception = do not consider tangent rays with crown
		double delta = coef_b * coef_b - 4 * coef_a * coef_c;

		double sol1 = std::numeric_limits<double>::quiet_NaN();
		double sol2 = std::numeric_limits<double>::quiet_NaN();

		if (delta > 0) {
			double delta_sqrt = sqrt(delta);
			sol1 = (-coef_b + delta_sqrt) / (2 * coef_a);
			sol2 = (-coef_b - delta_sqrt) / (2 * coef_a);
		}



		// COMPUTE COORDINATES OF INTERCEPTION POINTS WITH CROWN LIMITS AND TARGET POINT

			// Interception coords of the two roots
		double v1_x = sol1 * this->cosHeightAngle * this->cosAzimuth;
		double v1_y = sol1 * this->cosHeightAngle * this->sinAzimuth;
		double v1_z = sol1 * this->sinHeightAngle;

		double v2_x = sol2 * this->cosHeightAngle * this->cosAzimuth;
		double v2_y = sol2 * this->cosHeightAngle * this->sinAzimuth;
		double v2_z = sol2 * this->sinHeightAngle;

		// Interception with the base plane
		double solb = z / this->sinHeightAngle;

		double vb_x = solb * this->cosHeightAngle * this->cosAzimuth;
		double vb_y = solb * this->cosHeightAngle * this->sinAzimuth;
		double vb_z = z; // solb * this->sinHeightAngle



		// FIND SOLUTIONS

			// Get the limits of the box between shifted paraboloid center and crown limit
		double x_bbox_max = x + a;
		double y_bbox_max = y + b;
		double z_bbox_max = z + h;

		double x_bbox_min = x - a;
		double y_bbox_min = y - b;
		double z_bbox_min = z;

		// Find interception points within the crown (edge of ellipsoid or z axis if half ellipsoid = crown base plane)

		bool is_v1 =
			// interception 2 with crown is above the ground AND
			v1_z >= 0 &&
			v1_x >= x_bbox_min && v1_x <= x_bbox_max &&
			v1_y >= y_bbox_min && v1_y <= y_bbox_max &&
			// interception 2 with crown is in crown bbox AND
			v1_z >= z_bbox_min && v1_z <= z_bbox_max;


		bool is_v2 =
			// interception 2 with crown is above the ground AND
			v2_z >= 0 &&
			// interception 2 with crown is in crown bbox AND
			v2_x >= x_bbox_min && v2_x <= x_bbox_max &&
			v2_y >= y_bbox_min && v2_y <= y_bbox_max &&
			v2_z >= z_bbox_min && v2_z <= z_bbox_max;


		bool is_vb =
			// interception with crown base is above the ground AND
			vb_z >= 0 &&
			// interception with crown base is in crown bbox AND
			vb_x >= x_bbox_min && vb_x <= x_bbox_max &&
			vb_y >= y_bbox_min && vb_y <= y_bbox_max &&
			vb_z >= z_bbox_min && vb_z <= z_bbox_max &&
			// interception with crown base is in the paraboloid
			((vb_x - x) * (vb_x - x) / (a * a) +
				(vb_y - y) * (vb_y - y) / (b * b) -
				(vb_z - z) / (h * h)) <= 1;


		bool is_target =
			// interception with target point is in crown bbox AND
			0 >= x_bbox_min && 0 <= x_bbox_max &&
			0 >= y_bbox_min && 0 <= y_bbox_max &&
			0 >= z_bbox_min && 0 <= z_bbox_max &&
			// interception with target point is in the paraboloid
			(x * x / (a * a) +
				y * y / (b * b) -
				-z / (h * h)) <= 1;


		// FIND PATH LENGTH AND DISTANCE TO CELL CENTER
		// Distance between target cell and middle point between two interceptions
		// Length is the length of the ray path across the crown
		double distance = 0.0;
		double length = 0.0;

		// If 0 solution, no interception
		if (!(!std::isnan(sol1) && is_v1) && !(!std::isnan(sol2) && is_v2) && !is_vb && !is_target) {
			return;
		}

		std::stack<double> distance_tmp;
		std::stack<double> length_tmp;

		int nsols = 0;
		if (!std::isnan(sol1) && is_v1) {
			nsols++;
			distance_tmp.push(sol1);
			length_tmp.push(sol1);
		}

		if (!std::isnan(sol2) && is_v2) {
			nsols++;
			distance_tmp.push(sol2);
			length_tmp.push(sol2);
		}

		if (is_vb) {
			nsols++;
			distance_tmp.push(solb);
			length_tmp.push(solb);
		}

		if (is_target) {
			nsols++;
			distance_tmp.push(0.0);
			length_tmp.push(0.0);
		}


		// Compute length and distance if 2 solutions
		if (nsols == 2) {

			// Length
			double l1 = length_tmp.top();
			length_tmp.pop();
			double l2 = length_tmp.top();

			length = abs(l1 - l2);

			// Distance
			double d1 = distance_tmp.top();
			distance_tmp.pop();
			double d2 = distance_tmp.top();

			distance = (d1 + d2) / 2.0;
		}
		// Otherwise, there is a problem
		else {
			std::cout << "Not 0 or 2 solutions" << std::endl;
		}

		// Create and push interception to the vector of ray interceptions
		#ifdef _OPENMP
			#pragma omp critical
		#endif
		{
			this->interceptions.emplace_back(Interception(target_cell, tree, length, distance, shifted_cell.outside));
		}
	}

	void computeInterceptionCrownEllipsoid(Cell* target_cell, Tree* tree, ShiftedCell& shifted_cell) {

		// GET POSITION OF THE PARABOLOID
		double x = tree->getX();
		double y = tree->getY();
		double z = tree->getZ() + tree->getCrownBaseHeight() + (tree->getHeight() - tree->getCrownBaseHeight()) / 2.0;

		// SHIFT POSITION OF THE PARABOLOID
		// 2 shifts: the first one to set the target cell as origine, and the second to acoount for outside (*tree) when torus system
		x = x + shifted_cell.xShift - target_cell->getX();
		y = y + shifted_cell.yShift - target_cell->getY();
		z = z + shifted_cell.zShift - target_cell->getZ();

		// GET PARAMETERS OF THE PARABOLOID
		double a = tree->getCrownRadius();
		double b = tree->getCrownRadius();
		double c = (tree->getHeight() - tree->getCrownBaseHeight()) / 2.0;


		// FIND SOLUTION OF THE QUADRATIC EQUATION (A * x * x + B * x + C = 0)
		// Equation giving distance between ray intersection with (*tree) crown and target cell center

		   // Compute a, b and c coefficients
		double coef_a = this->cosHeightAngle * this->cosHeightAngle * this->cosAzimuth * this->cosAzimuth / (a * a) +
			this->cosHeightAngle * this->cosHeightAngle * this->sinAzimuth * this->sinAzimuth / (b * b) +
			this->sinHeightAngle * this->sinHeightAngle / (c * c);

		double coef_b = -2 * x * this->cosHeightAngle * this->cosAzimuth / (a * a) -
			2 * y * this->sinAzimuth * this->cosHeightAngle / (b * b) -
			2 * z * this->sinHeightAngle / (c * c);

		double coef_c = (x * x) / (a * a) + (y * y) / (b * b) + (z * z) / (c * c) - 1;


		// Find positive solutions of the above quadratic equations (distance to target cell)
		// if not 2 solutions, set to NaN
		// Because negative = no interception and null = 1 interception = do not consider tangent rays with crown
		double delta = coef_b * coef_b - 4 * coef_a * coef_c;

		double sol1 = std::numeric_limits<double>::quiet_NaN();
		double sol2 = std::numeric_limits<double>::quiet_NaN();

		if (delta > 0) {
			double delta_sqrt = sqrt(delta);
			sol1 = (-coef_b + delta_sqrt) / (2 * coef_a);
			sol2 = (-coef_b - delta_sqrt) / (2 * coef_a);
		}



		// COMPUTE COORDINATES OF INTERCEPTION POINTS WITH CROWN LIMITS AND TARGET POINT

			// Interception coords of the two roots
		double v1_x = sol1 * this->cosHeightAngle * this->cosAzimuth;
		double v1_y = sol1 * this->cosHeightAngle * this->sinAzimuth;
		double v1_z = sol1 * this->sinHeightAngle;

		double v2_x = sol2 * this->cosHeightAngle * this->cosAzimuth;
		double v2_y = sol2 * this->cosHeightAngle * this->sinAzimuth;
		double v2_z = sol2 * this->sinHeightAngle;

		// FIND SOLUTIONS

			// Get the limits of the box between shifted paraboloid center and crown limit
		double x_bbox_max = x + a;
		double y_bbox_max = y + b;
		double z_bbox_max = z + c;

		double x_bbox_min = x - a;
		double y_bbox_min = y - b;
		double z_bbox_min = z - c;

		// Find interception points within the crown (edge of ellipsoid or z axis if half ellipsoid = crown base plane)

		bool is_v1 =
			// interception 2 with crown is above the ground AND
			v1_z >= 0 &&
			v1_x >= x_bbox_min && v1_x <= x_bbox_max &&
			v1_y >= y_bbox_min && v1_y <= y_bbox_max &&
			// interception 2 with crown is in crown bbox AND
			v1_z >= z_bbox_min && v1_z <= z_bbox_max;


		bool is_v2 =
			// interception 2 with crown is above the ground AND
			v2_z >= 0 &&
			// interception 2 with crown is in crown bbox AND
			v2_x >= x_bbox_min && v2_x <= x_bbox_max &&
			v2_y >= y_bbox_min && v2_y <= y_bbox_max &&
			v2_z >= z_bbox_min && v2_z <= z_bbox_max;


		// FIND PATH LENGTH AND DISTANCE TO CELL CENTER
		// Distance between target cell and middle point between two interceptions
		// Length is the length of the ray path across the crown
		double distance = 0.0;
		double length = 0.0;

		// If 0 solution, no interception
		if (!(!std::isnan(sol1) && is_v1) && !(!std::isnan(sol2) && is_v2)) {
			return;
		}

		std::stack<double> distance_tmp;
		std::stack<double> length_tmp;

		int nsols = 0;
		if (!std::isnan(sol1) && is_v1) {
			nsols++;
			distance_tmp.push(sol1);
			length_tmp.push(sol1);
		}

		if (!std::isnan(sol2) && is_v2) {
			nsols++;
			distance_tmp.push(sol2);
			length_tmp.push(sol2);
		}

		// Compute length and distance if 2 solutions
		if (nsols == 2) {

			// Length
			double l1 = length_tmp.top();
			length_tmp.pop();
			double l2 = length_tmp.top();

			length = abs(l1 - l2);

			// Distance
			double d1 = distance_tmp.top();
			distance_tmp.pop();
			double d2 = distance_tmp.top();

			distance = (d1 + d2) / 2.0;
		}
		// Otherwise, there is a problem
		else {
			std::cout << "Not 0 or 2 solutions" << std::endl;
		}

		// Create and push interception to the vector of ray interceptions
		#ifdef _OPENMP
			#pragma omp critical
		#endif
		{
			this->interceptions.emplace_back(Interception(target_cell, tree, length, distance, shifted_cell.outside));
		}
	}

};

class Model {

private:
	// Stand object woth cells, trees and geometric variables
	Stand stand;

	// Vector of rays
	std::vector<Ray*> rays;

	// Light characteristics
	bool turbidMedium; // Choose light interception method by the crowns: either porous envelop (use "crownOpenness") or turbid medium (use "crownLAD")
	double totalEnergy; // Total energy coming from diffuse and direct rays above canopy (in MJ/m2)


public:

	Model(DataFrame trees, DataFrame cells, DataFrame rays,
		double e_above_m2,
		double slope, double north_to_x_cw, double aspect,
		double cell_size, double n_cells,
		bool use_torus, bool turbid_medium) :

		stand(cells, trees, e_above_m2, slope, north_to_x_cw, aspect, cell_size, n_cells, use_torus)
	{
		this->initRays(rays);
		this->turbidMedium = turbid_medium;
	}

	~Model() {
		// Delete rays pointers
		int n_rays = this->rays.size();
		for (int i = 0; i < n_rays; i++) {
			delete this->rays[i];
		}
	}

	// Getters
	double getMaximumTreeHeight() { return(this->stand.getMaximumTreeHeight()); }
	double getMaximumTreeCrownRadius() { return(this->stand.getMaximumTreeCrownRadius()); }
	int getNPotCells() {
		int n_potcells = 0;
		int n_rays = this->rays.size();
		for (int i = 0; i < n_rays; i++) {
			n_potcells += this->rays[i]->getNPotCells();
		}
		return(n_potcells);
	}

	// Methods
	void findPotentialCells() {

		int n_rays = this->rays.size();
		for (int i = 0; i < n_rays; i++) {

			this->rays[i]->findPotentialCells(this->stand);
		}
	}

	void printPotentialCells() {

		int n_rays = this->rays.size();
		for (int i = 0; i < n_rays; i++) {

			this->rays[i]->printPotentialCells();
		}

	}

	void computeInterceptions() {

		int n_rays = this->rays.size();

		// Get target cells
		std::map<int, Cell*> cells = this->stand.getCells();
		int n_cells = cells.size();

		// For each target cell
		#ifdef _OPENMP
			#pragma omp parallel for
		#endif
		for (int c = 0; c < n_cells; c++) {

			// Get interceptions for all rays X all trees
			for (int r = 0; r < n_rays; r++) {

				this->rays[r]->computeInterceptions(cells[c], this->stand);
			}
		}
	}

	void orderInterceptions() {

		int n_rays = this->rays.size();
		for (int r = 0; r < n_rays; r++) {

			this->rays[r]->orderInterceptions();
		}

	}

	void summarizeInterceptions() {

		int n_rays = this->rays.size();
		for (int i = 0; i < n_rays; i++) {

			this->rays[i]->summarizeInterceptions(this->stand, this->turbidMedium);
		}
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
				lci_tree[itree] = tree->getLci();

				itree++;
			}

			icell++;
		}

		// Create trees and cells RCPP DataFrames
		DataFrame output_trees = DataFrame::create(
			Named("id_tree") = id_tree,
			Named("epot") = epot_tree,
			Named("e") = e_tree,
			Named("lci") = lci_tree
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


private:
	void initRays(DataFrame rays) {
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
	// Initialize the model
	Model sl_model = Model(trees, cells,
		rays, total_energy_m2,
		slope, north_to_x_cw, aspect,
		cell_size, n_cells, use_torus, turbid_medium);

	// [OPTIMIZATION]: find the extend of possible interception of each ray
	// i.e. for each ray coming to an unknown target cell, find relative cells around the target containing tree that could possibly intercept the ray
	// We can whether or not consider a torus system
	sl_model.findPotentialCells();

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





