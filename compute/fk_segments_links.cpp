/**
 * @file fk_segments_links.cpp
 * @brief Gukov-Manolescu invariant computation for links
 * 
 * This file implements computation of the Gukov-Manolescu two-variable series F_K(x,q)
 * for knot complements using braid representations. The Gukov-Manolescu invariant
 * is obtained through parametric resurgence from asymptotic expansions of colored
 * Jones polynomials and satisfies recurrence relations given by quantum A-polynomials.
 * 
 * The computation processes crossing data, applies quantum binomial coefficients
 * and q-Pochhammer symbols, and outputs polynomial coefficients with their q-powers
 * in JSON format.
 * 
 * Reference: S. Gukov, C. Manolescu, "A two-variable series for knot complements"
 * 
 * TODO: Implement multithreading at top level of recursion
 */

#include <vector>
#include <cstddef>
#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <cmath>
#include <array>

#include "solution_pool_1a_double_links.hpp"
#include "string_to_int.hpp"
#include "linalg.hpp"
#include "bilvector.hpp"
#include "qalg_links.hpp"

/**
 * @class FK
 * @brief Main class for computing Gukov-Manolescu invariants of links
 * 
 * This class reads braid and crossing data from a CSV file, processes it through
 * quantum algebraic operations (q-binomial coefficients, q-Pochhammer symbols),
 * and outputs the resulting two-variable series F_K(x,q) coefficients in JSON format.
 * 
 * The computation follows the Gukov-Manolescu framework for knot complement invariants,
 * which connect to 3-manifold invariants through Dehn surgery relationships.
 */
class FK {
    private:
        // Core braid data
        std::vector<int> braid;                     ///< Braid word representation
        std::vector<int> signs;                     ///< Crossing signs
        std::vector<int> segments;                  ///< Braid segments
        
        // Matrix computation blocks - used for representing matrices as vectors
        // Each block stores products of preceding dimension lengths for vector indexing
        std::vector<int> acc_blocks = {1};          ///< Block structure: for each matrix dimension, product of all preceding dimension lengths
        
        // Link component data (0-indexed)  
        std::vector<int> closed_strand_components = {}; ///< Components of closed strands
        
        // Crossing geometry data
        // When viewing braid left-to-right: top/bottom components
        // When viewing braid bottom-to-top: left/right components  
        std::vector<int> top_crossing_components = {};   ///< Top components at crossings (left when viewed bottom-to-top)
        std::vector<int> bottom_crossing_components = {}; ///< Bottom components at crossings (right when viewed bottom-to-top)
        
        // Crossing matrix data - stores coordinates of segments adjacent to crossings
        // First coordinate component: strand-wise, Second: crossing-wise
        // Order for each crossing: i, j, i', j' (following Sunghyuk Kim's Inverted State Sums notation)
        std::vector<std::vector<std::array<int, 2>>> matrices = {}; ///< Crossing matrix coordinates [strand][crossing]
        
        // Crossing types - four combinations of crossing sign and inversion status (not Reidemeister moves)
        std::vector<int> r = {};                    ///< Crossing types: 1,2,3,4 for combinations of sign/inversion
        
        std::vector<int> nontrivial_map;           ///< Mapping for non-trivial elements
        
        // Main computation data
        std::vector<bilvector<int>> acc;           ///< Accumulator for the Gukov-Manolescu invariant values
        std::vector<std::vector<std::vector<int>>> assignment;    ///< Symbolic (vector) variable assignments for each braid segment  
        std::vector<std::vector<int>> numerical_assignment;       ///< Numerical variable assignments for each braid segment
        
        // Topological invariants
        int components;     ///< Number of link components (0-indexed)
        int writhe = 0;     ///< Writhe of the link
        int prefactors;     ///< Number of prefactors (equals number of closed strands)
        int crossings;      ///< Number of crossings
        int degree;         ///< Degree bound for the computation
        /**
         * @brief Computes polynomial contribution for a specific angle state assignment
         * @param angles Vector of angle parameters from ILP solver
         * 
         * This function is called for each feasible angle state assignment found by the ILP solver.
         * It creates a polynomial term contributing to the Gukov-Manolescu invariant, starting at:
         * - x_acc degrees in x for each component variable (x, y, z, ...)  
         * - q_acc degree in q
         * 
         * The computation applies quantum binomial coefficients and q-Pochhammer symbols
         * according to the R-matrix prescription from Sunghyuk Kim's framework.
         * 
         * Finally, the resulting polynomial term is added to the FK accumulator `acc`.
         */
        void f (const std::vector<int>& angles) {
            // Convert symbolic assignments to numerical values using dot product with angles
            for (int i = 0; i < crossings + 1; i++) {
                for (int j = 0; j < prefactors + 1; j++) {
                    numerical_assignment[i][j] = dot(assignment[i][j], angles);
                }
            }
            
            // Initialize q-power accumulator starting from writhe and prefactor contributions
            double q_acc_double = (writhe - prefactors) / 2.0;
            
            // Initialize x-power accumulators for each component (x, y, z, ...)
            std::vector<double> x_acc_double(components, 0);
            
            // Initialize sign factor
            int init = -1;
            
            // Apply strand closing contributions
            for (int i = 0; i < prefactors; i++) {
                q_acc_double -= numerical_assignment[0][i + 1];  // q-power contribution from closing strands
                x_acc_double[closed_strand_components[i]] -= 0.5; // x-power contribution from closing strands
            }
            
            // Normalize by FK of trivial link - subtract 1/2 from each component
            for (int comp = 0; comp < components; comp++) {
                x_acc_double[comp] -= 0.5;
            }

            // Process each crossing according to its type
            for (int ind = 0; ind < crossings; ind++) {
                // Extract segment values at crossing (i, j, i', j' notation from Kim's paper)
                // In Kim's notation: i, j, i'=k, j'=m
                int i = numerical_assignment[matrices[ind][0][0]][matrices[ind][0][1]];  // i
                int k = numerical_assignment[matrices[ind][2][0]][matrices[ind][2][1]];  // i' (Kim's notation)
                int j = numerical_assignment[matrices[ind][1][0]][matrices[ind][1][1]];  // j  
                int m = numerical_assignment[matrices[ind][3][0]][matrices[ind][3][1]];  // j' (Kim's notation)
                
                // Get component indices for this crossing
                int t = top_crossing_components[ind];    // Top component (left when viewed bottom-to-top)
                int b = bottom_crossing_components[ind]; // Bottom component (right when viewed bottom-to-top)
                
                // Apply R-matrix contributions based on crossing type (sign/inversion combination)
                if (r[ind] == 1 || r[ind] == 2) {  // Crossing types 1,2                               
                    q_acc_double += (j + m) / 2.0 + j * m;           // q-power contribution
                    x_acc_double[t] += ((j + k + 1) / 4.0);          // x-power for top component
                    x_acc_double[b] += ((3 * m - i + 1) / 4.0);      // x-power for bottom component
                }
                else if (r[ind] == 4) {  // Crossing type 4
                    q_acc_double -= (i + k + m * (m + 1) - i * (i + 1)) / 2.0 + i * k;  // q-power contribution
                    x_acc_double[t] -= ((3 * j - k + 1) / 4.0);      // x-power for top component
                    x_acc_double[b] -= ((i + m + 1) / 4.0);          // x-power for bottom component
                    if ((j - k) % 2 == 0) {
                        init *= -1;  // Sign change due to Pochhammer convention (positive x powers only)
                    }
                }
                else if (r[ind] == 3) {  // Crossing type 3
                    q_acc_double -= (i + k + m * (m + 1) - i * (i + 1)) / 2.0 + i * k;  // q-power contribution
                    x_acc_double[t] -= ((3 * j - k + 1) / 4.0);      // x-power for top component
                    x_acc_double[b] -= ((i + m + 1) / 4.0);          // x-power for bottom component
                    if ((k - j) % 2 == 1) {
                        init *= -1;  // Sign change due to Pochhammer convention (positive x powers only)
                    }
                }
            }

            // Debug output for accumulated powers
            std::cout << "q_acc_double : " <<  q_acc_double << "\n";
            std::cout << "x_acc_double : " <<  x_acc_double[0] << " " << x_acc_double[1] << "\n";

            // Convert double powers to integer powers (floor operation)
            int q_acc = static_cast<int>(std::floor(q_acc_double));
            
            // Set up integer x-powers and compute matrix block structure
            std::vector<int> x_acc(components);
            std::vector<int> X_MAX(components);        // Max indices for each component's matrix (degree - offset)
            std::vector<int> blocks(components);       // Block sizes for matrix-as-vector indexing
            blocks[0] = 1;
            
            for (int n = 0; n < components; n++) {
                x_acc[n] = x_acc_double[n];            // Convert to integer x-power (offset)
                X_MAX[n] = degree - x_acc[n];          // Max index = requested degree - offset (avoids computing past desired degree)
                if (n != 0) {
                    // Block size = product of all preceding dimension lengths
                    blocks[n] = (X_MAX[n - 1] + 1) * blocks[n - 1];
                }
            }
            
            // Total size of polynomial term storage
            int prod = blocks[components - 1] * (X_MAX[components - 1] + 1);

            // Initialize polynomial term storage with initial coefficient
            std::vector<bilvector<int>> term(prod,  bilvector<int>(0, 1, 20, 0));
            term[0][0] = init;  // Set initial coefficient with computed sign
            
            // Apply quantum binomial coefficients according to Kim's R-matrix prescription
            for (int ind = 0; ind < crossings; ind++) {
                if (r[ind] == 1) {  // Crossing type 1
                    int i = numerical_assignment[matrices[ind][0][0]][matrices[ind][0][1]];  // i
                    int m = numerical_assignment[matrices[ind][3][0]][matrices[ind][3][1]];  // j' (m)
                    if (i > 0) {
                        pp_q_binom(term, i, i - m, false);  // Positive-positive q-binomial
                    }
                    else {
                        np_q_binom(term, i, i - m, false);  // Negative-positive q-binomial  
                    }
                }
                else if (r[ind] == 2) {  // Crossing type 2
                    int i = numerical_assignment[matrices[ind][0][0]][matrices[ind][0][1]];  // i
                    int m = numerical_assignment[matrices[ind][3][0]][matrices[ind][3][1]];  // j' (m)
                    np_q_binom(term, i, m, false);  // Negative-positive q-binomial
                }
                else if (r[ind] == 3) {  // Crossing type 3    
                    int j = numerical_assignment[matrices[ind][1][0]][matrices[ind][1][1]];  // j
                    int k = numerical_assignment[matrices[ind][2][0]][matrices[ind][2][1]];  // i' (k)
                    np_q_binom(term, j, k, true);  // Negative-positive q-binomial (inverted)
                }
                else {  // Crossing type 4
                    int j = numerical_assignment[matrices[ind][1][0]][matrices[ind][1][1]];  // j
                    int k = numerical_assignment[matrices[ind][2][0]][matrices[ind][2][1]];  // i' (k)
                    if (j > 0) {
                        pp_q_binom(term, j, j - k, true);  // Positive-positive q-binomial (inverted)
                    }
                    else {
                        np_q_binom(term, j, j - k, true);  // Negative-positive q-binomial (inverted)
                    }
                }
            }
            
            // Apply q-Pochhammer symbols according to Kim's R-matrix prescription  
            for (int ind = 0; ind < crossings; ind++) {
                if (r[ind] == 1) {  // Crossing type 1
                    int b = bottom_crossing_components[ind];  // Bottom component
                    int j = numerical_assignment[matrices[ind][1][0]][matrices[ind][1][1]];  // j
                    int k = numerical_assignment[matrices[ind][2][0]][matrices[ind][2][1]];  // i' (k)
                    x_q_pochhammer(term, k, j + 1, b, components, X_MAX, blocks);  // (x_b*q)_k / (x_b*q)_{j+1}
                }
                else if (r[ind] == 2) {  // Crossing type 2
                    int b = bottom_crossing_components[ind];  // Bottom component
                    int j = numerical_assignment[matrices[ind][1][0]][matrices[ind][1][1]];  // j
                    int k = numerical_assignment[matrices[ind][2][0]][matrices[ind][2][1]];  // i' (k)
                    x_q_inv_pochhammer(term, j, k + 1, b, components, X_MAX, blocks);  // Inverse Pochhammer
                }
                else if (r[ind] == 3) {  // Crossing type 3
                    int t = top_crossing_components[ind];  // Top component
                    int i = numerical_assignment[matrices[ind][0][0]][matrices[ind][0][1]];  // i
                    int m = numerical_assignment[matrices[ind][3][0]][matrices[ind][3][1]];  // j' (m)
                    x_q_inv_pochhammer(term, i, m + 1, t, components, X_MAX, blocks);  // Inverse Pochhammer
                }
                else {  // Crossing type 4
                    int t = top_crossing_components[ind];  // Top component
                    int i = numerical_assignment[matrices[ind][0][0]][matrices[ind][0][1]];  // i
                    int m = numerical_assignment[matrices[ind][3][0]][matrices[ind][3][1]];  // j' (m)
                    x_q_pochhammer(term, m, i + 1, t, components, X_MAX, blocks);  // (x_t*q)_m / (x_t*q)_{i+1}
                }
            } 
            
            // Add computed polynomial term to main accumulator
            offset_addition(acc, term, x_acc, q_acc, components, X_MAX, 1, acc_blocks, blocks);
        }
        /**
         * @brief Writes the computed Gukov-Manolescu invariant to JSON file
         * @param file_ Output filename (without extension)
         * 
         * Outputs the accumulated polynomial coefficients in JSON format as
         * "coefficient_q_powers" array, where each entry contains [q_power, coefficient] pairs
         * for the two-variable series F_K(x,q).
         */
        void write (std::string file_) {
            // Open JSON output file
            std::ofstream file;
            file.open(file_ + ".json");
            
            // Write JSON header for coefficient-q_power pairs
            file << "{\n\t\"coefficient_q_powers\":[\n";
            
            // Iterate through each component combination in the accumulator
            for (int i = 0; i < acc.size(); i++) {
                file << "\t\t[";  // Start array for this component combination
                
                bool first_write = true;
                
                // Write all non-zero coefficient-q_power pairs for this component combination
                for (int j = acc[i].get_max_nindex(); j <= acc[i].get_max_pindex(); j++) {
                    if (acc[i][j] != 0) {  // Only write non-zero coefficients
                        if (!first_write) {
                            file << ",[";  // Comma separator for subsequent pairs
                        }
                        else {
                            file << "[";   // First pair needs no comma
                            first_write = false;
                        }
                        file << j << "," << acc[i][j] << "]";  // [q_power, coefficient]
                    }
                }
                
                // End array for this component combination
                if (i < acc.size() - 1) {
                    file << "],\n";  // Comma for more combinations
                }
                else {
                    file << "]\n";   // No comma for last combination
                }
            }
            
            // Close JSON structure
            file << "\t]\n}";
            file.close();
        }
    public:
        std::vector<std::vector<double>> inequalities;  ///< Linear inequalities read from ILP data file
        std::vector<std::vector<double>> criteria;      ///< Criteria constraints read from ILP data file  
        std::string metadata;                           ///< Metadata read from ILP file (degree, braid, r matrix types, number of crossings)
        /**
         * @brief Constructor that reads input data and computes Gukov-Manolescu invariant
         * @param infile_ Input CSV filename (without extension)  
         * @param outfile_ Output JSON filename (without extension)
         * 
         * Reads braid and crossing data from CSV, sets up the computation framework,
         * calls the ILP solver to find feasible angle assignments, processes each
         * assignment through function f, and writes the result to JSON.
         */
        FK (std::string infile_, std::string outfile_) {

            /**
             * we're going to read in the following data from the infile csv file:
             * - degree
             * - number of components
             * - writhe
             * - 2 * n_crossings numbers, the absolute values of the crossing generators and the R-matrix types (gen1, r1, gen2, r2, ...)
             * - the components of the closed strands (0-indexed)
             * - the top and bottom components at each crossing (0-indexed)
             * - criteria for the ILP solver (one per line, comma-separated)
             * - /
             * - inequalities for the ILP solver (one per line, comma-separated)
             * - /
             * - (vector-)symbolic variable assignments for each braid segment (one line per segment, comma-separated, first prefactors + 1 for first segment, then crossings lines of prefactors + 1 each)
            */

            std::ifstream infile;
            infile.open(infile_ + ".csv");
            if (infile.is_open())
            {
                std::string line;

                // Read degree bound for computation
                std::getline (infile,line,'\n');
                int index = line.find(",");
                degree = string_to_int(line.substr(0, index));

                // Read number of link components
                std::getline (infile,line,'\n');
                index = line.find(",");
                components = string_to_int(line.substr(0, index));

                // Read writhe of the link
                std::getline (infile,line,'\n');
                index = line.find(",");
                writhe = string_to_int(line.substr(0, index));

                // Initialize main accumulator with appropriate size
                acc.resize(std::pow(degree + 1, components), bilvector<int>(0, 1, 20, 0));

                // Read crossing data: generators and R-matrix types
                std::getline (infile,line,'\n');
                int height = 0;
                index = line.find(",");
                bool done = false;
                while (true) {
                    if (index == -1) {
                        break;
                    }
                    
                    // Read crossing generator (absolute value)
                    int c = string_to_int(line.substr(0, index));
                    
                    // Store matrix coordinates for segments adjacent to this crossing
                    // Following Kim's notation: i, j, i', j' at positions [height, c-1], [height, c], [height+1, c-1], [height+1, c]
                    matrices.push_back({
                        {height, c - 1},      // i position
                        {height, c},          // j position  
                        {height + 1, c - 1},  // i' position
                        {height + 1, c}       // j' position
                    });
                    line = line.substr(index + 1, line.size() - index - 1);
                    index = line.find(",");

                    // Read R-matrix type (1,2,3,4 for sign/inversion combinations)
                    r.push_back(string_to_int(line.substr(0, index)));
                    line = line.substr(index + 1, line.size() - index - 1);

                    height++;  // Move to next crossing level
                    index = line.find(",");
                }
                // Read closed strand component assignments
                std::getline (infile,line,'\n');
                index = line.find(",");
                done = false;
                while (true) {
                    if (index == -1) {
                        break;
                    }
                    closed_strand_components.push_back(string_to_int(line.substr(0, index)));
                    line = line.substr(index + 1, line.size() - index - 1);
                    index = line.find(",");
                }
                
                // Set derived quantities
                prefactors = closed_strand_components.size();  // Number of prefactors = number of closed strands
                crossings = r.size();                          // Number of crossings
                
                // Initialize assignment matrices
                assignment.resize(crossings + 1, std::vector<std::vector<int>>(prefactors + 1));
                numerical_assignment.resize(crossings + 1, std::vector<int>(prefactors + 1));
                // Read top and bottom component assignments for each crossing
                std::getline (infile,line,'\n');
                index = line.find(",");
                done = false;
                while (true) {
                    if (index == -1) {
                        break;
                    }
                    // Read top component (left component when viewed bottom-to-top)
                    top_crossing_components.push_back(string_to_int(line.substr(0, index)));
                    line = line.substr(index + 1, line.size() - index - 1);
                    index = line.find(",");
                    
                    // Read bottom component (right component when viewed bottom-to-top)  
                    bottom_crossing_components.push_back(string_to_int(line.substr(0, index)));
                    line = line.substr(index + 1, line.size() - index - 1);
                    index = line.find(",");
                }
                // Parse remaining data in three stages separated by "/" lines
                int stage = 0;
                int criteria_index = 0;
                int inequality_index = 0;
                int extension_index = 0;
                
                while ( std::getline (infile,line,'\n') )
                {
                    if (line[0] == '/') {
                        stage++;  // Move to next parsing stage
                    }
                    else if (stage == 0) {  // Stage 0: Read criteria for ILP solver
                        criteria.push_back({});
                        done = false;
                        index = line.find(",");
                        while (true) {
                            if (index == -1) {
                                break;
                            }
                            criteria[criteria_index].push_back(string_to_double(line.substr(0, index)));
                            line = line.substr(index + 1, line.size() - index - 1);
                            index = line.find(",");
                        }
                        criteria_index++;
                    }
                    else if (stage == 1) {  // Stage 1: Read linear inequalities for ILP solver
                        inequalities.push_back({});
                        done = false;
                        index = line.find(",");
                        while (true) {
                            if (index == -1) {
                                break;
                            }
                            inequalities[inequality_index].push_back(string_to_int(line.substr(0, index)));
                            line = line.substr(index + 1, line.size() - index - 1);
                            index = line.find(",");
                        }
                        inequality_index++;
                    }
                    else if (stage == 2) {  // Stage 2: Read symbolic variable assignments for each braid segment
                        done = false;
                        index = line.find(",");
                        while (true) {
                            if (index == -1) {
                                break;
                            }
                            // Store assignment coefficients in row-major order (that is, the assignment coefficient rows are stored by keeping the crossing level the same and moving across all strands, then moving to the next crossing level, et cetera)
                            assignment[extension_index % (crossings + 1)][extension_index / (crossings + 1)].push_back(string_to_int(line.substr(0, index)));
                            line = line.substr(index + 1, line.size() - index - 1);
                            index = line.find(",");
                        }
                        extension_index++;
                    }
                }
            }
            else {
                std::cout << "ERROR: Unable to open file '" + infile_ + ".csv'!";
                exit(0);
            }
            
            // Set up accumulator block structure for matrix-as-vector operations
            for (int i = 1; i < components; i++) {
                acc_blocks.push_back(acc_blocks[i - 1] * (degree + 1));
            }
            
            // Create function wrapper for ILP solver callback
            std::function<void(const std::vector<int>&)> f_wrapper = [this](const std::vector<int>& v) { f(v); };
            
            // Run ILP solver to find feasible angle assignments and process each one
            pooling(criteria, inequalities, f_wrapper);
            
            // Apply final offset to accumulator (normalization step)
            std::vector<int> increment_offset(components);
            increment_offset[0] = 1;
            std::vector<int> maxima(components, degree - 1);
            offset_addition(acc, acc, increment_offset, 0, components, maxima, -1, acc_blocks, acc_blocks);
            
            // Write results to output file
            write(outfile_);
        }
};

/**
 * @brief Main entry point for Gukov-Manolescu invariant computation
 * @param argc Number of command line arguments (should be 3)
 * @param argv Command line arguments: [program_name] [input_ILP_file] [output_FK_file]
 * @return 0 on success, 1 on error
 * 
 * Usage: ./program input_ILP_filename output_FK_filename
 * - input_ILP_filename: CSV file containing braid and ILP data (without .csv extension)  
 * - output_FK_filename: JSON file to write FK results (without .json extension)
 */
int main(int argc, char* argv[]) {
    if (argc == 3) {
        FK(argv[1], argv[2]);
        return 0;
    }
    else {
        std::cout << "ERROR: Incorrect number of arguments! Please provide input and output file names without extensions." << std::endl;
        return 1;
    }
}
