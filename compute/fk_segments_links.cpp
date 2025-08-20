// implement multithreading at top level of recursion

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

class FK {
    private:
        std::vector<int> angles;
        std::vector<int> angle_signs;
        std::vector<int> braid;
        std::vector<int> signs;
        std::vector<int> segments;
        std::vector<int> acc_blocks = {1};
        std::vector<int> closed_strand_components = {};
        std::vector<std::vector<std::array<int, 2>>> matrices = {};
        std::vector<int> r = {};
        std::vector<std::vector<int>> transitions;
        std::vector<bool> trivial_angles_;
        std::vector<int> nontrivial_map;
        std::vector<int> inversion_data;
        std::vector<bilvector<int>> acc;
        std::vector<std::vector<std::vector<int>>> assignment;
        std::vector<std::vector<int>> numerical_assignment;
        std::vector<int> top_crossing_components = {};
        std::vector<int> bottom_crossing_components = {};
        int components;
        int writhe = 0;
        int prefactors;
        int crossings;
        int degree;
        void f (const std::vector<int>& angles) {
            for (int i = 0; i < crossings + 1; i++) {
                for (int j = 0; j < prefactors + 1; j++) {
                    numerical_assignment[i][j] = dot(assignment[i][j], angles);
                }
            }
            // std::cout << "ININININ\n";
            // std::vector<int> point = {};
            // for (auto x : angles) {
            //     std::cout << x << " ";
            //     point.push_back(x);
            // }
            // point[0] = 1;
            // point[1] = 0;
            // std::cout << "\n" << crossings << "\n";
            // std::cout << prefactors << "\n";
            // for (int i = 0; i < crossings + 1; i++) {
            //     for (int j = 0; j < prefactors + 1; j++) {
            //         std::cout << dot(assignment[i][j], point) << " ";
            //     }
            //     std::cout << "separatus";
            //     std::cout << std::endl;
            // }
            // exit(0);
            // std::cout << "\n";
            // for (auto x : numerical_assignment) {
            //     for (auto y : x) {
            //         std::cout << y << " ";
            //     }
            //     std::cout << std::endl;
            // }
            // std::cout << std::endl;
            // for (auto x : matrices) {
            //     std::cout << x[1][1] << " ";
            // }
            // std::cout << std::endl;
            // for (auto x : r) {
            //     std::cout << x << " ";
            // }
            // std::cout << std::endl;
            // // exit(0);
            // std::cout << "The below " << "\n";
            // for (auto x : assignment[0][1]) {
            //     std::cout << x << " ";
            // }
            // std::cout << "\n";
            // std::cout << numerical_assignment[0][1] << "\n";

            // for (auto x : assignment) {
            //     for (auto y : x) {
            //         for (auto z : y) {
            //             std::cout << z << " ";
            //         }
            //         std::cout << std::endl;
            //     }
            //     std::cout << std::endl;
            // }
            // exit(0);

            double q_acc_double = (writhe - prefactors) / 2.0;
            std::vector<double> x_acc_double(components, 0);
            int init = 1;
            for (int i = 0; i < prefactors; i++) {
                q_acc_double -= numerical_assignment[0][i + 1];
                x_acc_double[closed_strand_components[i]] -= 0.5;
            }

            for (int ind = 0; ind < crossings; ind++) {
                int i = numerical_assignment[matrices[ind][0][0]][matrices[ind][0][1]];
                int k = numerical_assignment[matrices[ind][2][0]][matrices[ind][2][1]];
                int j = numerical_assignment[matrices[ind][1][0]][matrices[ind][1][1]];  
                int m = numerical_assignment[matrices[ind][3][0]][matrices[ind][3][1]];  
                int t = top_crossing_components[ind];
                int b = bottom_crossing_components[ind];
                if (r[ind] == 1 || r[ind] == 2) {                                  
                    q_acc_double += (j + m) / 2.0 + j * m;
                    x_acc_double[t] += ((j + k + 1) / 4.0);
                    x_acc_double[b] += ((3 * m - i + 1) / 4.0);
                }
                else if (r[ind] == 4) { 
                    q_acc_double -= (i + k + m * (m + 1) - i * (i + 1)) / 2.0 + i * k;
                    x_acc_double[t] -= ((2 * j - i + m + 1) / 4.0);
                    x_acc_double[b] -= ((i + m + 1) / 4.0);
                    if ((j - k) % 2 == 0) {
                        init *= -1;
                    }
                }
                else if (r[ind] == 3) {  
                    q_acc_double -= (i + k + m * (m + 1) - i * (i + 1)) / 2.0 + i * k;
                    x_acc_double[t] -= ((3 * j - k + 1) / 4.0);
                    x_acc_double[b] -= ((i + m + 1) / 4.0);
                    if ((k - j) % 2 == 1) {
                        init *= -1;
                    }
                }
            }
            int q_acc = q_acc_double;
            std::vector<int> x_acc(components);
            std::vector<int> X_MAX(components);
            std::vector<int> blocks(components);
            blocks[0] = 1;
            for (int n = 0; n < components; n++) {
                x_acc[n] = x_acc_double[n];
                X_MAX[n] = degree - x_acc[n];
                if (n != 0) {
                    blocks[n] = (X_MAX[n - 1] + 1) * blocks[n - 1];
                }
            }
            int prod = blocks[components - 1] * (X_MAX[components - 1] + 1);
            // std::cout << "x_acc:" << x_acc_double[0] << "\n\n";
            std::cout << "x_acc:" << x_acc_double[0] << " "<< x_acc_double[1] << "\n\n"; // modifying: x_acc's are coming out to half-integers in general
            std::vector<bilvector<int>> term(prod,  bilvector<int>(0, 1, 20, 0)); // error is the degree issue; modifying
            term[0][0] = init;
            for (int ind = 0; ind < crossings; ind++) {
                if (r[ind] == 1) {
                    int i = numerical_assignment[matrices[ind][0][0]][matrices[ind][0][1]];
                    int m = numerical_assignment[matrices[ind][3][0]][matrices[ind][3][1]];
                    if (i > 0) {
                        pp_q_binom(term, i, i - m, false);
                    }
                    else {
                        np_q_binom(term, i, i - m, false);
                    }
                }
                else if (r[ind] == 2) {
                    int i = numerical_assignment[matrices[ind][0][0]][matrices[ind][0][1]];
                    int m = numerical_assignment[matrices[ind][3][0]][matrices[ind][3][1]];
                    np_q_binom(term, i, m, false);
                }
                else if (r[ind] == 3) {      
                    int j = numerical_assignment[matrices[ind][1][0]][matrices[ind][1][1]];
                    int k = numerical_assignment[matrices[ind][2][0]][matrices[ind][2][1]];
                    np_q_binom(term, j, k, true);
                }
                else {
                    int j = numerical_assignment[matrices[ind][1][0]][matrices[ind][1][1]];
                    int k = numerical_assignment[matrices[ind][2][0]][matrices[ind][2][1]];
                    if (j > 0) {
                        pp_q_binom(term, j, j - k, true);
                    }
                    else {
                        np_q_binom(term, j, j - k, true);
                    }
                }
            }
            for (int ind = 0; ind < crossings; ind++) {
                if (r[ind] == 1) {
                    int b = bottom_crossing_components[ind];
                    int j = numerical_assignment[matrices[ind][1][0]][matrices[ind][1][1]];
                    int k = numerical_assignment[matrices[ind][2][0]][matrices[ind][2][1]];
                    x_q_pochhammer(term, k, j + 1, b, components, X_MAX, blocks);
                }
                else if (r[ind] == 2) {
                    int b = bottom_crossing_components[ind];
                    int j = numerical_assignment[matrices[ind][1][0]][matrices[ind][1][1]];
                    int k = numerical_assignment[matrices[ind][2][0]][matrices[ind][2][1]];
                    x_q_inv_pochhammer(term, j, k + 1, b, components, X_MAX, blocks);
                }
                else if (r[ind] == 3) {  
                    int t = top_crossing_components[ind];
                    int i = numerical_assignment[matrices[ind][0][0]][matrices[ind][0][1]];
                    int m = numerical_assignment[matrices[ind][3][0]][matrices[ind][3][1]];    
                    x_q_inv_pochhammer(term, i, m + 1, t, components, X_MAX, blocks);
                }
                else {
                    int t = top_crossing_components[ind];
                    int i = numerical_assignment[matrices[ind][0][0]][matrices[ind][0][1]];
                    int m = numerical_assignment[matrices[ind][3][0]][matrices[ind][3][1]];
                    x_q_pochhammer(term, m, i + 1, t, components, X_MAX, blocks);
                }
            } 
            // for (int w = 0; w < X_MAX[0] + 1; w++) {
            //     for (int j = term[w].get_max_nindex(); j <= term[w].get_max_pindex(); j++) {
            //         std::cout << term[w][j] << " ";
            //     }
            //     std::cout << std::endl;
            // } 
            // std::cout << "q_acc : " << q_acc << "\n";
            offset_addition(acc, term, x_acc, q_acc, components, X_MAX, 1, acc_blocks, blocks);
            // std::cout << "acc : \n";
            // for (int w = 0; w < degree + 1; w++) {
            //     for (int j = acc[w].get_max_nindex(); j <= acc[w].get_max_pindex(); j++) {
            //         std::cout << acc[w][j] << " ";
            //     }
            //     std::cout << std::endl;
            // } 
            // std::cout << "acc completo\n";
        }
        void write (std::string file_) {
            std::ofstream file;
            file.open(file_ + ".json");
            file << "{\n\t\"coefficient_q_powers\":[\n";
            for (int i = 0; i < acc.size(); i++) {
                file << "\t\t[";
                bool first_write = true;
                for (int j = acc[i].get_max_nindex(); j <= acc[i].get_max_pindex(); j++) {
                    if (acc[i][j] != 0) {
                        if (!first_write) {
                            file << ",[";
                        }
                        else {
                            file << "[";
                            first_write = false;
                        }
                        file << j << "," << acc[i][j] << "]";
                    }
                }
                if (i < acc.size() - 1) {
                    file << "],\n";
                }
                else {
                    file << "]\n";
                }
            }
            file << "\t]\n}";
            file.close();
        }
    public:
        std::vector<std::vector<double>> inequalities;
        std::vector<std::vector<double>> criteria;
        std::vector<std::vector<int>> extensions;
        std::string metadata;
        FK (std::string infile_, std::string outfile_) {

            std::ifstream infile;
            infile.open(infile_ + ".csv");
            if (infile.is_open())
            {
                std::string line;

                std::getline (infile,line,'\n');
                int index = line.find(",");
                degree = string_to_int(line.substr(0, index));

                std::getline (infile,line,'\n');
                index = line.find(",");
                components = string_to_int(line.substr(0, index));

                std::getline (infile,line,'\n');
                index = line.find(",");
                writhe = string_to_int(line.substr(0, index));

                acc.resize(std::pow(degree + 1, components), bilvector<int>(0, 1, 20, 0));

                std::getline (infile,line,'\n');
                int height = 0;
                index = line.find(",");
                bool done = false;
                while (true) {
                    if (index == -1) {
                        break;
                    }
                    
                    int c = string_to_int(line.substr(0, index));
                    matrices.push_back({
                        {height, c - 1}, 
                        {height, c},
                        {height + 1, c - 1},
                        {height + 1, c}
                    });
                    line = line.substr(index + 1, line.size() - index - 1);
                    index = line.find(",");

                    r.push_back(string_to_int(line.substr(0, index)));
                    line = line.substr(index + 1, line.size() - index - 1);

                    height++;
                    index = line.find(",");
                }
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
                prefactors = closed_strand_components.size();
                crossings = r.size();
                assignment.resize(crossings + 1, std::vector<std::vector<int>>(prefactors + 1));
                numerical_assignment.resize(crossings + 1, std::vector<int>(prefactors + 1));
                std::getline (infile,line,'\n');
                index = line.find(",");
                done = false;
                while (true) {
                    if (index == -1) {
                        break;
                    }
                    top_crossing_components.push_back(string_to_int(line.substr(0, index)));
                    line = line.substr(index + 1, line.size() - index - 1);
                    index = line.find(",");
                    bottom_crossing_components.push_back(string_to_int(line.substr(0, index)));
                    line = line.substr(index + 1, line.size() - index - 1);
                    index = line.find(",");
                }
                int stage = 0;
                int criteria_index = 0;
                int inequality_index = 0;
                int extension_index = 0;
                while ( std::getline (infile,line,'\n') )
                {
                    if (line[0] == '/') {
                        stage++;
                    }
                    else if (stage == 0) {
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
                    else if (stage == 1) {
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
                    else if (stage == 2) {
                        done = false;
                        index = line.find(",");
                        while (true) {
                            if (index == -1) {
                                break;
                            }
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
            for (int i = 1; i < components; i++) {
                acc_blocks.push_back(acc_blocks[i - 1] * (degree + 1));
            }
            std::function<void(const std::vector<int>&)> f_wrapper = [this](const std::vector<int>& v) { f(v); };
            pooling(criteria, inequalities, f_wrapper);
            std::vector<int> increment_offset(components);
            increment_offset[0] = 1;
            std::vector<int> maxima(components, degree - 1);
            offset_addition(acc, acc, increment_offset, 0, components, maxima, -1, acc_blocks, acc_blocks);
            write(outfile_);
        }
};

int main(int argc, char* argv[]) {

    // std::cout << argv[1] << "\n";

    // std::string braid_string = argv[1]; 
    // std::vector<int> braid = {};
    // int index;
    // while (true) {
    //     index = braid_string.find(",");
    //     if (index == -1) {
    //         break;
    //     }
    //     braid.push_back(string_to_int(braid_string.substr(0, index)));
    //     braid_string = braid_string.substr(index + 1, braid_string.size() - index - 1);
    // }
    // std::string inversion_data_string = argv[2]; 
    // std::vector<int> inversion_data = {};
    // index = 0;
    // while (true) {
    //     if (index == -1) {
    //         break;
    //     }
    //     index = inversion_data_string.find(",");
    //     inversion_data.push_back(string_to_int(inversion_data_string.substr(0, index)));
    //     inversion_data_string = inversion_data_string.substr(index + 1, inversion_data_string.size() - index - 1);
    // }
    // int degree = string_to_int(argv[3]);
    // std::string file = argv[4];
    // for (int x: braid) {
    //     std::cout << x << "\n";
    // }
    // FK(braid, inversion_data, degree, file);
    FK("davide/A", "davide/A");
    FK("davide/B", "davide/B");
    FK("davide/C", "davide/C");
    FK("davide/D", "davide/D");
}

// implement multithreading at top level of recursion