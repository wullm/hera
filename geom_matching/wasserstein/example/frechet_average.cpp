/*
Original work Copyright (c) 2015, M. Kerber, D. Morozov, A. Nigmetov
Modified work Copyright (c) 2019, Willem Elbers
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
(Enhancements) to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to copyright holder,
without imposing a separate written license agreement for such Enhancements,
then you hereby grant the following license: a  non-exclusive, royalty-free
perpetual license to install, use, modify, prepare derivative works, incorporate
into other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.

*/

#include <iostream>
#include <locale>
#include <iomanip>
#include <vector>

#include <experimental/filesystem> //to list files in directory

#include "opts/opts.h"

#include "wasserstein.h"

using IdxType = int;

namespace fs = std::experimental::filesystem;

int main(int argc, char* argv[])
{
    //Type used for diagrams
    using PairVector = std::vector<std::pair<double, double>>;

    //The vector of diagrams
    std::vector<PairVector> diagrams;

    //Parameters for the bidding algorithm
    hera::AuctionParams<double> params;

    //Maximum number of bidding rounds
    params.max_num_phases = 800;
    //Use the 2-Wasserstein distance
    params.wasserstein_power = 2;
    //Use Euclidean norm for the distances pApB with pA, pB points in persistence diagrams
    params.internal_p = 2;
    //The Frechet algorithm justifies a very small delta in the Wasserstein step, necessary for the gradient descent
    params.delta = 1e-7;

    //The number of random starting points to try, in finding the best Frechet average
    int starting_points = 1;
    //Maximum number of iterations in the Frechet average diagram
    int max_ITER = 1000;

    //Load user options
    opts::Options ops(argc, argv);
    ops >> opts::Option('q', "degree", params.wasserstein_power, "Wasserstein degree")
        >> opts::Option('d', "error", params.delta, "Relative error")
        >> opts::Option('p', "internal-p", params.internal_p, "Internal norm")
        >> opts::Option("initial-epsilon", params.initial_epsilon, "Initial epsilon")
        >> opts::Option("epsilon-factor", params.epsilon_common_ratio, "Epsilon factor")
        >> opts::Option("max-bids-per-round", params.max_bids_per_round, "Maximal number of bids per round")
        >> opts::Option('m', "max-rounds", params.max_num_phases, "Maximal number of iterations")
        >> opts::Option('s', "starting-points", starting_points, "Number of random starting points to try for the Frechet average")
        >> opts::Option("Fr-max-iter", max_ITER, "Maximum number of iterations in the Frechet averaging step");

    //This option applies both to the Wasserstein bidding & to the Frechet averaging
    params.tolerate_max_iter_exceeded  = ops >> opts::Present('t', "tolerate", "Suppress max-iterations-exceeded error and print the best result.");
    //Output steps?
    bool verbose  = ops >> opts::Present('v', "verbose", "Print step information in addition to the average diagram.");

    //Help message
    if (ops >> opts::Present('h', "help", "show help message")) {
        std::cerr << "Usage: " << argv[0] << " path " << std::endl;
        std::cerr << "Compute Frechet average of diagrams contained in directory at path.\n";
        std::cout << ops << std::endl;
        return 1;
    }

    //Print options
    if (verbose) {
        std::cout << "#options: q = " << params.wasserstein_power << ", delta = " << params.delta << ", p = " << params.internal_p << ", max_round = " << params.max_num_phases << ", starting_points = " << starting_points <<  std::endl;
    }

    //The path to the directory containing the diagrams
    std::string dir_path;
    bool dir_path_given = (ops >> opts::PosOption(dir_path));

    if (!dir_path_given) {
        std::cerr << "Usage: " << argv[0] << " path " << std::endl;
        std::cerr << "Compute Frechet average of diagrams contained in directory at path.\n";
        std::cerr << ops << std::endl;
        return 1;
    }

    //Open all the files in the input directory
    for (const auto & entry : fs::directory_iterator(dir_path)) {
        if (verbose) {
            std::cout << "#Reading " << entry.path() << std::endl;
        }
        PairVector read_diagram;
        hera::read_diagram_point_set(entry.path(), read_diagram);

        if (read_diagram.size() > 0) {
            diagrams.push_back(read_diagram);
        } else {
            std::cerr << "Ignoring file " << entry.path() << "." << std::endl;
        }
    }

    //Manual way of reading the diagrams
    // hera::read_diagram_point_set("tests/data/sample2/test_100_A", diagrams[0]);

    //The number of diagrams
    const int diag_num = diagrams.size();

    if (diag_num == 0) {
        std::cerr << "No non-empty diagrams." << std::endl;
        return 0;
    } else if (diag_num == 1) {
        std::cerr << "Only one non-empty diagram, average not needed." << std::endl;
        return 0;
    }


    if (verbose) {
        std::cout << "#Read " << diag_num << " diagrams from files in " << dir_path << std::endl;
    }

    PairVector best_Fravg;
    double lowest_cost = -1;

    //Seed the rng
    srand (time(NULL));

    /**
    * We can run the gradient descent-inspired algorithm for multiple random
    * starting points, to find the best minimum.
    **/

    for (int sptry=0; sptry<starting_points; sptry++) {

        /**
        * We want to initialize the average diagram with some guess,
        * we choose the midpoint of two random diagrams.
        **/

        //Choose two random diagrams
        int one = rand() % diag_num;
        int two = rand() % diag_num;

        //Make sure that the two diagrams are different
        while (one == two) {
            two = rand() % diag_num;
        }

        if (verbose) {
            std::cout << "#Intializing from midpoint of diagrams " << one << " and " << two << "." << std::endl;
        }

        //Set the size of the average diagram equal to that of the first random diagram
        PairVector Fravg;
        int diag_size = diagrams[one].size();
        if (diag_size == 0) {
            assert(false);
        }
        Fravg.resize(diag_size);

        //Calculate the pairing with the second random diagram
        std::vector<IdxType> ini_bidders_to_items;
        std::vector<double> ini_edge_costs;
        double ini_distance = hera::wasserstein_dist_with_pairings(diagrams[one], diagrams[two], params, ini_bidders_to_items, ini_edge_costs);

        //Update the points, setting them to the midpoints of the two matches
        for (int j = 0; j < diag_size; j++) {
            //If the point is matched with an off-diagonal point in diagrams[i]
            if (ini_bidders_to_items[j] >= diag_size) {
                Fravg[j].first = 0.5*(diagrams[one][j].first + diagrams[two][ini_bidders_to_items[j]-diag_size].first);
                Fravg[j].second = 0.5*(diagrams[one][j].second + diagrams[two][ini_bidders_to_items[j]-diag_size].second);
            } else {
                //Matched with projection onto diagonal of this point
                double off_x = diagrams[one][ini_bidders_to_items[j]].first;
                double off_y = diagrams[one][ini_bidders_to_items[j]].second;

                //The projection
                double diag_x = 0.5*(off_x + off_y);
                double diag_y = 0.5*(off_x + off_y);

                if (off_x != diagrams[one][j].first || off_y != diagrams[one][j].second) {
                    assert(false);
                }

                //Set at the midpoint
                Fravg[j].first = 0.5*(diagrams[one][j].first + diag_x);
                Fravg[j].second = 0.5*(diagrams[one][j].second + diag_y);
            }
        }


        /**
        * End of initialization
        **/


        //Initialize the vectors used to calculate the arithmetic means of the off-
        //diagonal matches
        std::vector<int> off_diagonal_matches;
        std::vector<std::pair<double, double>> arithmetic_mean;

        //Resize the vectors
        arithmetic_mean.resize(diag_size);
        off_diagonal_matches.resize(diag_size);

        //Do the main loop
        bool done = false;
        int ITER = 0;
        double total_cost;

        while(!done) {
            //Reset the vectors
            for (int i = 0; i < diag_size; i++) {
                arithmetic_mean[i].first = 0;
                arithmetic_mean[i].second = 0;
                off_diagonal_matches[i] = 0;
            }

            total_cost = 0;

            //Calculate the optimal matchings
            for (int i = 0; i < diag_num; i++) {
                //Vectors to store the pairing and costs
                std::vector<IdxType> bidders_to_items;
                std::vector<double> edge_costs;

                /**
                * The edges are recorded in the vector bidders_to_items.
                * Each edge connects a point in X = D_avg U D_i' with a point in
                * Y = D_avg' U D_i. Here, D_avg and D_i are the two persistence diagrams
                * with m,n off-diagonal points respectively. And D_avg' and D_i'
                * are the sets of orthogonal projections of points (x,y) onto
                * the diagonal at 0.5*(x+y, x+y), also with m,n points respectively.
                * Hence, both X and Y have m+n points and there are m+n edges.
                *
                * We only need to record the edges between off-diagonal points in
                * D_avg and off-diagonal points in D_i, for each i. These correspond
                * to the edges (j,k) with j < m and k = bidders_to_items[j] >= m.
                **/

                double distance = hera::wasserstein_dist_with_pairings(Fravg, diagrams[i], params, bidders_to_items, edge_costs);

                //For each off-diagonal point in D_avg
                for (int j = 0; j < diag_size; j++) {
                    //If the point is matched with an off-diagonal point in diagrams[i]
                    if (bidders_to_items[j] >= diag_size) {
                        off_diagonal_matches[j]++;
                        arithmetic_mean[j].first += diagrams[i][bidders_to_items[j]-diag_size].first;
                        arithmetic_mean[j].second += diagrams[i][bidders_to_items[j]-diag_size].second;
                    }
                }

                total_cost += distance*distance;
            }

            //Calculate the total distance moved
            double move_x = 0;
            double move_y = 0;

            //Update the points
            for (int i = 0; i < diag_size; i++) {
                double old_x = Fravg[i].first;
                double old_y = Fravg[i].second;

                int k = off_diagonal_matches[i];
                if (k > 0) {
                    //The arithmetic mean of the off-diagonal matches
                    double w_x = arithmetic_mean[i].first / k;
                    double w_y = arithmetic_mean[i].second / k;

                    //The closest diagonal point to (w_x, w_y)
                    double delta_x = 0.5*(w_x + w_y);
                    double delta_y = 0.5*(w_x + w_y);

                    //Update the point location to the weighted average of w and delta
                    Fravg[i].first = (w_x * k + delta_x * (diag_num-k)) / diag_num;
                    Fravg[i].second = (w_y * k + delta_y * (diag_num-k)) / diag_num;

                    //Maintain the move distance
                    move_x += Fravg[i].first - old_x;
                    move_y += Fravg[i].second - old_y;
                } else { //only diagonal matches
                    //Just put the point equal to a diagonal point
                    Fravg[i].second = Fravg[i].first ;

                    //Maintain the move distance
                    move_y += Fravg[i].second - old_y;
                }
            }

            double dist_moved = sqrt(move_x*move_x + move_y*move_y);

            if (dist_moved == 0) {
                done = true;
            }

            if (verbose) {
                std::cout << ITER << ". Total distance moved = " << dist_moved << ", total cost = " << total_cost << ", error = " << params.delta << std::endl;
            }

            ITER++;
            if (ITER>max_ITER && !params.tolerate_max_iter_exceeded) {
                done = true;
                std::cerr << "Exceeded max number of loops." << std::endl;
                return 0;
            }
        }

        if (sptry == 0 || total_cost < lowest_cost) {
            lowest_cost = total_cost;
            best_Fravg = Fravg;
        }
    }



    /**
    * Determine the Frechet variances associated with each point
    **/
    std::vector<double> variances;

    double total_variance = 0;

    //Initialize
    int diag_size = best_Fravg.size();
    variances.resize(diag_size);
    for (int i = 0; i < diag_size; i++) {
        variances[i] = 0;
    }

    //Calculate the optimal matchings
    for (int i = 0; i < diag_num; i++) {
        //Vectors to store the pairing and costs
        std::vector<IdxType> bidders_to_items;
        std::vector<double> edge_costs;

        double distance = hera::wasserstein_dist_with_pairings(best_Fravg, diagrams[i], params, bidders_to_items, edge_costs);

        total_variance += distance*distance;

        for (int j = 0; j < diag_size; j++) {
            double x = best_Fravg[j].first;
            double y = best_Fravg[j].second;

            //If matched with off-diagonal point
            if (bidders_to_items[j] >= diag_size) {
                double x2 = diagrams[i][bidders_to_items[j]-diag_size].first;
                double y2 = diagrams[i][bidders_to_items[j]-diag_size].second;

                variances[j] += ((x-x2)*(x-x2) + (y-y2)*(y-y2))/diag_num;
            } else {
                //Otherwise find the closest point on the diagonal
                double x2 = 0.5*(x+y);
                double y2 = 0.5*(x+y);

                variances[j] += ((x-x2)*(x-x2) + (y-y2)*(y-y2))/diag_num;
            }
        }
    }

    //Print the results
    std::cout << "#total variance = " << total_variance << std::endl;
    std::cout << "#x\ty\tsigma" << std::endl;

    for (int i = 0; i < diag_size; i++) {
        std::cout << best_Fravg[i].first << "\t" << best_Fravg[i].second << "\t" << sqrt(variances[i]) << std::endl;
    }

    return 0;

}
