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

#include "opts/opts.h"

//#define LOG_AUCTION

//#include "auction_runner_fr.h"
//#include "auction_runner_fr.hpp"

#include "wasserstein.h"

// any container of pairs of doubles can be used,
// we use vector in this example.


using IdxType = int;

int main(int argc, char* argv[])
{
    using PairVector = std::vector<std::pair<double, double>>;
    PairVector diagramA, diagramB;

    std::vector<PairVector> diagrams(3);

    hera::AuctionParams<double> params;
    params.max_num_phases = 800;

    opts::Options ops(argc, argv);
    ops >> opts::Option('q', "degree", params.wasserstein_power, "Wasserstein degree")
        >> opts::Option('d', "error", params.delta, "Relative error")
        >> opts::Option('p', "internal-p", params.internal_p, "Internal norm")
        >> opts::Option("initial-epsilon", params.initial_epsilon, "Initial epsilon")
        >> opts::Option("epsilon-factor", params.epsilon_common_ratio, "Epsilon factor")
        >> opts::Option("max-bids-per-round", params.max_bids_per_round, "Maximal number of bids per round")
        >> opts::Option('m', "max-rounds", params.max_num_phases, "Maximal number of iterations");


    //Use the 2-Wasserstein distance
    params.wasserstein_power = 2;
    //Use Euclidean norm for the distances pApB with pA, pB points in persistence diagrams
    params.internal_p = 2;


    bool print_relative_error = ops >> opts::Present('e', "--print-error", "Print real relative error");

    params.tolerate_max_iter_exceeded  = ops >> opts::Present('t', "tolerate", "Suppress max-iterations-exceeded error and print the best result.");

    std::cout << "q = " << params.wasserstein_power << ", delta = " << params.delta << ", p = " << params.internal_p << ", max_round = " << params.max_num_phases <<  std::endl;
    std::cout << "print relative error: " << print_relative_error << std::endl;



    //Read the diagrams
    hera::read_diagram_point_set("tests/data/sample/test_5_A", diagrams[0]);
    hera::read_diagram_point_set("tests/data/sample/test_5_B", diagrams[1]);
    hera::read_diagram_point_set("tests/data/sample/test_5_C", diagrams[2]);

    //Initialize the average diagram
    PairVector Fravg = diagrams[0];

    //The size of the diagram
    int diag_size = Fravg.size();

    if (diag_size == 0) {
        std::cout << "Empty first diagram." << std::endl;
        return 0;
    }

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
    int max_ITER = 10;

    while(!done) {
        //Reset the vectors
        for (int i = 0; i < diag_size; i++) {
            arithmetic_mean[i].first = 0;
            arithmetic_mean[i].second = 0;
            off_diagonal_matches[i] = 0;
        }


        //Calculate the optimal matchings
        for (int i = 0; i < diagrams.size(); i++) {
            //Vectors to store the pairing and costs
            std::vector<IdxType> bidders_to_items;
            std::vector<double> edge_costs;

            double distance = hera::wasserstein_dist_with_pairings(Fravg, diagrams[i], params, bidders_to_items, edge_costs);

            for (int j = 0; j < diag_size; j++) {
                //If matched with off-diagonal point
                if (bidders_to_items[j] >= diag_size) {
                    off_diagonal_matches[j]++;
                    arithmetic_mean[j].first += diagrams[i][bidders_to_items[j]-diag_size].first;
                    arithmetic_mean[j].second += diagrams[i][bidders_to_items[j]-diag_size].second;
                }
            }
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
                Fravg[i].first = (w_x * k + delta_x * (diag_size-k)) / diag_size;
                Fravg[i].second = (w_y * k + delta_y * (diag_size-k)) / diag_size;

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

        std::cout << "Total distance moved = " << dist_moved << std::endl;

        ITER++;
        if (ITER>max_ITER) {
            done = true;
            std::cout << "Exceeded max number of loops." << std::endl;
            return 0;
        }
    }


    //Determine the Frechet variances associated with each point
    std::vector<double> variances;

    //Initialize
    variances.resize(diag_size);
    for (int i = 0; i < diag_size; i++) {
        variances[i] = 0;
    }

    //Calculate the optimal matchings
    for (int i = 0; i < diagrams.size(); i++) {
        //Vectors to store the pairing and costs
        std::vector<IdxType> bidders_to_items;
        std::vector<double> edge_costs;

        double distance = hera::wasserstein_dist_with_pairings(Fravg, diagrams[i], params, bidders_to_items, edge_costs);

        for (int j = 0; j < diag_size; j++) {
            double x = Fravg[j].first;
            double y = Fravg[j].second;

            //If matched with off-diagonal point
            if (bidders_to_items[j] >= diag_size) {
                double x2 = diagrams[i][bidders_to_items[j]-diag_size].first;
                double y2 = diagrams[i][bidders_to_items[j]-diag_size].second;

                variances[j] += ((x-x2)*(x-x2) + (y-y2)*(y-y2))/diag_size;
            } else {
                //Otherwise find the closest point on the diagonal
                double x2 = 0.5*(x+y);
                double y2 = 0.5*(x+y);


                variances[j] += ((x-x2)*(x-x2) + (y-y2)*(y-y2))/diag_size;
            }
        }

                std::cout << bidders_to_items.size() << std::endl;
    }

    //Print the results
    std::cout << "#x\ty\tsigma" << std::endl;

    for (int i = 0; i < diag_size; i++) {
        std::cout << Fravg[i].first << "\t" << Fravg[i].second << "\t" << sqrt(variances[i]) << std::endl;
    }

    return 0;

}
