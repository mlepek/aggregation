/******************************************
The code for the numerical simulations of coagulating systems
Author: Michal Lepek
Date: 02 Dec 2019
If using this code for your research, please refer to the main work:
"Combinatorial solutions to coagulation kernel for linear chains"
by M. Lepek, A. Fronczak & P. Fronczak, Physica D: Nonlinear Phenomena, Volume 415, January 2021, 132756
******************************************/

/// *********** Useful notes:
// To simulate the aggregating system, you need:
// 1. Define number of simulations to be performed and averaged;
// 2. Define number of monomeric units in the system (N);
// 3. Define time of observation (t_max);
// 4. and define the kernel (rule). All this you set in the main function.
//
// The rules can be chosen from the list of my predefined kernels, or you can define you own kernels.
// The list of kernels is given in the main function.
//
// This code will produce two output files: one with averaged number of clusters of different sizes, <ns>,
// and the second, with standard deviation of these numbers derived from simulated data.
// Standard deviation was not tested for the number of simulations higher than 10^6.
//
// "Mer" is a synonyme of "cluster".

#include<iostream>
#include<vector>
#include<cstdlib>
#include<cstdio>
#include<ctime>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<valarray>
#include<random>
#include<string>

#include "functions.hpp"

using namespace std;

typedef double (*PointerToKernel)(int,int);

int aggregation(const int N, const int t_max, const int number_of_simulations, const string rule, const int);

int main()
{
    // Define initial conditions
    int initial_conditions = 1; // 1 - monodisperse, 2 - exponential, 3 - power-law

    // Define the number of simulations to be averaged.
    int number_of_simulations = 100000;

    // Define the kernel.
    string rule = "constant_fast";
    /* You can choose several rules listed here:

     constant_fast
     additive_fast
     product_fast
     product
     constant
     additive
     electrorheological
     condensation
     electrorheological_alpha (*)
     additive_with_one
     antiproduct
     antiadditive
     kernel_X (K=ij/(i+j)+ij)

    Rules that are called "_fast" are simulated as choosing appropriate elements from vectors.
    This method is fast but not applicable to arbitrary kernel.

    Rules with no "_fast" identifier are simulated using the "generalized real number vector".
    You can read about this in the Appendix of the main paper.
    This method is slower but it can be applied to any form of kernel.
    Thus, you can define your own kernel, and do not worry of "how to simulate it".
    Kernels that I defined before, are in lines 94-135.

    (*) In case of the electrorheological_alpha kernel, to specify alpha you need to change the float number
    in the function definition (line 118).
    */

    // using this for loop you can perform several simulations at once, one by one
    for(int i=1; i<=1; i++)
    {
        int N = 1000; // number of monomeric units in the system
        int t_max = 950; // time of simulation - no more than N; k = N - t_max

        //aggregation(N, t_max, number_of_simulations, rule, initial_conditions);
        //rule = "additive_fast";
        aggregation(N, t_max, number_of_simulations, rule, 1);
        aggregation(N, t_max, number_of_simulations, rule, 2);
        aggregation(N, t_max, number_of_simulations, rule, 3);
        //rule = "product_fast";
        //aggregation(N, t_max, number_of_simulations, rule, 1);
        aggregation(N, t_max, number_of_simulations, rule, 2);
        aggregation(N, t_max, number_of_simulations, rule, 3);

    }

    return 0;
}

// Here I defined several kernel functions.

inline double constant(int i, int j)
{
    return 1;
}
inline double product(int i, int j)
{
    return i * j;
}
inline double additive(int i, int j)
{
    return i + j;
}
inline double electrorheological(int i, int j)
{
    return 1./i + 1./j;
}
inline double condensation(int i, int j)
{
    return (i + 1) * (j + 1);
}
inline double electrorheological_alpha(int i, int j)
{
    return pow(1./i + 1./j, -1.0);
}
inline double additive_with_one(int i, int j)
{
    return 1 + i + j;
}
inline double antiproduct(int i, int j)
{
    return 1./(i*j);
}
inline double kernel_X(int i, int j)
{
    return i*j*( 1./(i+j) + 1 );
}
inline double antiadditive(int i, int j)
{
    return 1./(0+i+j);
}

// This function simulates the aggregation process
int aggregation(const int N, const int t_max, const int number_of_simulations, const string rule, const int initial_conditions)
{
    clock_t begin = clock();

    PointerToKernel pKernel;

    // Here, I choose the kernel that was specified by the user.
    // If one of the "fast" rules specified, then I do nothing.
    // If other rule specified, then I set pointer to the specified kernel function.

    if( rule == "constant_fast" or rule == "product_fast" or rule == "additive_fast" )
        ;
    else if( rule == "constant" )
        pKernel = constant;
    else if( rule == "product" )
        pKernel = product;
    else if( rule == "additive" )
        pKernel = additive;
    else if( rule == "electrorheological" )
        pKernel = electrorheological;
    else if( rule == "condensation" )
        pKernel = condensation;
    else if( rule == "electrorheological_alpha" )
        pKernel = electrorheological_alpha;
    else if( rule == "additive_with_one" )
        pKernel = additive_with_one;
    else if( rule == "antiproduct" )
        pKernel = antiproduct;
    else if( rule == "kernel_X" )
        pKernel = kernel_X;
    else if( rule == "antiadditive" )
        pKernel = antiadditive;
    else
    {
        cout << "Kernel you have chosen is not implemented" << endl;
        return 1;
    }

    // strings for output files
    string output_file_string = "N=" + to_string(N) + "_t=" + to_string(t_max) + "_ic=" + to_string(initial_conditions) +
                                + "_s=" + to_string(number_of_simulations) + "_" + rule;
    string output_file_avns_string = output_file_string + "_av_ns.txt";
    string output_file_stddev_string = output_file_string + "_std_dev.txt";

	// initialization of the random number generator
	default_random_engine generator;
	int maxInt = 2147483647;
	uniform_int_distribution<int> distribution(0, maxInt);

	// for random seed
	generator.seed( time(NULL) );

	// container for saving cluster distributions obtained from simulations
	vector<IVector> container_cluster_distributions;

	// container for partially averaged numbers of clusters
	// this is used to not to overload RAM if you set number of simulations > 10^6
    vector<DVector> container_averaged_number_of_clusters;


	// This loop performs specified number of simulations.
	for (int simulation = 0; simulation < number_of_simulations; simulation++)
	{


        if (simulation % 100 == 0)
			cout << "N = " << N << ", simulation " << simulation + 1 << ":  " << double(simulation) / number_of_simulations * 100. << " %\t"
                 << "Time elapsed: " << double(clock() - begin) / CLOCKS_PER_SEC << " s" << endl;


        // This vector stores labels (ordinal numbers) of mers (clusters).
        // Size of this vector decreases in time as the number of mers decreases.
        IVector labels_of_mers;

        // This vector contains masses (sizes) of mers. Its length decreases in time.
        IVector masses_of_mers;

        // This vector contains the information on: which cluster (mer) does this monomer belongs to.
        // This vector is of constant length.
        IVector membership_of_monomers;

        /////////////////////////////////////// initial conditions ////////////////////////////////////////

        if(initial_conditions == 1) // monodisperse initial conditions
        {
            // labels from 1 to N
            labels_of_mers = IVector(N);
            for (int i = 0; i<N; ++i)
                labels_of_mers[i] = i + 1;

            // all masses of size 1
            masses_of_mers = IVector(N, 1);

            // all are distinct clusters
            membership_of_monomers = IVector(N);
            for (int i = 0; i<N; ++i)
                membership_of_monomers[i] = i + 1;
        }
        else if(initial_conditions == 2 or initial_conditions == 3 ) // exponential or power-law initial conditions
        {
            // let's prepare distribution for random choosing

            valarray<double> values_x(N/100+1);
            double step = 1./(N/100+1);
            for(int i=0; i<values_x.size(); i++)
                values_x[i] = step + step*i;

            double lambda = 15.0; // for exponential initial conditions, 5, 15
            double alpha = 2.0; // for power-law initial conditions, 1, 2

            valarray<double> values_y;

            // exponential initial conditions
            if(initial_conditions == 2)
                values_y = exp( -lambda * values_x ); // <<< exponential function
            // power-law initial conditions
            else if(initial_conditions == 3)
                values_y = pow(values_x, -alpha); // <<< power function


            // normalization so the first element shall be 1
            double a = 1. / values_y[0];
            values_y = values_y * a;

            //for(int i=0; i<values_y.size(); i++)
            //    cout << values_y[i] << endl;

            // no we have distribution ready (span_y = 1, span_x = 1)

            int sum_of_masses = 0;

            // now, we have to randomly choose a number of clusters, their sizes from the distribution, the total number of monomers equal to N
            while(true)
            {
                // choosing two random floats from 0 to 1
                int cluster_size = distribution(generator) % (N/100+1) + 1; // random int from 1 to N/100
                double random_double = double( distribution(generator) ) / maxInt;

                //cout << random_double << "  " << cluster_size << "  " << values_y[cluster_size] << endl;
                //getchar();

                // testing the distribution
                if( random_double < values_y[cluster_size] )
                    ; // accept
                else
                    continue; // reject

                // lets count how many aggregation events (steps) were needed to construct such a configuration
                int time_steps_needed = 0;
                for(int i=0; i<masses_of_mers.size(); i++)
                        time_steps_needed  += ( masses_of_mers.at(i) - 1 );
                time_steps_needed += (cluster_size-1);


                // I have to add mass of this new cluster to know the sum of masses (it must be equal to N).
                // I must make sure that the sum of masses of clusters will not be greater than N
                // and that number of steps needed to create such a configuartion will not be greater than lets say N/4.
                // The latter is important to compare different simulation together.
                if( (sum_of_masses + cluster_size > N) or (time_steps_needed > N*0.7 ) )
                    continue;
                else
                {
                    sum_of_masses += cluster_size;

                    // adding new cluster in labels_of_mers
                    int current_ordinal_number_of_cluster = labels_of_mers.size() + 1;
                    labels_of_mers.push_back( current_ordinal_number_of_cluster );

                    // adding new cluster in masses_of_mers
                    masses_of_mers.push_back(cluster_size);

                    // adding new cluster in membership_of_monomers;
                    for(int i=0; i<cluster_size; i++)
                        membership_of_monomers.push_back( current_ordinal_number_of_cluster );
                }

                // if sum of the masses is N then stop, if not then go again and generate
                if(sum_of_masses == N)
                    break;
            } // end of while

            /*cout << "labels_of_mers.size = " << labels_of_mers.size() << endl;
            cout << "masses_of_mers.size = " << masses_of_mers.size() << endl;
            cout << "membership_of_monomers.size = " << membership_of_monomers.size() << endl;
            for(int i=0; i<membership_of_monomers.size(); i++)
                cout << membership_of_monomers.at(i) << " " ;
            cout << endl;
            cout << endl;
            for(int i=0; i<masses_of_mers.size(); i++)
                cout << masses_of_mers.at(i) << " " ;
            cout << endl;*/
        }

        // test output of initial conditions
        /*IVector initial_cluster_distribution = histogram(masses_of_mers, N);
        ofstream test_output("test_output.txt");
        for(int i=0; i<initial_cluster_distribution.size(); i++)
            test_output << initial_cluster_distribution.at(i) << endl;
        test_output.close();
        getchar();/**/

        /////////////////////////////////////// end of initial conditions ////////////////////////////////////////


        // setting preliminary t's value
        int starting_t_value = 0;
        switch(initial_conditions)
        {
            case 1: starting_t_value = 1;
                    break;

            case 2: for(int i=0; i<masses_of_mers.size(); i++)
                        starting_t_value += ( masses_of_mers.at(i) - 1 );
                    break;

            case 3: for(int i=0; i<masses_of_mers.size(); i++)
                        starting_t_value += ( masses_of_mers.at(i) - 1 );
                    break;
        }

        //cout << "initial conditions " << initial_conditions << endl;

		// This loop performs a single simulation.
		for (int t = starting_t_value ; ; t++)
		{

		    //cout << "t = " << t  << "     labels_of_mers.size = " << labels_of_mers.size() << endl;
		    //getchar();

			if (t > t_max)
			{
				// save the histogram of the system (cluster size distribution)
				IVector cluster_distribution = histogram(masses_of_mers, N);
				container_cluster_distributions.push_back(cluster_distribution);

				// and break this single simulation
				break;
			}

			// if everything is aggregated to one single giant cluster, then stop the simulation
			if( labels_of_mers.size() <= 1 )
			{
			    cout << "Everything is aggregated into a single giant cluster!" << endl;
			    break;
			}

			int ind1;
			int ind2;

            // --------------------- Product Rule Fast ----------------------
			if(rule == "product_fast")
            {
                // choosing randomly two clusters
                do
                {
                    // I randomly choose clusters from the vector of available monomers

                    int which_monomer_1 = distribution(generator) % N;
                    int which_monomer_2 = distribution(generator) % N;

                    ind1 = membership_of_monomers.at(which_monomer_1);
                    ind2 = membership_of_monomers.at(which_monomer_2);

                    // If I picked up two different clusters its ok and go on, if not, repeat choosing.
                    if( ind1 != ind2 )
                        break;
                }
                while( true );
            }
            // ---------------------- Constant Rule Fast -----------------------------------
            else if(rule == "constant_fast")
            {
                // choosing randomly two clusters
                do
                {
                    // I randomly choose clusters from the vector of labels of existing clusters.

                    int which_mer_1 = distribution(generator) % labels_of_mers.size();
                    int which_mer_2 = distribution(generator) % labels_of_mers.size();

                    ind1 = labels_of_mers.at( which_mer_1 );
                    ind2 = labels_of_mers.at( which_mer_2 );

                    // If I picked up two different clusters its ok and go on, otherwise, repeat choosing.
                    if( ind1 != ind2 )
                        break;

                }
                while( true );
            }
            // ---------------------- Additive Rule Fast -----------------------------------
            else if(rule == "additive_fast")
            {
                // choosing randomly two clusters
                do
                {
                    // The first cluster is chosen from the labels of existing clusters.

                    int which_mer_1 = distribution(generator) % labels_of_mers.size();
                    ind1 = labels_of_mers.at( which_mer_1 );

                    // The second cluster is chosen from the vector of available monomers.

                    int which_monomer_2 = distribution(generator) % N;
                    ind2 = membership_of_monomers.at(which_monomer_2);

                    // If I picked up two different clusters its ok and go on, otherwise, repeat choosing.
                    if( ind1 != ind2 )
                        break;
                }
                while( true );
            }
            // ---------------------- Other rules - generalized method -----------------------------------
            else
            {
                // Here, I create the vector that contain all possible combinations of two clusters merging.
                // In the single field of this vector, there will be 3 numbers: cumulative propability K, mass i, and mass j.

                struct Field
                {
                    double value;
                    int i;
                    int j;
                } field;

                // vector of cumulative sums of probability
                vector<Field> vec_cum_sum_prob;

                double sum_probab = 0.;

                for(int i=0; i<masses_of_mers.size(); i++)
                    for(int j=i+1; j<masses_of_mers.size(); j++)
                    {
                        sum_probab += pKernel( masses_of_mers.at(i) , masses_of_mers.at(j) );

                        field.value = sum_probab;
                        field.i = i;
                        field.j = j;
                        vec_cum_sum_prob.push_back( field );
                    }

                // Now, we will randomly choose one configuration from the vector vec_cum_sum_prob which now holds
                // all possible coagulations acts that can occur in this moment.

                do
                {
                    // I take a random number from the range [1, sum_probab], which will tell us, which combination
                    // shall we take from the all available configurations
                    double random_number = sum_probab * double(distribution(generator))/maxInt;

                    int configuration = 0;

                    // Now I scan the vector of available configurations, starting from left to right.
                    // I scan it until k.value is bigger than my random number. If k.value is bigger, it means that
                    // the previous configuration from the vector was the configuration of my interest.
                    for(auto k : vec_cum_sum_prob)
                        if( k.value > random_number )
                            break;
                        else
                            configuration++;

                    ind1 = labels_of_mers.at( vec_cum_sum_prob.at(configuration).i );
                    ind2 = labels_of_mers.at( vec_cum_sum_prob.at(configuration).j );

                    // If I picked up two different clusters its ok and go on, otherwise, repeat choosing.
                    // To be honest, this condition should be always true, taking into account the method with
                    // the configuration vector.
                    if( ind1 != ind2 )
                        break;
                }
                while( true );
            }

            // Now, we have two randomly chosen mers (clusters) that are going to be merged.
            // In fact, ind1 and ind2 are the labels of mers.

            // Guarantee that ind2 will be always lower than ind1
            if ( ind2 > ind1 )
            {
                  int help = ind1;
                  ind1 = ind2;
                  ind2 = help;
            }

            int this_will_grow = ind1;
            int this_will_vanish = ind2;

            int ind_mer_1, ind_mer_2;

            // I have labels of mers, so here, I have to find the indexes of these mers in my memory vectors
            for (unsigned int e = 0; e<labels_of_mers.size(); e++)
                if (labels_of_mers.at(e) == this_will_grow)
                    ind_mer_1 = e;
                else if (labels_of_mers.at(e) == this_will_vanish)
                    ind_mer_2 = e;

            // In the vector of labels of mers, I delete the mer which is anihilated.
            labels_of_mers.erase(labels_of_mers.begin() + ind_mer_2);

            // I increase the mass of the mer which is growing
            // and I delete the mass of the mer which is anihilated.
            masses_of_mers.at(ind_mer_1) += masses_of_mers.at(ind_mer_2);
            masses_of_mers.erase(masses_of_mers.begin() + ind_mer_2);

            // Here, I change the membership of monomers. The monomers which were assigned to mer which vanished
            // will be now assigned to the mer which is grown.
            for (auto& a : membership_of_monomers)
                if (a == this_will_vanish)
                    a = this_will_grow;

		} // the end of single simulation loop


		// Here, I perform partial averaging of cluster size histograms, in case on simulations > 1e6.
		if(simulation > 0 and simulation % 1000000 == 0)
        {
            DVector average_partial_cluster_distribution(N);

            for(unsigned int i = 0; i<container_cluster_distributions.size(); i++)
                for (unsigned int j = 0; j<average_partial_cluster_distribution.size(); j++)
                    average_partial_cluster_distribution.at(j) += container_cluster_distributions.at(i).at(j);

            for(auto& s : average_partial_cluster_distribution)
                s /= 1000000;

            container_averaged_number_of_clusters.push_back(average_partial_cluster_distribution);

            // clearing the distributions that were averaged above
            container_cluster_distributions.clear();
        }

	} // end of loop performing many simulations


	// Here, I average all cluster size distributions. The distributions that were averaged ("partially")
	// during simulation, they have weight of 1000000, and the distributions which were not "partially" averaged,
	// they have weight of 1.

	DVector average_numbers_of_clusters(N);

	for(unsigned int i = 0; i<container_averaged_number_of_clusters.size(); i++)
        for (unsigned int j = 0; j<average_numbers_of_clusters.size(); j++)
            average_numbers_of_clusters.at(j) += 1000000*container_averaged_number_of_clusters.at(i).at(j);

    for(unsigned int i = 0; i<container_cluster_distributions.size(); i++)
        for (unsigned int j = 0; j<average_numbers_of_clusters.size(); j++)
            average_numbers_of_clusters.at(j) += container_cluster_distributions.at(i).at(j);

    for(auto& s : average_numbers_of_clusters)
        s /= number_of_simulations;

	cout << "Container of cluster size distributions is of length of: " << container_cluster_distributions.size() << endl;
	cout << "Single cluster size distribution is of lenght of: " << container_cluster_distributions.at(0).size() << endl;

	// Here, I calculate standard deviation of data.

	DVector stddev_numbers_of_clusters(N);

    for(unsigned int i = 0; i<container_cluster_distributions.size(); i++)
        for (unsigned int j = 0; j<stddev_numbers_of_clusters.size(); j++)
        {
            stddev_numbers_of_clusters.at(j) += ( container_cluster_distributions.at(i).at(j) - average_numbers_of_clusters.at(j) )
                                         * ( container_cluster_distributions.at(i).at(j) - average_numbers_of_clusters.at(j) );
        }

    for(auto& s : stddev_numbers_of_clusters)
        s = std::sqrt( s / ( number_of_simulations - 1 ) );


	// Saving results to files

	ofstream out(output_file_avns_string.c_str());
	out << setprecision(12);
	for (int i = 0; i < N; ++i)
        out << i + 1 << " " << average_numbers_of_clusters.at(i) << endl;
	out.close();

	cout << endl << "<ns> saved as " << endl
        << output_file_avns_string.c_str() << endl;

	ofstream out_stddev(output_file_stddev_string.c_str());
	out_stddev << setprecision(12);
	for (int i = 0; i < N; ++i)
        out_stddev << i + 1 << " " << stddev_numbers_of_clusters.at(i) << endl;
	out_stddev.close();

	cout << "Std. dev. saved as " << endl
        << output_file_stddev_string.c_str() << endl << endl;

	return 0;
}
