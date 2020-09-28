void PrintMatrixOfPresence(const boost::filesystem::path& path_output_folder,
	const unsigned int ndata,
	const Eigen::VectorXi& number_important_species,
	const Eigen::MatrixXi& important_species)
{
	boost::filesystem::path filename = path_output_folder / "presence.out";
	std::ofstream fOut(filename.c_str(), std::ios::out);
	for (int j = 0; j < ndata; j++)
	{
		fOut << std::setw(5) << std::left << number_important_species(j);
		for (int k = 0; k < important_species.cols(); k++)
			fOut << std::setw(2) << std::left << important_species(j, k);
		fOut << std::endl;
	}
	fOut.close();
}

void PrintErrorAnalysis(const boost::filesystem::path& path_output_folder,
	const unsigned int ndata, const unsigned int nclusters,
	const Eigen::VectorXi& number_important_species,
	const Eigen::MatrixXi& important_species,
	const std::vector< std::vector<int> > belonging,
	const std::vector<std::string> species_names,
	const double retained_threshold,
	Eigen::MatrixXi& retained_species)
{
	const unsigned int ns = important_species.cols();

	Eigen::MatrixXi sums(nclusters, ns); sums.setZero();
	Eigen::MatrixXd ratios(nclusters, ns); ratios.setZero();
	Eigen::MatrixXd errors(nclusters, ns); errors.setZero();
	retained_species.resize(nclusters, ns); retained_species.setZero();

	Eigen::VectorXi retained(nclusters); retained.setZero();
	Eigen::VectorXi retained_000(nclusters); retained_000.setZero();
	Eigen::VectorXi retained_001(nclusters); retained_001.setZero();
	Eigen::VectorXi retained_002(nclusters); retained_002.setZero();
	Eigen::VectorXi retained_005(nclusters); retained_005.setZero();
	Eigen::VectorXi retained_010(nclusters); retained_010.setZero();
	Eigen::VectorXi retained_020(nclusters); retained_020.setZero();

	// Analysis of errors
	for (int i = 0; i < nclusters; i++)
	{
		for (int j = 0; j < belonging[i].size(); j++)
		{
			const int g = belonging[i][j];
			for (int k = 0; k < ns; k++)
				sums(i, k) += important_species(g, k);
		}

		for (int k = 0; k < ns; k++)
			ratios(i, k) = double(sums(i, k)) / double(belonging[i].size());

		for (int k = 0; k < ns; k++)
			if (ratios(i, k) > 0.) errors(i, k) = std::fabs(1. - ratios(i, k));

		// Retain species only if the frequency of presence is larger than a threshold
		for (int k = 0; k < ns; k++)
		{
			if (ratios(i, k) > 0.00)	retained_000(i)++;
			if (ratios(i, k) > 0.01)	retained_001(i)++;
			if (ratios(i, k) > 0.02)	retained_002(i)++;
			if (ratios(i, k) > 0.05)	retained_005(i)++;
			if (ratios(i, k) > 0.10)	retained_010(i)++;
			if (ratios(i, k) > 0.20)	retained_020(i)++;

			if (ratios(i, k) > retained_threshold)
			{
				retained(i)++;
				retained_species(i, k) = 1;
			}
		}
	}

	{
		boost::filesystem::path filename = path_output_folder / "sums.out";
		std::ofstream fSums(filename.c_str(), std::ios::out);

		fSums << std::left << std::setw(7) << "#";
		fSums << std::left << std::setw(12) << "Items";

		unsigned int count = 3;
		for (int k = 0; k < ns; k++)
		{
			std::stringstream index; index << count++;
			std::string label = species_names[k] + "(" + index.str() + ")";
			fSums << std::left << std::setw(20) << label;
		}
		fSums << std::endl;

		// Local
		for (int i = 0; i < nclusters; i++)
		{
			fSums << std::left << std::setw(7) << i;
			fSums << std::left << std::setw(12) << belonging[i].size();
			for (int k = 0; k < ns; k++)
				fSums << std::left << std::setw(20) << sums(i, k);
			fSums << std::endl;
		}

		// Global
		{
			fSums << std::setw(7) << std::left << "Tot";
			fSums << std::setw(12) << std::left << ndata;
			for (int k = 0; k < ns; k++)
			{
				unsigned int sum = 0;
				for (int i = 0; i < nclusters; i++)
					sum += sums(i, k);
				fSums << std::setw(20) << std::left << sum;
			}
			fSums << std::endl;
		}

		fSums.close();
	}

	// Errors
	{
		Eigen::VectorXd sum_errors(nclusters); sum_errors.setZero();
		for (int i = 0; i < nclusters; i++)
		{
			for (int k = 0; k < ns; k++)
				sum_errors(i) += errors(i, k);
			sum_errors(i) *= 100. / double(ns);
		}

		{
			boost::filesystem::path filename = path_output_folder / "errors.out";
			std::ofstream fErrors(filename.c_str(), std::ios::out);

			fErrors << std::left << std::fixed << std::setw(5) << "#";
			fErrors << std::left << std::fixed << std::setw(12) << "error(%)";
			fErrors << std::left << std::fixed << std::setw(7) << "S(0%)";
			fErrors << std::left << std::fixed << std::setw(7) << "S(1%)";
			fErrors << std::left << std::fixed << std::setw(7) << "S(2%)";
			fErrors << std::left << std::fixed << std::setw(7) << "S(5%)";
			fErrors << std::left << std::fixed << std::setw(7) << "S(10%)";
			fErrors << std::left << std::fixed << std::setw(7) << "S(20%)";
			fErrors << std::endl;

			for (int i = 0; i < nclusters; i++)
			{
				fErrors << std::left << std::fixed << std::setw(5) << i;
				fErrors << std::left << std::fixed << std::setprecision(3) << std::setw(12) << sum_errors(i);
				fErrors << std::left << std::fixed << std::setw(7) << retained_000(i);
				fErrors << std::left << std::fixed << std::setw(7) << retained_001(i);
				fErrors << std::left << std::fixed << std::setw(7) << retained_002(i);
				fErrors << std::left << std::fixed << std::setw(7) << retained_005(i);
				fErrors << std::left << std::fixed << std::setw(7) << retained_010(i);
				fErrors << std::left << std::fixed << std::setw(7) << retained_020(i);
				fErrors << std::endl;
			}
			fErrors.close();
		}

		std::cout << "Errors(%)" << std::endl;
		for (int i = 0; i < nclusters; i++)
		{
			std::cout << std::left << std::fixed << std::setw(5) << i;
			std::cout << std::left << std::fixed << std::setprecision(3) << std::setw(12) << sum_errors(i);
			std::cout << std::left << std::fixed << std::setw(7) << retained_000(i);
		}
	}

	// Similarities between groups
	{
		std::cout << " * Writing similarity file..." << std::endl;

		Eigen::MatrixXd similarities(nclusters, nclusters);
		similarities.setZero();
		{
			for (int i = 0; i < nclusters; i++)
				for (int j = 0; j < nclusters; j++)
				{
					double sum = 0.;
					for (int k = 0; k < ns; k++)
						sum += std::fabs(retained_species(i, k) - retained_species(j, k));
					sum /= double(ns);
					similarities(i, j) = 1. - sum;
				}
		}

		boost::filesystem::path filename = path_output_folder / "similarities.out";
		std::ofstream fSimilarities(filename.c_str(), std::ios::out);

		fSimilarities << std::setw(12) << std::left << "Clusters";
		for (int i = 0; i < nclusters; i++)
			fSimilarities << std::setw(6) << std::left << i;
		fSimilarities << std::endl;

		for (int i = 0; i < nclusters; i++)
		{
			fSimilarities << std::left << std::fixed << std::setw(12) << i;
			for (int j = 0; j < nclusters; j++)
				fSimilarities << std::setw(6) << std::setprecision(3) << std::left << similarities(i, j);
			fSimilarities << std::endl;
		}

		fSimilarities.close();
	}

	// Uniformity coefficients
	{
		std::cout << " * Writing similarity file..." << std::endl;

		boost::filesystem::path filename = path_output_folder / "uniformity.out";
		std::ofstream fUniformity(filename.c_str(), std::ios::out);

		fUniformity << std::setw(10) << std::setprecision(3) << std::left << "Cluster";
		fUniformity << std::setw(10) << std::setprecision(3) << std::left << "Samples";
		fUniformity << std::setw(10) << std::setprecision(3) << std::left << "Species";
		fUniformity << std::setw(10) << std::setprecision(6) << std::left << "lambda";
		fUniformity << std::endl;

		Eigen::VectorXi nspecies(nclusters); nspecies.setZero();
		Eigen::VectorXd lambda(nclusters); lambda.setZero();

		for (int i = 0; i < nclusters; i++)
		{
			// Participation index
			Eigen::VectorXd x(ns); x.setZero();
			for (int j = 0; j < belonging[i].size(); j++)
			{
				const int g = belonging[i][j];
				for (int k = 0; k < ns; k++)
					x(k) += important_species(g, k);
			}

			if (belonging[i].size() != 0)
				x /= static_cast<double>(belonging[i].size());

			// Uniformity coefficient
			for (int k = 0; k < ns; k++)
			{
				if (x(k) != 0)
				{
					nspecies(i)++;
					lambda(i) += std::pow(x(k) - 1., 2.);
				}
			}

			lambda(i) = std::sqrt(lambda(i));
			if (nspecies(i) != 0)
				lambda(i) /= static_cast<double>(nspecies(i));

			fUniformity << std::setw(10) << std::setprecision(3) << std::left << i;
			fUniformity << std::setw(10) << std::setprecision(3) << std::left << belonging[i].size();
			fUniformity << std::setw(10) << std::setprecision(3) << std::left << nspecies(i);
			fUniformity << std::setw(10) << std::setprecision(6) << std::left << std::fixed << lambda(i);
			fUniformity << std::endl;
		}

		double mean_nspecies = 0.;
		double mean_lambda = 0.;
		unsigned int n_feasible = 0;
		unsigned int min_nspecies = ns;
		double min_lambda = 1.;
		double sigma2_nspecies = 0.;
		double sigma2_lambda = 0.;
		for (int i = 0; i < nclusters; i++)
		{
			if (nspecies(i) != 0)
			{
				mean_nspecies += nspecies(i);
				mean_lambda += lambda(i);

				if (nspecies(i) < min_nspecies)
					min_nspecies = nspecies(i);

				if (lambda(i) < min_lambda)
					min_lambda = lambda(i);

				sigma2_nspecies += nspecies(i) * nspecies(i);
				sigma2_lambda += lambda(i) * lambda(i);

				n_feasible++;
			}
		}
		mean_nspecies /= static_cast<double>(n_feasible);
		mean_lambda /= static_cast<double>(n_feasible);
		sigma2_nspecies -= n_feasible * mean_nspecies * mean_nspecies;
		sigma2_nspecies /= static_cast<double>(n_feasible);
		sigma2_lambda -= n_feasible * mean_lambda * mean_lambda;
		sigma2_lambda /= static_cast<double>(n_feasible);
		const double sigma_nspecies = std::sqrt(sigma2_nspecies);
		const double sigma_lambda = std::sqrt(sigma2_lambda);


		fUniformity << std::setw(10) << std::setprecision(3) << std::left << "Zeros";
		fUniformity << std::setw(10) << std::setprecision(3) << std::left << ndata;
		fUniformity << std::setw(10) << std::setprecision(3) << std::left << nclusters - n_feasible;
		fUniformity << std::setw(10) << std::setprecision(6) << std::left << std::fixed << (nclusters - n_feasible) / static_cast<double>(nclusters);
		fUniformity << std::endl;

		fUniformity << std::setw(10) << std::setprecision(3) << std::left << "Mean";
		fUniformity << std::setw(10) << std::setprecision(3) << std::left << ndata;
		fUniformity << std::setw(10) << std::setprecision(3) << std::left << mean_nspecies;
		fUniformity << std::setw(10) << std::setprecision(6) << std::left << std::fixed << mean_lambda;
		fUniformity << std::endl;

		fUniformity << std::setw(10) << std::setprecision(3) << std::left << "Sigma";
		fUniformity << std::setw(10) << std::setprecision(3) << std::left << ndata;
		fUniformity << std::setw(10) << std::setprecision(3) << std::left << sigma_nspecies;
		fUniformity << std::setw(10) << std::setprecision(6) << std::left << std::fixed << sigma_lambda;
		fUniformity << std::endl;

		fUniformity << std::setw(10) << std::setprecision(3) << std::left << "Min";
		fUniformity << std::setw(10) << std::setprecision(3) << std::left << ndata;
		fUniformity << std::setw(10) << std::setprecision(3) << std::left << min_nspecies;
		fUniformity << std::setw(10) << std::setprecision(6) << std::left << std::fixed << min_lambda;
		fUniformity << std::endl;

		fUniformity << std::setw(10) << std::setprecision(3) << std::left << "Max";
		fUniformity << std::setw(10) << std::setprecision(3) << std::left << ndata;
		fUniformity << std::setw(10) << std::setprecision(3) << std::left << nspecies.maxCoeff();
		fUniformity << std::setw(10) << std::setprecision(6) << std::left << std::fixed << lambda.maxCoeff();
		fUniformity << std::endl;

		fUniformity.close();
	}
}

void SelectImportantReactions(OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
	Eigen::MatrixXi& retained_species,
	Eigen::MatrixXi& retained_reactions)
{
	std::cout << " * Selecting important reactions..." << std::endl;

	unsigned int nclusters = retained_species.rows();
	unsigned int ns = retained_species.cols();
	unsigned int nr = kineticsMap.NumberOfReactions();
	Eigen::SparseMatrix<int> delta_sparse(ns, nr);

	{
		// Build a full matrix of net stoichiometric coefficients nu = nuB - nuF
		Eigen::MatrixXi nu(nr, ns);
		{
			// Be careful: eigen vectors and matrices are 0-index based
			// Be careful: if the kinetic scheme is large, this matrix, since it is full, can be very memory expensive
			//             Example: 10^3 species, 10^4 reactions = size of the matrix 10^7 elements!
			//             This is the reason why we store stoichiometric matrices in sparse format.
			//             Of course te advantage of having a full matrix, is that you access the elements directly, without
			//             using iterators and pointers, as reported above
			nu.setZero();

			// Loop over all the reactions (product side)
			for (int k = 0; k < kineticsMap.stoichiometry().stoichiometric_matrix_products().outerSize(); ++k)
			{
				// Loop over all the non-zero stoichiometric coefficients (product side) of reaction k
				for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMap.stoichiometry().stoichiometric_matrix_products(), k); it; ++it)
				{
					nu(it.row(), it.col()) = 1;
				}
			}

			// Loop over all the reactions (product side)
			for (int k = 0; k < kineticsMap.stoichiometry().stoichiometric_matrix_reactants().outerSize(); ++k)
			{
				// Loop over all the non-zero stoichiometric coefficients (product side) of reaction k
				for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMap.stoichiometry().stoichiometric_matrix_reactants(), k); it; ++it)
				{
					nu(it.row(), it.col()) = 1;
				}
			}
		}

		// Sparse delta matrix, 1 means the species is involved in the reaction, (NR x NS)
		{
			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList;
			tripletList.reserve(nr * 4);
			for (unsigned int i = 0; i < nr; i++)
				for (unsigned int j = 0; j < ns; j++)
				{
					if (nu(i, j) != 0)
						tripletList.push_back(T(j, i, 1));
				}

			delta_sparse.setFromTriplets(tripletList.begin(), tripletList.end());
		}
	}

	// Important reactions
	retained_reactions.resize(nclusters, kineticsMap.NumberOfReactions());
	retained_reactions.setConstant(1);
	for (unsigned int i = 0; i < nclusters; i++)
	{
		for (int k = 0; k < delta_sparse.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<int>::InnerIterator it(delta_sparse, k); it; ++it)
			{
				if (retained_species(i, it.row()) == 0)
				{
					retained_reactions(i, k) = 0;
				}
			}
		}
	}

	// Number of retained reactions
	for (unsigned int i = 0; i < nclusters; i++)
	{
		const unsigned int number_retained_reactions = retained_reactions.row(i).sum();
		std::cout << i << " " << number_retained_reactions << std::endl;
	}

	// Include (if needed) third-body species
	{
		std::vector<unsigned int> third_body_reactions = kineticsMap.IndicesOfThirdbodyReactions();
		const std::vector< std::vector<unsigned int> > third_body_species = kineticsMap.IndicesOfThirdbodySpecies();

		for (unsigned int i = 0; i < nclusters; i++)
		{
			// First loop
			for (unsigned int k = 0; k < third_body_reactions.size(); k++)
			{
				unsigned int index_reaction = third_body_reactions[k] - 1;
				if (retained_reactions(i, index_reaction) == 1)
				{
					std::cout << "Looking for third body reaction: " << index_reaction << std::endl;
					for (unsigned int j = 0; j < third_body_species[k].size(); j++)
					{
						unsigned int index_species = third_body_species[k][j] - 1;
						if (retained_species(i, index_species) == 0)
						{
							std::cout << "Adding species: " << index_species << std::endl;
							retained_species(i, index_species) = 1;
						}
					}
				}
			}

			// Falloff reactions
			for (unsigned int k = 0; k < kineticsMap.NumberOfFallOffReactions(); k++)
			{
				const unsigned int index_reaction = kineticsMap.IndicesOfFalloffReactions()[k] - 1;

				if (retained_reactions(i, index_reaction) == 1)
				{
					std::cout << "Looking for fall-off reaction: " << index_reaction << std::endl;

					if (kineticsMap.FallOffIndexOfSingleThirdbodySpecies()[k] == 0)
					{
						for (unsigned int j = 0; j < kineticsMap.FallOffIndicesOfThirdbodySpecies()[k].size(); j++)
						{
							unsigned int index_species = kineticsMap.FallOffIndicesOfThirdbodySpecies()[k][j] - 1;
							if (retained_species(i, index_species) == 0)
							{
								std::cout << "Adding species: " << index_species << std::endl;
								retained_species(i, index_species) = 1;
							}
						}
					}
				}
			}
		}
	}
}

void WriteKineticMechanisms(const OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
	const OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
	const boost::filesystem::path& path_output_folder,
	const boost::filesystem::path& chemkin_thermodynamic_file,
	const boost::filesystem::path& chemkin_kinetics_file,
	const Eigen::MatrixXi& retained_species,
	const Eigen::MatrixXi& retained_reactions)
{
	std::cout << " * Writing mechanisms..." << std::endl;

	boost::filesystem::path path_folder = path_output_folder / "Mechanisms";
	if (!boost::filesystem::exists(path_folder))
		OpenSMOKE::CreateDirectory(path_folder);

	// Log file
	std::ofstream flog;
	{
		boost::filesystem::path file_name = path_folder / "log";
		flog.open(file_name.c_str(), std::ios::out);
		flog.setf(std::ios::scientific);
	}

	//Preprocessing kinetic file
	typedef OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > > ThermoReader_CHEMKIN;
	ThermoReader_CHEMKIN* thermoreader;

	// Reading thermodynamic database
	thermoreader = new OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > >;
	thermoreader->ReadFromFile(chemkin_thermodynamic_file.string());

	//Preprocessing Chemical Kinetics
	typedef OpenSMOKE::PreProcessorKinetics< OpenSMOKE::PreProcessorKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_CHEMKIN> > PreProcessorKinetics_CHEMKIN;
	PreProcessorKinetics_CHEMKIN* preprocessor_kinetics;
	preprocessor_kinetics = new PreProcessorKinetics_CHEMKIN(flog);
	preprocessor_kinetics->ReadFromASCIIFile(chemkin_kinetics_file.string());

	// Preprocessing the thermodynamics
	typedef OpenSMOKE::Species< OpenSMOKE::ThermoPolicy_CHEMKIN, OpenSMOKE::TransportPolicy_CHEMKIN > SpeciesCHEMKIN;
	typedef OpenSMOKE::PreProcessorSpecies< OpenSMOKE::PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<SpeciesCHEMKIN> > PreProcessorSpecies_CHEMKIN_WithoutTransport;
	PreProcessorSpecies_CHEMKIN_WithoutTransport* preprocessor_species_without_transport;
	preprocessor_species_without_transport = new PreProcessorSpecies_CHEMKIN_WithoutTransport(*thermoreader, *preprocessor_kinetics, flog);
	CheckForFatalError(preprocessor_species_without_transport->Setup());

	// Read kinetics from file
	CheckForFatalError(preprocessor_kinetics->ReadKineticsFromASCIIFile(preprocessor_species_without_transport->AtomicTable()));

	delete thermoreader;
	delete preprocessor_species_without_transport;

	std::cout << " * Writing mechanisms..." << std::endl;
	unsigned int nclusters = retained_species.rows();
	unsigned int ns = retained_species.cols();
	for (unsigned int i = 0; i < nclusters; i++)
	{
		std::vector<bool> is_reduced_species(ns);
		for (unsigned int k = 0; k < ns; k++)
			is_reduced_species[k] = true;

		std::stringstream label; label << i;
		std::string name = "kinetics." + label.str() + ".CKI";
		boost::filesystem::path filename = path_folder / name;

		std::ofstream fKinetics(filename.c_str(), std::ios::out);

		fKinetics << "ELEMENTS" << std::endl;
		for (unsigned int k = 0; k < thermodynamicsMap.elements().size(); k++)
			fKinetics << thermodynamicsMap.elements()[k] << std::endl;
		fKinetics << "END" << std::endl;
		fKinetics << std::endl;

		fKinetics << "SPECIES" << std::endl;
		unsigned int count = 0;
		for (unsigned int k = 0; k < ns; k++)
		{
			if (retained_species(i, k) == 1)
			{
				count++;
				fKinetics << thermodynamicsMap.NamesOfSpecies()[k] << "  ";
				if (count % 6 == 0)	fKinetics << std::endl;
			}
		}
		fKinetics << std::endl;
		fKinetics << "END" << std::endl;
		fKinetics << std::endl;

		fKinetics << "REACTIONS" << std::endl;

		for (unsigned int k = 0; k < preprocessor_kinetics->reactions().size(); k++)
		{
			if (retained_reactions(i, k) == 1)
			{
				std::stringstream reaction_data;
				std::string reaction_string;
				preprocessor_kinetics->reactions()[k].GetReactionStringCHEMKIN(thermodynamicsMap.NamesOfSpecies(), reaction_data, is_reduced_species);
				fKinetics << reaction_data.str();
			}
		}
		fKinetics << "END" << std::endl;
		fKinetics << std::endl;

		fKinetics.close();
	}
}

void PreprocessKineticMechanisms(
	const unsigned int nclusters,
	const boost::filesystem::path& path_output_folder,
	const boost::filesystem::path& chemkin_thermodynamic_file,
	const boost::filesystem::path& chemkin_transport_file)
{
	typedef OpenSMOKE::Species< OpenSMOKE::ThermoPolicy_CHEMKIN, OpenSMOKE::TransportPolicy_CHEMKIN > SpeciesCHEMKIN;
	typedef OpenSMOKE::PreProcessorSpecies< OpenSMOKE::PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<SpeciesCHEMKIN>  > PreProcessorSpecies_CHEMKIN;
	typedef OpenSMOKE::PreProcessorKinetics< OpenSMOKE::PreProcessorKineticsPolicy_CHEMKIN<OpenSMOKE::ReactionPolicy_CHEMKIN> > PreProcessorKinetics_CHEMKIN;
	typedef OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > > ThermoReader_CHEMKIN;

	std::cout << " * Writing mechanisms..." << std::endl;

	boost::filesystem::path input_path_folder = path_output_folder / "Mechanisms";
	boost::filesystem::path output_path_folder = path_output_folder / "PreprocessedMechanisms";
	if (!boost::filesystem::exists(output_path_folder))
		OpenSMOKE::CreateDirectory(output_path_folder);

	for (unsigned int i = 0; i < nclusters; i++)
	{
		std::cout << "Preprocessing mechanism " << i + 1 << " over " << nclusters << std::endl;

		// Log file
		boost::filesystem::path file_name = output_path_folder / "log";
		std::ofstream flog;
		flog.open(file_name.c_str(), std::ios::out);
		flog.setf(std::ios::scientific);

		std::cout << " * Creating thermo reader..." << std::endl;

		// Readers
		ThermoReader_CHEMKIN* thermoreader;
		OpenSMOKE::TransportReader< OpenSMOKE::TransportReaderPolicy_CHEMKIN<OpenSMOKE::TransportPolicy_CHEMKIN > >* transportreader;

		// Reading thermodynamic database
		thermoreader = new OpenSMOKE::ThermoReader< OpenSMOKE::ThermoReaderPolicy_CHEMKIN< OpenSMOKE::ThermoPolicy_CHEMKIN > >;
		CheckForFatalError(thermoreader->ReadFromFile(chemkin_thermodynamic_file.string()));

		std::cout << " * Creating kinetics reader..." << std::endl;

		// Preprocessing the kinetic mechanism
		std::stringstream label; label << i;
		const std::string local_file_name = "kinetics." + label.str() + ".CKI";
		const boost::filesystem::path chemkin_kinetics_file = input_path_folder / local_file_name;
		PreProcessorKinetics_CHEMKIN preprocessor_kinetics(flog);
		CheckForFatalError(preprocessor_kinetics.ReadFromASCIIFile(chemkin_kinetics_file.string()));

		std::cout << " * Creating transport reader..." << std::endl;

		// Preprocessors
		PreProcessorSpecies_CHEMKIN* preprocessor_species_with_transport = nullptr;

		// Reading transport database
		transportreader = new OpenSMOKE::TransportReader< OpenSMOKE::TransportReaderPolicy_CHEMKIN<OpenSMOKE::TransportPolicy_CHEMKIN > >;
		CheckForFatalError(transportreader->ReadFromFile(chemkin_transport_file.string()));

		// Preprocessing the thermodynamic and transport properties
		preprocessor_species_with_transport = new PreProcessorSpecies_CHEMKIN(*thermoreader, *transportreader, preprocessor_kinetics, flog);

		// Preprocessing thermodynamics
		CheckForFatalError(preprocessor_species_with_transport->Setup());

		// Read kinetics from file
		CheckForFatalError(preprocessor_kinetics.ReadKineticsFromASCIIFile(preprocessor_species_with_transport->AtomicTable()));

		// Fit transport data
		CheckForFatalError(preprocessor_species_with_transport->Fitting());

		// Write XML files
		std::cout << " * Writing XML files..." << std::endl;
		{
			std::string author_name("undefined");
			std::string place_name("undefined");
			std::string comments("no comments");
			std::string preprocessing_date(OpenSMOKE::GetCurrentDate());
			std::string preprocessing_time(OpenSMOKE::GetCurrentTime());

			std::stringstream xml_string;
			xml_string << std::setprecision(8);
			xml_string.setf(std::ios::scientific);

			xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
			xml_string << "<opensmoke version=\"0.1a\">" << std::endl;

			xml_string << "<Properties>" << std::endl;
			xml_string << "  <Author>" << author_name << "</Author>" << std::endl;
			xml_string << "  <Place>" << place_name << "</Place>" << std::endl;
			xml_string << "  <Date>" << preprocessing_date << "</Date>" << std::endl;
			xml_string << "  <Time>" << preprocessing_time << "</Time>" << std::endl;
			xml_string << "  <Comments>" << "\n" << OpenSMOKE::SplitStringIntoSeveralLines(comments, 80, "\t\r ") << "\n" << "</Comments>" << std::endl;
			xml_string << "</Properties>" << std::endl;

			// Thermodynamics and transport properties
			preprocessor_species_with_transport->WriteXMLFile(xml_string);

			// Kinetic mechanism
			preprocessor_kinetics.WriteXMLFile(xml_string);

			xml_string << "</opensmoke>" << std::endl;

			// Write file
			std::stringstream label; label << i;
			const std::string folder_name = "kinetics." + label.str();
			const boost::filesystem::path local_folder_path = output_path_folder / folder_name;
			if (!boost::filesystem::exists(local_folder_path))
				OpenSMOKE::CreateDirectory(local_folder_path);

			boost::filesystem::path kinetics_xml = local_folder_path / "kinetics.xml";
			std::ofstream fOutput;
			fOutput.open(std::string(kinetics_xml.string()).c_str(), std::ios::out);
			fOutput.setf(std::ios::scientific);
			fOutput << xml_string.str();
			fOutput.close();
		}

		flog.close();
		delete transportreader;
		delete preprocessor_species_with_transport;
		delete thermoreader;
	}
}

double EpsilonDRG(const double T)
{
	const double w = 300.;
	const double epsMax = 0.05;
	const double epsMin = 0.01;
	const double Tm = 1000.;

	return epsMax - (epsMax - epsMin) * 0.50 * (1. + std::tanh((T - Tm) / w));
}

double EpsilonDRGEP(const double T)
{
	const double w = 300.;
	const double epsMax = 0.020;
	const double epsMin = 0.005;
	const double Tm = 1000.;

	return epsMax - (epsMax - epsMin) * 0.50 * (1. + std::tanh((T - Tm) / w));
}
