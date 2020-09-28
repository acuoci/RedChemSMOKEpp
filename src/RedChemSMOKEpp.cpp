//#include <iostream>
//#include <fstream>
//#include <Eigen/Dense>
//#include "neural/myNeuralNetworkBasedOnPCA.h"
//#include "fitctree/classifyPoint.h"
//#include "fitctree/classifyPoint_terminate.h"
//#include "fitctree/classifyPoint_initialize.h"

// Boost
//#include "boost / program_options.hpp

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// Thermodynamics
#include "kernel/thermo/Species.h"
#include "kernel/thermo/ThermoPolicy_CHEMKIN.h"
#include "kernel/thermo/ThermoReader.h"
#include "kernel/thermo/ThermoReaderPolicy_CHEMKIN.h"

// Kinetics
#include "kernel/kinetics/ReactionPolicy_CHEMKIN.h"

// Preprocessing
#include "preprocessing/PreProcessorSpecies.h"
#include "preprocessing/PreProcessorKinetics.h"
#include "preprocessing/PreProcessorKineticsPolicy_CHEMKIN.h"
#include "preprocessing/PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport.h"
#include "preprocessing/PreProcessorSpeciesPolicy_CHEMKIN_WithTransport.h"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"
#include "DRG.h"
#include "DRGEP.h"

// Grammar
#include "Grammar_RedChemSMOKEpp.h"

// Utilities
#include "Utilities.h"


int main(int argc, char** argv)
{
	boost::filesystem::path executable_file = OpenSMOKE::GetExecutableFileName(argv);
	boost::filesystem::path executable_folder = executable_file.parent_path();

	OpenSMOKE::OpenSMOKE_logo("RedChemSMOKE++", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	//unsigned int max_number_allowed_species = 100000;
	//OpenSMOKE::OpenSMOKE_CheckLicense(executable_folder, "RedChemSMOKE++", max_number_allowed_species);

	std::string input_file_name_ = "input.dic";
	std::string main_dictionary_name_ = "RedChemSMOKE++";

	// Program options from command line
	{
		namespace po = boost::program_options;
		po::options_description description("Options for the RedChemSMOKE++");
		description.add_options()
			("help", "print help messages")
			("input", po::value<std::string>(), "name of the file containing the main dictionary (default \"input.dic\")")
			("dictionary", po::value<std::string>(), "name of the main dictionary to be used (default \"RedChemSMOKE++\")");

		po::variables_map vm;
		try
		{
			po::store(po::parse_command_line(argc, argv, description), vm); // can throw 

			if (vm.count("help"))
			{
				std::cout << "Basic Command Line Parameters" << std::endl;
				std::cout << description << std::endl;
				return OPENSMOKE_SUCCESSFULL_EXIT;
			}

			if (vm.count("input"))
				input_file_name_ = vm["input"].as<std::string>();

			if (vm.count("dictionary"))
				main_dictionary_name_ = vm["dictionary"].as<std::string>();

			po::notify(vm); // throws on error, so do after help in case  there are any problems 
		}
		catch (po::error& e)
		{
			std::cerr << "Fatal error: " << e.what() << std::endl << std::endl;
			std::cerr << description << std::endl;
			return OPENSMOKE_FATAL_ERROR_EXIT;
		}
	}

	// Defines the grammar rules
	OpenSMOKE::Grammar_RedChemSMOKEpp grammar_chemredsmokepp;

	// Define the dictionaries
	OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
	dictionaries.ReadDictionariesFromFile(input_file_name_);
	dictionaries(main_dictionary_name_).SetGrammar(grammar_chemredsmokepp);

	// Kinetic folder
	boost::filesystem::path path_kinetics_folder;
	if (dictionaries(main_dictionary_name_).CheckOption("@KineticsFolder") == true)
	{
		dictionaries(main_dictionary_name_).ReadPath("@KineticsFolder", path_kinetics_folder);
		OpenSMOKE::CheckKineticsFolder(path_kinetics_folder);
	}

	// Kinetic file
	boost::filesystem::path chemkin_kinetics_file;
	if (dictionaries(main_dictionary_name_).CheckOption("@Kinetics") == true)
		dictionaries(main_dictionary_name_).ReadPath("@Kinetics", chemkin_kinetics_file);

	// Thermodynamic file
	boost::filesystem::path chemkin_thermodynamics_file;
	if (dictionaries(main_dictionary_name_).CheckOption("@Thermodynamics") == true)
		dictionaries(main_dictionary_name_).ReadPath("@Thermodynamics", chemkin_thermodynamics_file);

	// Thermodynamic file
	bool iTransport = false;
	boost::filesystem::path chemkin_transport_file;
	if (dictionaries(main_dictionary_name_).CheckOption("@Transport") == true)
	{
		iTransport = true;
		dictionaries(main_dictionary_name_).ReadPath("@Transport", chemkin_transport_file);
	}

	// Input file
	boost::filesystem::path path_xml_input_file;
	if (dictionaries(main_dictionary_name_).CheckOption("@XMLInput") == true)
		dictionaries(main_dictionary_name_).ReadPath("@XMLInput", path_xml_input_file);

	// Output folder
	boost::filesystem::path path_output_folder;
	if (dictionaries(main_dictionary_name_).CheckOption("@Output") == true)
	{
		dictionaries(main_dictionary_name_).ReadPath("@Output", path_output_folder);

		if (!boost::filesystem::exists(path_output_folder))
			OpenSMOKE::CreateDirectory(path_output_folder);
	}

	// DRG Analysis
	bool iDRG = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@DRG") == true)
		dictionaries(main_dictionary_name_).ReadBool("@DRG", iDRG);

	// DRG-EP Analysis
	bool iDRGEP = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@DRGEP") == true)
		dictionaries(main_dictionary_name_).ReadBool("@DRGEP", iDRGEP);

	// Testing neural network
	bool iTestingNeuralNetwork = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@TestingNeuralNetwork") == true)
		dictionaries(main_dictionary_name_).ReadBool("@TestingNeuralNetwork", iTestingNeuralNetwork);

	// Testing tree
	bool iTestingFitCTree = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@TestingFitCTree") == true)
		dictionaries(main_dictionary_name_).ReadBool("@TestingFitCTree", iTestingFitCTree);

	// Threshold
	double epsilon = 0.01;
	if (dictionaries(main_dictionary_name_).CheckOption("@Epsilon") == true)
		dictionaries(main_dictionary_name_).ReadDouble("@Epsilon", epsilon);

	// Temperature threshold
	double T_threshold = 300;
	if (dictionaries(main_dictionary_name_).CheckOption("@TemperatureThreshold") == true)
		dictionaries(main_dictionary_name_).ReadDouble("@TemperatureThreshold", T_threshold);

	// Concentration threshold
	double c_threshold = 1.e-12;
	if (dictionaries(main_dictionary_name_).CheckOption("@ConcentrationThreshold") == true)
		dictionaries(main_dictionary_name_).ReadDouble("@ConcentrationThreshold", c_threshold);

	// Retained threshold
	double retained_threshold = 0.01;
	if (dictionaries(main_dictionary_name_).CheckOption("@RetainedThreshold") == true)
		dictionaries(main_dictionary_name_).ReadDouble("@RetainedThreshold", retained_threshold);

	// Pressure 
	double P = 101325.;
	if (dictionaries(main_dictionary_name_).CheckOption("@Pressure") == true)
	{
		std::string units;
		dictionaries(main_dictionary_name_).ReadMeasure("@Pressure", P, units);
		if (units == "Pa")			P *= 1.;
		else if (units == "atm")	P *= 101325.;
		else if (units == "bar")	P *= 100000.;
		else OpenSMOKE::FatalErrorMessage("Wrong Pressure units");
	}

	// List of key species
	std::vector<std::string> key_species;
	if (dictionaries(main_dictionary_name_).CheckOption("@KeySpecies") == true)
		dictionaries(main_dictionary_name_).ReadOption("@KeySpecies", key_species);

	// Diffusion map analysis
	bool iDiffusionMapAnalysis = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@DiffusionMapAnalysis") == true)
		dictionaries(main_dictionary_name_).ReadBool("@DiffusionMapAnalysis", iDiffusionMapAnalysis);


	// Applying Diffusion Map Analysis
	if (iDiffusionMapAnalysis == true)
	{
		std::cout << "Applying Diffusion Map Analysis..." << std::endl;

		OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMap;
		OpenSMOKE::KineticsMap_CHEMKIN* kineticsMap;

		// Reading thermodynamic and kinetic files
		{
			boost::property_tree::ptree ptree;
			boost::property_tree::read_xml((path_kinetics_folder / "kinetics.xml").string(), ptree);

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

			thermodynamicsMap = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
			kineticsMap = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMap, ptree);

			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			std::cout << " * Time to read XML file: " << tEnd - tStart << std::endl;
		}

		// Diffusion maps: weight matrix
		{
			std::cout << "Diffusion map: weight matrix" << std::endl;

			const int NS = thermodynamicsMap->NumberOfSpecies();
			const int NR = kineticsMap->NumberOfReactions();

			// Diffusion Map Temperature 
			double T = 101325.;
			if (dictionaries(main_dictionary_name_).CheckOption("@DiffusionMapTemperature") == true)
			{
				std::string units;
				dictionaries(main_dictionary_name_).ReadMeasure("@DiffusionMapTemperature", T, units);
				if (units == "K")			T = T;
				else if (units == "C")		T = T + 273.15;
				else OpenSMOKE::FatalErrorMessage("Wrong Temperature units");
			}

			// Diffusion Map Pressure 
			double P_Pa = 101325.;
			if (dictionaries(main_dictionary_name_).CheckOption("@DiffusionMapPressure") == true)
			{
				std::string units;
				dictionaries(main_dictionary_name_).ReadMeasure("@DiffusionMapPressure", P_Pa, units);
				if (units == "Pa")			P_Pa *= 1.;
				else if (units == "atm")	P_Pa *= 101325.;
				else if (units == "bar")	P_Pa *= 100000.;
				else OpenSMOKE::FatalErrorMessage("Wrong Pressure units");
			}

			OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMap->NumberOfSpecies());
			x = 1. / static_cast<double>(NS);

			// Input file
			boost::filesystem::path path_input_file;
			if (dictionaries(main_dictionary_name_).CheckOption("@DiffusionMapMassFractions") == true)
			{
				dictionaries(main_dictionary_name_).ReadPath("@DiffusionMapMassFractions", path_input_file);

				OpenSMOKE::OpenSMOKEVectorDouble y(thermodynamicsMap->NumberOfSpecies());
				std::ifstream fIn(path_input_file.c_str(), std::ios::in);
				for (int i = 0; i < NS; i++)
					fIn >> y[i + 1];
				fIn.close();

				double MW;
				thermodynamicsMap->MoleFractions_From_MassFractions(x.GetHandle(), MW, y.GetHandle());
			}

			const double cTot = P_Pa / 8314 / T;
			OpenSMOKE::OpenSMOKEVectorDouble c(thermodynamicsMap->NumberOfSpecies()); c = x; c *= cTot;
			OpenSMOKE::OpenSMOKEVectorDouble R(kineticsMap->NumberOfReactions());

			// Now we know T, P and composition.
			// We have to pass those data to the thermodynamic and kinetic maps
			kineticsMap->SetTemperature(T);
			kineticsMap->SetPressure(P_Pa);
			thermodynamicsMap->SetTemperature(T);
			thermodynamicsMap->SetPressure(P_Pa);

			// Now we can calculate (internally) the reaction rates concentrations are needed
			kineticsMap->ReactionRates(c.GetHandle());

			// Build a full matrix of net stoichiometric coefficients nu = nuB - nuF
			std::cout << "Build stoichiometrix matrix" << std::endl;
			Eigen::MatrixXd nu(NR, NS);
			{
				nu.setZero();

				// Loop over all the reactions (product side)
				for (int k = 0; k < kineticsMap->stoichiometry().stoichiometric_matrix_products().outerSize(); ++k)
				{
					// Loop over all the non-zero stoichiometric coefficients (product side) of reaction k
					for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMap->stoichiometry().stoichiometric_matrix_products(), k); it; ++it)
					{
						nu(it.row(), it.col()) += it.value();
					}
				}

				// Loop over all the reactions (product side)
				for (int k = 0; k < kineticsMap->stoichiometry().stoichiometric_matrix_reactants().outerSize(); ++k)
				{
					// Loop over all the non-zero stoichiometric coefficients (product side) of reaction k
					for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMap->stoichiometry().stoichiometric_matrix_reactants(), k); it; ++it)
					{
						nu(it.row(), it.col()) -= it.value();
					}
				}
			}

			std::vector<double> kappa = kineticsMap->KArrheniusModified();

			std::cout << "Build weight matrix" << std::endl;
			Eigen::MatrixXd W(NS, NS);
			W.setZero();
			for (int i = 0; i < NS; i++)
				for (int j = 0; j < NS; j++)
					if (i != j)
					{
						Eigen::VectorXd ri = nu.col(i);
						Eigen::VectorXd rj = nu.col(j);

						double maxk = 0.;
						for (int k = 0; k < NR; k++)
							if (ri(k) != 0 && rj(k) != 0)
								if (kappa[k] > maxk) maxk = kappa[k];
						W(i, j) = maxk;
					}
			for (int i = 0; i < NS; i++)
				W(i, i) = W.row(i).maxCoeff();

			const double epsilon = 0;
			for (int i = 0; i < NS; i++)
				for (int j = 0; j < NS; j++)
					if (W(i, j) == 0.)
						W(i, j) = epsilon;
			{
				std::cout << "Write stoichiometrix matrix" << std::endl;
				std::ofstream fOut("W.out", std::ios::out);
				fOut.setf(std::ios::scientific);
				for (int i = 0; i < NS; i++)
				{
					for (int j = 0; j < NS; j++)
						fOut << W(i, j) << " ";
					fOut << std::endl;
				}
				fOut.close();
			}

			{
				std::cout << "Write labels" << std::endl;
				std::ofstream fOut("labels.out", std::ios::out);
				fOut.setf(std::ios::scientific);
				fOut << "labels = { ";
				for (int i = 0; i < NS; i++)
					fOut << "'" << thermodynamicsMap->NamesOfSpecies()[i] << "', ";
				fOut << " };" << std::endl;
				fOut.close();
			}
		}

		std::cout << "Diffusion map analysis completed..." << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(0);
	}


	// Applying Static Reduction
	{
		std::cout << "Applying static reduction..." << std::endl;

		OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMap;
		OpenSMOKE::KineticsMap_CHEMKIN* 		kineticsMap;

		// Reading thermodynamic and kinetic files
		{ 
			boost::property_tree::ptree ptree;
			boost::property_tree::read_xml((path_kinetics_folder / "kinetics.xml").string(), ptree);

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

			thermodynamicsMap = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
			kineticsMap = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMap, ptree);

			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			std::cout << " * Time to read XML file: " << tEnd - tStart << std::endl;
		}

		std::cout << "Reading input data..." << std::endl;

		bool read_pca_data = false;
		unsigned int nclusters = 0;					// total number of clusters
		unsigned int ndata = 0;						// total number of items
		unsigned int ns = 0;						// number of original species
		unsigned int nf = 0;						// number of filtered species
		unsigned int npca = 0;						// number of principal components
		unsigned int nretspecies = 0;				// number of retained species
		Eigen::VectorXi group;						// group to which each item belongs
		Eigen::VectorXd csi;						// mixture fraction
		Eigen::VectorXd T;							// temperature [K]
		Eigen::MatrixXd omega;						// mass fractions (original)
		std::vector< std::vector<int> > belonging;	// list of items belonging to each group
		Eigen::VectorXd mu;							// means of filtered variables
		Eigen::VectorXd sigma;						// std deviations of filtered variables
		Eigen::MatrixXd w;							// PCA weights
		Eigen::MatrixXd pca;						// PCA data
		std::vector<int> listretspecies;			// list of retained species (1-based)

		{
			std::cout << "Opening XML file..." << std::endl;

			// Open the XML file
			boost::property_tree::ptree ptree;
			boost::property_tree::read_xml(path_xml_input_file.string(), ptree);
		
			// Read number of clusters
			std::cout << "Reading number of clusters [classes]" << std::endl;
			nclusters = ptree.get<unsigned int>("opensmoke.classes");
			std::cout << " * Number of clusters: " << nclusters << std::endl;

			// Read number of items
			std::cout << "Reading number of items [items]" << std::endl;
			ndata = ptree.get<unsigned int>("opensmoke.items");
			std::cout << " * Number of items: " << ndata << std::endl;
			
			// Read number of original species
			std::cout << "Reading number of original species [original-components]" << std::endl;
			ns = ptree.get<unsigned int>("opensmoke.original-components");
			ns -= 1; // one of the original components is always the temperature
			if (ns != thermodynamicsMap->NumberOfSpecies())
				OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "It seems that the preprocessed file is not consistent with the kinetic mechanism you are using, since the number of species is not the same");
			std::cout << " * Number of original species: " << ns << std::endl;

			// Read number of filtered species
			std::cout << "Reading number of filtered species [filtered-components]" << std::endl;
			nf = ptree.get<unsigned int>("opensmoke.filtered-components");
			nf -= 1; // one of the filtered components is always the temperature
			std::cout << " * Number of filtered species: " << nf << std::endl;
			std::cout << " * Number of removed species:  " << ns - nf << std::endl;

			if (read_pca_data == true)
			{
				// Read number of principal components
				std::cout << "Reading number of principal components [principal-components]" << std::endl;
				npca = ptree.get<unsigned int>("opensmoke.principal-components");
				std::cout << " * Number of principal components: " << npca << std::endl;
			}

			// Read number retained species 
			std::cout << "Reading number of retained species [number-retained-species]" << std::endl;
			nretspecies = ptree.get<unsigned int>("opensmoke.number-retained-species");
			std::cout << " * Number of retained species: " << nretspecies << std::endl;

			// Read list retained species 
			if (nretspecies != thermodynamicsMap->NumberOfSpecies())
			{
				std::cout << "Reading list of retained species [list-retained-species]" << std::endl;

				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.list-retained-species"));
				listretspecies.resize(nretspecies);

				for (unsigned int j = 0; j < nretspecies; j++)
					stream >> listretspecies[j];	// (1-index based)

				std::cout << std::endl;
				std::cout << " * Retained species: " << nretspecies << std::endl;
				for (unsigned int j = 0; j < nretspecies; j++)
					std::cout << "  " << thermodynamicsMap->NamesOfSpecies()[listretspecies[j] - 1] << std::endl;
				std::cout << std::endl;
			}

			group.resize(ndata);
			csi.resize(ndata);
			T.resize(ndata);
			omega.resize(ndata, ns);
			belonging.resize(nclusters);

			if (read_pca_data == true)
			{
				mu.resize(nf + 1);
				sigma.resize(nf + 1);
				w.resize(nf + 1, npca);
				pca.resize(ndata, npca);
			}

			// Read original data
			{
				std::cout << "Reading original data [data-original]" << std::endl;
				std::stringstream stream;
				stream.str(ptree.get< std::string >("opensmoke.data-original"));

				for (unsigned int j = 0; j < ndata; j++)
				{
					stream >> group(j);
					stream >> csi(j);
					stream >> T(j);

					for (int i = 0; i < ns; i++)
						stream >> omega(j, i);

					belonging[group(j) - 1].push_back(j);
				}

				std::cout << " * Distribution of items among groups: " << nf << std::endl;
				for (unsigned int j = 0; j < nclusters; j++)
					std::cout << "    + "<<  j << " " << belonging[j].size() << " " << belonging[j].size() / double(ndata)*100. << "%" << std::endl;
			}

			if (read_pca_data == true)
			{
				// Read mean values of filtered variables
				{
					std::cout << "Reading mean values [mu]" << std::endl;
					std::stringstream stream;
					stream.str(ptree.get< std::string >("opensmoke.mu"));

					for (unsigned int j = 0; j < nf + 1; j++)
						stream >> mu(j);
				}
			
				// Read std deviations of filtered variables
				{
					std::cout << "Reading std deviations [sigma]" << std::endl;
					std::stringstream stream;
					stream.str(ptree.get< std::string >("opensmoke.sigma"));

					for (unsigned int j = 0; j < nf + 1; j++)
						stream >> sigma(j);
				}

				// Read PCA weights
				{
					std::cout << "Reading PCA weights [weights]" << std::endl;
					std::stringstream stream;
					stream.str(ptree.get< std::string >("opensmoke.weights"));

					for (unsigned int j = 0; j < nf + 1; j++)
						for (unsigned int i = 0; i < npca; i++)
							stream >> w(j, i);
				}

				// Read PCA
				{
					std::cout << "Reading PCA data [data-pca]" << std::endl;
					std::stringstream stream;
					stream.str(ptree.get< std::string >("opensmoke.data-pca"));

					for (unsigned int j = 0; j < ndata; j++)
					{
						unsigned int dummy;
						stream >> dummy;
						for (int i = 0; i < npca; i++)
							stream >> pca(j, i);
					}
				}
			}

			// Check PCA reconstruction (TOADJUST)
			/*
			{
				Eigen::VectorXd original(nf + 1);
				Eigen::VectorXd transformed(npca);

				for (unsigned int k = 0; k < 1000; k += 5)
				{
					std::cout << "Point " << k << std::endl;

					original(0) = (T(k)-mu(0))/sigma(0);
					for (unsigned int i = 0; i < ns; i++)
						original(i + 1) = (omega(k, i)-mu(i+1))/sigma(i+1);
					transformed = original.transpose()*w;

					for (unsigned int i = 0; i < npca; i++)
						std::cout << transformed(i) << " " << pca(k, i) << " " << transformed(i) - pca(k, i) << std::endl;

					getchar();
				}
			}
			*/
		}

		if (iDRG == true)
		{
			std::cout << "Applying DRG..." << std::endl;

			boost::filesystem::path path_drg_output_folder = path_output_folder / "DRG";
			if (!boost::filesystem::exists(path_drg_output_folder))
				OpenSMOKE::CreateDirectory(path_drg_output_folder);

			// Preparing Analysis
			Eigen::VectorXi number_important_species(ndata);
			Eigen::MatrixXi important_species(ndata, ns); important_species.setZero();
			
			// Analysis
			OpenSMOKE::OpenSMOKEVectorDouble y(ns);
			OpenSMOKE::OpenSMOKEVectorDouble x(ns);
			OpenSMOKE::OpenSMOKEVectorDouble c(ns);

			OpenSMOKE::DRG drg(thermodynamicsMap, kineticsMap);
			for (int j = 0; j < ndata; j++)
			{
				drg.SetKeySpecies(key_species);
				drg.SetEpsilon(epsilon);
				//drg.SetEpsilon(EpsilonDRG(T(j)));		// TODO
				drg.SetTemperatureThreshold(T_threshold);

				for (int k = 0; k < ns; k++)
					y[k + 1] = omega(j, k);

				double MW;
				thermodynamicsMap->MoleFractions_From_MassFractions(x.GetHandle(), MW, y.GetHandle());

				const double cTot = P / PhysicalConstants::R_J_kmol / T(j);
				OpenSMOKE::Product(cTot, x, &c);
				for (int k = 1; k <= ns; k++)
					if (c(k) < c_threshold)	c(k) = 0.;

				drg.Analysis(T(j), P, c);

				number_important_species(j) = drg.number_important_species();
				for (int i = 0; i < number_important_species(j); ++i)
				{
					const unsigned int k = drg.indices_important_species()[i];
					important_species(j, k) = 1;
				}

				if ((j + 1) % 1000 == 1)
					std::cout << j << "/" << ndata << " T: " << T(j) << " Species: " << number_important_species(j) << " Eps: " << epsilon << std::endl;
	//				std::cout << j << "/" << ndata << " T: " << T(j) << " Species: " << number_important_species(j) << " Eps: " << EpsilonDRG(T(j)) << std::endl;
			}

			// Post processing analysis

			// Matrix of presence of species
			PrintMatrixOfPresence(	path_drg_output_folder, ndata, 
									number_important_species, important_species);

			// Error analysis
			Eigen::MatrixXi retained_species;
			PrintErrorAnalysis(	path_drg_output_folder, ndata, nclusters,
								number_important_species,important_species,
								belonging, thermodynamicsMap->NamesOfSpecies(),
								retained_threshold,
								retained_species);

			// Select reactions
			Eigen::MatrixXi retained_reactions;
			SelectImportantReactions(*kineticsMap, retained_species, retained_reactions);

			// Write kinetic mechanisms
			WriteKineticMechanisms(*thermodynamicsMap, *kineticsMap,
									path_drg_output_folder, chemkin_thermodynamics_file, chemkin_kinetics_file,
									retained_species, retained_reactions);

			// Preprocess kinetic mechanisms
			if (iTransport == true)
			{
				PreprocessKineticMechanisms(nclusters, path_drg_output_folder,
											chemkin_thermodynamics_file, 
											chemkin_transport_file);
			}
		}

		if (iDRGEP == true)
		{
			std::cout << "Applying DRGEP..." << std::endl;

			boost::filesystem::path path_drgep_output_folder = path_output_folder / "DRGEP";
			if (!boost::filesystem::exists(path_drgep_output_folder))
				OpenSMOKE::CreateDirectory(path_drgep_output_folder);

			Eigen::VectorXi number_important_species(ndata);
			Eigen::MatrixXi important_species(ndata, ns); important_species.setZero();

			OpenSMOKE::OpenSMOKEVectorDouble y(ns);
			OpenSMOKE::OpenSMOKEVectorDouble x(ns);
			OpenSMOKE::OpenSMOKEVectorDouble c(ns);

			for (int jj = 1; jj <= nclusters; jj++)
			{
				unsigned int jLocal = 0;
				for (int j = 0; j < ndata; j++)
				{
					if (group(j) == jj)
					{
						OpenSMOKE::DRGEP drgep(thermodynamicsMap, kineticsMap);
						drgep.SetKeySpecies(key_species);
						drgep.SetEpsilon(epsilon);
						//drgep.SetEpsilon(EpsilonDRGEP(T(j)));	// TODO
						drgep.SetTemperatureThreshold(T_threshold);
						drgep.PrepareKineticGraph(key_species);

						for (int k = 0; k < ns; k++)
							y[k + 1] = omega(j, k);

						double MW;
						thermodynamicsMap->MoleFractions_From_MassFractions(x.GetHandle(), MW, y.GetHandle());

						const double cTot = P / PhysicalConstants::R_J_kmol / T(j);
						OpenSMOKE::Product(cTot, x, &c);
						for (int k = 1; k <= ns; k++)
							if (c(k) < c_threshold)	c(k) = 0;

						drgep.Analysis(T(j), P, c);

						number_important_species(j) = drgep.number_important_species();
						for (int i = 0; i < number_important_species(j); ++i)
						{
							const unsigned int k = drgep.indices_important_species()[i];
							important_species(j, k) = 1;
						}

						if (jLocal == 0)
							std::cout	<< "Group: " << jj << "/" << nclusters 
										<< " T: " << T(j) 
										<< " Species: " << number_important_species(j) 
										<< " Eps: " << epsilon << std::endl;
						
						if (jLocal % 1000 == 0)
							std::cout << " - processing " << jLocal + 1 << " out of " << belonging[jj - 1].size() << std::endl;

						jLocal++;
					}
				}
			}

			// Post processing analysis

			// Matrix of presence of species
			PrintMatrixOfPresence(	path_drgep_output_folder, ndata,
									number_important_species, important_species);

			// Error analysis
			Eigen::MatrixXi retained_species;
			PrintErrorAnalysis(	path_drgep_output_folder, ndata, nclusters,
								number_important_species, important_species,
								belonging, thermodynamicsMap->NamesOfSpecies(),
								retained_threshold,
								retained_species);

			// Select reactions
			Eigen::MatrixXi retained_reactions;
			SelectImportantReactions(*kineticsMap, retained_species, retained_reactions);

			// Write kinetic mechanisms
			WriteKineticMechanisms(*thermodynamicsMap, *kineticsMap,
									path_drgep_output_folder, chemkin_thermodynamics_file, chemkin_kinetics_file,
									retained_species, retained_reactions);

			// Preprocess kinetic mechanisms
			if (iTransport == true)
			{
				PreprocessKineticMechanisms(nclusters, path_drgep_output_folder,
											chemkin_thermodynamics_file,
											chemkin_transport_file);
			}
		}

		// Testing the network
		/*
		if (iTestingNeuralNetwork == true)
		{
			Eigen::VectorXd y(nclusters);

			int success = 0;
			for (unsigned int i = 0; i < ndata; i++)
			{
				Eigen::VectorXd row = pca.row(i);
				myNeuralNetworkBasedOnPCA(row.data(), y.data());
				
				int k;  
				y.maxCoeff(&k);
				if (k+1 == group(i))
					success++;

				// Example
				if (i % 1000 == 0)
				{
					std::cout << "Testing point: " << i << " - Group: " << group(i) << std::endl;
					std::cout << " * Returned group: " << k + 1 << std::endl;
					std::cout << " * Complete list " << std::endl;
					for (int j = 0; j < nclusters; j++)
						std::cout << j + 1 << " " << y(j) << std::endl;
				}
			}

			std::cout << "Summary" << std::endl;
			std::cout << success << " / " << ndata << " / " << double(success) / double(ndata)*100. << std::endl;
		}
		*/
	}

	bool iPostProcessFoam = false;
	if (iPostProcessFoam == true)
	{
		OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMap;

		// Reading thermodynamic and kinetic files
		{
			boost::property_tree::ptree ptree;
			boost::property_tree::read_xml((path_kinetics_folder / "kinetics.xml").string(), ptree);

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

			thermodynamicsMap = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);

			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			std::cout << " * Time to read XML file: " << tEnd - tStart << std::endl;
		}

		std::vector<std::string> list_labels;
		list_labels.push_back("zMix");
		list_labels.push_back("T");
		for (unsigned int i = 0; i < thermodynamicsMap->NumberOfSpecies(); i++)
			list_labels.push_back(thermodynamicsMap->NamesOfSpecies()[i]);

		std::vector<std::string> list_times;
		list_times.push_back("0");
		list_times.push_back("0.1");
		list_times.push_back("0.11");
		list_times.push_back("0.12");
		list_times.push_back("0.13");
		list_times.push_back("0.14");
		list_times.push_back("0.15");

		const unsigned int skip_lines = 22;
		const unsigned int npoints = 21334;
		const double T_threshold = 299.;

		boost::filesystem::path global_folder_name = "C:/Users/acuoci/Desktop/PCI 2018/Simulations/Unsteady/Detailed/";
		for (unsigned int k = 0; k < list_times.size(); k++)
		{
			std::cout << "Processing time: " << list_times[k] << std::endl;

			boost::filesystem::path folder_name = global_folder_name / list_times[k];
			
			Eigen::MatrixXd values(npoints, list_labels.size());
			for (unsigned int i = 0; i < list_labels.size(); i++)
			{
				boost::filesystem::path file_name = folder_name / list_labels[i];
				std::ifstream fInput(file_name.c_str(), std::ios::in);

				std::string dummy;
				for(unsigned j=0;j<skip_lines;j++)
					std::getline(fInput, dummy);

				for (unsigned j = 0; j < npoints; j++)
					fInput >> values(j, i);

				fInput.close();
			}

			boost::filesystem::path file_name = global_folder_name / ("output." + list_times[k]);
			std::ofstream fOutput(file_name.c_str(), std::ios::out);
			fOutput << "Data from OpenFOAM simulation" << std::endl;
			for (unsigned j = 0; j < npoints; j++)
			{
				if (values(j, 1) > T_threshold)
				{
					fOutput << std::setw(3) << std::fixed << 0;	// dummy col
					fOutput << std::setw(3) << std::fixed << 0;	// dummy col
					fOutput << std::setw(3) << std::fixed << 0;	// dummy col
					for (unsigned int i = 0; i < list_labels.size(); i++)
						fOutput << std::setw(15) << std::scientific << std::setprecision(6) << values(j, i);
					fOutput << std::endl;
				}
			}
			fOutput.close();
		}

	}

	std::cout << "Press enter to exit..." << std::endl;
	getchar();
	return 0;

}

