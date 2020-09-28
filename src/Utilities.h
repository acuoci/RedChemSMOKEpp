#ifndef RedChemSMOKEpp_Utilities
#define RedChemSMOKEpp_Utilities

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

double EpsilonDRG(const double T);

double EpsilonDRGEP(const double T);

void PrintMatrixOfPresence(const boost::filesystem::path& path_output_folder,
	const unsigned int ndata,
	const Eigen::VectorXi& number_important_species,
	const Eigen::MatrixXi& important_species);

void PrintErrorAnalysis(const boost::filesystem::path& path_output_folder,
	const unsigned int ndata, const unsigned int nclusters,
	const Eigen::VectorXi& number_important_species,
	const Eigen::MatrixXi& important_species,
	const std::vector< std::vector<int> > belonging,
	const std::vector<std::string> species_names,
	const double retained_threshold,
	Eigen::MatrixXi& retained_species);

void WriteKineticMechanisms(const OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
	const OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
	const boost::filesystem::path& path_output_folder,
	const boost::filesystem::path& chemkin_thermodynamic_file,
	const boost::filesystem::path& chemkin_kinetics_file,
	const Eigen::MatrixXi& retained_species,
	const Eigen::MatrixXi& retained_reactions);

void SelectImportantReactions(OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
	Eigen::MatrixXi& retained_species,
	Eigen::MatrixXi& retained_reactions);

void PreprocessKineticMechanisms(const unsigned int nclusters,
	const boost::filesystem::path& path_output_folder,
	const boost::filesystem::path& chemkin_thermodynamic_file,
	const boost::filesystem::path& chemkin_transport_file);

#include "Utilities.hpp"

#endif /* RedChemSMOKEpp_Utilities */