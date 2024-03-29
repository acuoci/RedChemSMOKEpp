/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_PolimiSoot_Analyzer_H
#define OpenSMOKE_PolimiSoot_Analyzer_H

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"
#include "preprocessing/PolimiSootClasses.h"

// Dictionary
#include "Grammar_PolimiSoot_Analyzer.h"

namespace OpenSMOKE
{
	//!  A class for analyzing soot formation and distribution based on the POLIMI kinetic mechanism
	/*!
	 A class for analyzing soot formation and distribution based on the POLIMI kinetic mechanism
	*/

	class PolimiSoot_Analyzer
	{

	public:

		enum class ThermophoreticEffectInEnthalpyFluxes { EXCLUDE_TOTAL_SOOT_CONTRIBUTION, EXCLUDE_THERMOPHORETIC_SOOT_CONTRIBUTION, DO_NOT_EXCLUDE_SOOT_CONTRIBUTION };

		enum class SootPlanckCoefficient { SOOT_PLANCK_COEFFICIENT_NONE, SOOT_PLANCK_COEFFICIENT_SMOOKE, SOOT_PLANCK_COEFFICIENT_KENT, SOOT_PLANCK_COEFFICIENT_SAZHIN, SOOT_PLANCK_COEFFICIENT_HUBBARD };

		/**
		*@brief Constructor based on the thermodynamic map
		*@param thermodynamicsMap map containing the thermodynamic data
		*/
		PolimiSoot_Analyzer(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML);

		/**
		*@brief Constructor based on the thermodynamic map and input from a user-defined dictionary
		*@param thermodynamicsMap map containing the thermodynamic data
		*@param dictionary the dictionary defined by the user
		*/
		PolimiSoot_Analyzer(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML, OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Prepares the object according to what defined in the external dictionary
		*@param dictionary the dictionary defined by the user
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Initializes the default values
		*/
		void Initialize();

		/**
		*@brief Reads the soot classes (if available) from the kinetics.xml file
		*@param XML file name
		*/
		void ClassesFromXMLFile(const boost::filesystem::path& file_name);

		/**
		*@brief Prepares the object with default options
		*/
		void Setup();

		/**
		*@brief Analysis of soot
		*@param T temperature [K]
		*@param P_Pa pressure [Pa]
		*@param rhoGas density of gaseous mixture [kg/m3]
		*@param omegaGas mass fractions of gaseous species
		*@param xGas mole fractions of gaseous species
		*/
		void Analysis(const double T, const double P_Pa, const double rhoGas, const Eigen::VectorXd &omegaGas, const Eigen::VectorXd &xGas);

		/**
		*@brief Calculations of the Particle Size Distribution Function (PSDF)
		*/
		void Distribution();
	

	public:		// query functions

		/**
		*@brief Returns the number of BINs
		*/
		unsigned int number_sections() const { return static_cast<unsigned int>(bin_indices_.size()); }

		/**
		*@brief Returns true if soot particles are really available in the kinetic mechanism
		*/
		bool is_active() const { return (bin_indices_.size() != 0); }

		/**
		*@brief Returns true if the PSDF (Particle Size Distribution Function) has to be written on the file
		*/
		bool write_psdf() const { return write_psdf_; }

		/**
		*@brief Returns the minimum soot volume fraction required for writing the PSDF on the file
		*/
		double threshold_for_psdf() const { return threshold_for_psdf_; }


	public:		// set options

		/**
		*@brief Sets the label to identify the soot particles
		*@param label the label identifying the soot particles
		*/
		void SetLabel(const std::string label);

		/**
		*@brief Sets the fractal diameter
		*@param Df the fractal diameter (default is 1.8)
		*/
		void SetFractalDiameter(const double Df);

		/**
		*@brief Sets the minimum section size for spherical soot particles
		*@param minimum_section minimum section size for spherical soot particles (default is 5)
		*/
		void SetMinimumSectionSphericalParticles(const int minimum_section_spheres);

		/**
		*@brief Sets the minimum section size for spherical soot particles
		*@param bin_minimumn minimum section size for spherical soot particles (default is BIN5)
		*/
		void SetMinimumSectionSphericalParticles(const std::string bin_minimum_spheres);

		/**
		*@brief Sets the minimum section size for soot aggregates
		*@param minimum_section minimum section size for soot aggregates (default is 13)
		*/
		void SetMinimumSectionAggregates(const int minimum_section_aggregates);

		/**
		*@brief Sets the minimum section size for soot aggregates
		*@param bin_minimumn minimum section size for soot aggregates (default is BIN13)
		*/
		void SetMinimumSectionAggregates(const std::string bin_minimum_aggrgates);

		/**
		*@brief Sets the linear law for calculating the density of soot particles
		*@param bin_index_zero initial bin
		*@param bin_index_final final bin
		*@param bin_density_zero density of initial bin [kg/m3]
		*@param bin_density_final density of final bin [kg/m3]
		*/
		void SetDensity(const int bin_index_zero, const int bin_index_final, const double bin_density_zero, const double bin_density_final);
		
		/**
		*@brief Sets the maximum BIN index at which the diffusivity reduction coefficient has to be applied
		*@param physical_diffusion_bin_to_cut the BIN index
		*/
		void SetPhysicalDiffusionBinToCut(const unsigned int physical_diffusion_bin_to_cut);

		/**
		*@brief Sets the minimum BIN index at which the diffusivity reduction coefficient has to be applied
		*@param physical_diffusion_bin_to_cut the BIN index
		*/
		void SetPhysicalDiffusionBinToStart(const unsigned int physical_diffusion_bin_to_start);

		/**
		*@brief Sets the minimum BIN index for which the thermophoretic effect has to be applied
		*@param thermophoretic_effect_minimum_bin the minimum BIN index
		*/
		void SetThermophoreticEffectMinimumBin(const unsigned int thermophoretic_effect_minimum_bin);

	public:		// soot and radiative heat transfer

		/**
		*@brief Sets the law for calculating the soot Planck mean absorption coefficient
		*@param soot_planck_coefficient the law for calculating the soot Planck mean absorption coefficient
		*/
		void SetPlanckAbsorptionCoefficient(const SootPlanckCoefficient soot_planck_coefficient);

		/**
		*@brief Sets the law for calculating the soot Planck mean absorption coefficient
		*@param soot_planck_coefficient the law for calculating the soot Planck 
		        mean absorption coefficient: Smooke (default) | Kent | Sazhin | Hubbard
		*/
		void SetPlanckAbsorptionCoefficient(const std::string soot_planck_coefficient);

		/**
		*@brief Returns the soot Planck mean absorption coefficient
		*@param rhoGas density of gaseous mixture [kg/m3]
		*@param T temperature of gaseous mixture [K]
		*@param omegaGas mass fractions of gaseous mixture
		*/
		double planck_coefficient(const double rhoGas, const double T, const Eigen::VectorXd &omegaGas) const;
	
		/**
		*@brief Returns true is soot has to be accounted for in radiative heat transfer
		*/
		bool radiative_heat_transfer() const { return radiative_heat_transfer_; }


	public:		// soot and thermophoretic effect 

		/**
		*@brief Returns true is the thermophoretic effect is turned on
		*/
		bool thermophoretic_effect() const { return thermophoretic_effect_; }

		/**
		*@brief Returns true if the thermophoretic velocity has to be included in the estimation of
		        correction velocity (diffusion fluxes)
		*/
		bool thermophoretic_effect_included_in_correction() const { return thermophoretic_effect_included_in_correction_;  }
		
		/**
		*@brief Returns the policy to adopt for managing the thermophoretic velocity in enthalpy fluxes
		*/
		ThermophoreticEffectInEnthalpyFluxes thermophoretic_effect_in_enthalpy_fluxes() const { return thermophoretic_effect_in_enthalpy_fluxes_; }
		
		/**
		*@brief Returns the smoothing time (in s) for thermophoretic velocity
		*/
		double thermophoretic_effect_smoothing_time() const { return thermophoretic_effect_smoothing_time_; }
		
		/**
		*@brief Returns the amplification factor for calculating the thermophoretic velocity
		*/
		double thermophoretic_effect_amplification_factor() const {return thermophoretic_effect_amplification_factor_; }
		
	public:		// soot and physical diffusion

		/**
		*@brief Returns true if soot particle diffusion coefficients are calculated using the 
		        kinetic theory gas gases with proper extrapolation
		*/
		bool physical_diffusion() const { return physical_diffusion_; };

		/**
		*@brief Returns the reduction coefficient for diffusion reduction
		*/
		double physical_diffusion_reduction_coefficient() const { return physical_diffusion_reduction_coefficient_; }

		/**
		*@brief Returns the BIN index where to cut for diffusion reduction
		*/
		int physical_diffusion_bin_to_cut() const { return physical_diffusion_bin_to_cut_; }

		/**
		*@brief Returns the BIN index where to start for diffusion reduction
		*/
		int physical_diffusion_bin_to_start() const { return physical_diffusion_bin_to_start_; }

		/**
		*@brief Returns the exponent in the reduction law
		*/
		double physical_diffusion_exp() const { return physical_diffusion_exp_; }

		/**
		*@brief Returns the correction factors for diffusivity of soot particles
		*/
		const std::vector<double>& bin_physical_diffusion_correction_factors() const { return bin_physical_diffusion_correction_factors_; }
		
		/**
		*@brief Returns the reference BIN adopted for correcting the diffusivities of BIN particles
		*/
		unsigned int bin_physical_diffusion_reference_species() const { return bin_physical_diffusion_reference_species_; }


	public:		// small particles

		/**
		*@brief Returns the volume fraction of small soot particles
		*/
		double fv_small() const { return fv_small_; }

		/**
		*@brief Returns the partial density (in kg/m3) of small soot particles
		*/
		double rho_small() const { return rho_small_; }

		/**
		*@brief Returns the number of particles (in #/m3) of small soot particles
		*/
		double N_small() const { return N_small_; }

		/**
		*@brief Returns the mass fraction of small soot particles
		*/
		double omega_small() const { return omega_small_; }

		/**
		*@brief Returns the mole fraction of small soot particles
		*/
		double x_small() const { return x_small_; }

		/**
		*@brief Returns the H/C ratio of small soot particles
		*/
		double h_over_c_small() const { return h_over_c_small_; }

		/**
		*@brief Returns the O/C ratio of small soot particles
		*/
		double o_over_c_small() const { return o_over_c_small_; }

		/**
		*@brief Returns the O/H ratio of small soot particles
		*/
		double o_over_h_small() const { return o_over_h_small_; }


	public:		// large particles

		/**
		*@brief Returns the volume fraction of large soot particles
		*/
		double fv_large() const { return fv_large_; }

		/**
		*@brief Returns the partial density (in kg/m3) of large soot particles
		*/
		double rho_large() const { return rho_large_; }

		/**
		*@brief Returns the number density (in #/m3) of large soot particles
		*/
		double N_large() const { return N_large_; }

		/**
		*@brief Returns the mass fraction of large soot particles
		*/
		double omega_large() const { return omega_large_; }

		/**
		*@brief Returns the mole fraction of large soot particles
		*/
		double x_large() const { return x_large_; }

		/**
		*@brief Returns the H/C ratio of large soot particles
		*/
		double h_over_c_large() const { return h_over_c_large_; }

		/**
		*@brief Returns the O/C ratio of large soot particles
		*/
		double o_over_c_large() const { return o_over_c_large_; }

		/**
		*@brief Returns the O/H ratio of large soot particles
		*/
		double o_over_h_large() const { return o_over_h_large_; }


	public:		// large spherical particles (i.e. from BIN5 to BIN11 by default)

		/**
		*@brief Returns the volume fraction for large particles with spherical shape (i.e. from BIN5 to BIN11 by default)
		*/
		double fv_large_spherical() const { return fv_large_spherical_; }

		/**
		*@brief Returns the partial density (in kg/m3) for large particles with spherical shape (i.e. from BIN5 to BIN11 by default)
		*/
		double rho_large_spherical() const { return rho_large_spherical_; }

		/**
		*@brief Returns the number density (in #/m3) for large particles with spherical shape (i.e. from BIN5 to BIN11 by default)
		*/
		double N_large_spherical() const { return N_large_spherical_; }

		/**
		*@brief Returns the mass fraction for large particles with spherical shape (i.e. from BIN5 to BIN11 by default)
		*/
		double omega_large_spherical() const { return omega_large_spherical_; }

		/**
		*@brief Returns the mole fraction for large particles with spherical shape (i.e. from BIN5 to BIN11 by default)
		*/
		double x_large_spherical() const { return x_large_spherical_; }

		/**
		*@brief Returns the H/C ratio of spherical soot particles
		*/
		double h_over_c_large_spherical() const { return h_over_c_large_spherical_; }

		/**
		*@brief Returns the O/C ratio of large spherical soot particles
		*/
		double o_over_c_large_spherical() const { return o_over_c_large_spherical_; }

		/**
		*@brief Returns the O/H ratio of large spherical soot particles
		*/
		double o_over_h_large_spherical() const { return o_over_h_large_spherical_; }


	public:		// large aggregates (i.e. from BIN12 by default)

		/**
		*@brief Returns the volume fraction for large aggregates (i.e. from BIN12 by default)
		*/
		double fv_large_aggregates() const { return fv_large_aggregates_; }

		/**
		*@brief Returns the partial density (in kg/m3) for large aggregates (i.e. from BIN12 by default)
		*/
		double rho_large_aggregates() const { return rho_large_aggregates_; }

		/**
		*@brief Returns the number density (#/m3) for large aggregates (i.e. from BIN12 by default)
		*/
		double N_large_aggregates() const { return N_large_aggregates_; }

		/**
		*@brief Returns the mass fraction for large aggregates (i.e. from BIN12 by default)
		*/
		double omega_large_aggregates() const { return omega_large_aggregates_; }

		/**
		*@brief Returns the mole fraction for large aggregates (i.e. from BIN12 by default)
		*/
		double x_large_aggregates() const { return x_large_aggregates_; }

		/**
		*@brief Returns the H/C ratio of aggregates
		*/
		double h_over_c_large_aggregates() const { return h_over_c_large_aggregates_; }

		/**
		*@brief Returns the O/C ratio of aggregates
		*/
		double o_over_c_large_aggregates() const { return o_over_c_large_aggregates_; }

		/**
		*@brief Returns the O/H ratio of aggregates
		*/
		double o_over_h_large_aggregates() const { return o_over_h_large_aggregates_; }


	public:		// PAHs

		/**
		*@brief Returns the mass fraction of PAHs with 1 or 2 aromatic rings
		*/
		double omega_pah_1_2_rings() const { return omega_pah_1_2_rings_; }

		/**
		*@brief Returns the mass fraction of PAHs with 3 or 4 aromatic rings
		*/
		double omega_pah_3_4_rings() const { return omega_pah_3_4_rings_; }
		
		/**
		*@brief Returns the mass fraction of PAHs with more than 4 aromatic rings
		*/
		double omega_pah_large() const { return omega_pah_large_; }


		/**
		*@brief Returns the mole fraction of PAHs with 1 or 2 aromatic rings
		*/
		double x_pah_1_2_rings() const { return x_pah_1_2_rings_; }

		/**
		*@brief Returns the mole fraction of PAHs with 3 or 4 aromatic rings
		*/
		double x_pah_3_4_rings() const { return x_pah_3_4_rings_; }

		/**
		*@brief Returns the mole fraction of PAHs with more than 4 aromatic rings
		*/
		double x_pah_large() const { return x_pah_large_; }

		/**
		*@brief Returns the volume fraction of PAHs with 1 or 2 aromatic rings
		*/
		double fv_pah_1_2_rings() const { return fv_pah_1_2_rings_; }

		/**
		*@brief Returns the volume fraction of PAHs with 3 or 4 aromatic rings
		*/
		double fv_pah_3_4_rings() const { return fv_pah_3_4_rings_; }

		/**
		*@brief Returns the volume fraction of PAHs with more than 4 aromatic rings
		*/
		double fv_pah_large() const { return fv_pah_large_; }

		/**
		*@brief Returns the partial density (in kg/m3) of PAHs with 1-2 rings
		*/
		double rho_pah_1_2_rings() const { return rho_pah_1_2_rings_; }

		/**
		*@brief Returns the number density (in #/m3) of PAHs with 1-2 rings
		*/
		double N_pah_1_2_rings() const { return N_pah_1_2_rings_; }

		/**
		*@brief Returns the partial density (in kg/m3) of PAHs with 3-4 rings
		*/
		double rho_pah_3_4_rings() const { return rho_pah_3_4_rings_; }

		/**
		*@brief Returns the number density (in #/m3) of PAHs with 3-4 rings
		*/
		double N_pah_3_4_rings() const { return N_pah_3_4_rings_; }

	public:		// on-the-fly evaluation

		/**
		*@brief Returns the soot volume fraction of large soot particles
		*@param rhoGas gas density [kg/m3]
		*@param omegaGas mass fractions of gaseous species
		*/
		double fv_large(const double rhoGas, const Eigen::VectorXd &omegaGas) const;

		/**
		*@brief Returns the soot volume fraction of small soot particles
		*@param rhoGas gas density [kg/m3]
		*@param omegaGas mass fractions of gaseous species
		*/
		double fv_small(const double rhoGas, const Eigen::VectorXd &omegaGas) const;

		/**
		*@brief Returns the total (small+large) soot volume fraction
		*@param rhoGas gas density [kg/m3]
		*@param omegaGas mass fractions of gaseous species
		*/
		double fv(const double rhoGas, const Eigen::VectorXd &omegaGas) const;

		/**
		*@brief Returns the partial density (in kg/m3) of large soot particles
		*@param rhoGas gas density [kg/m3]
		*@param omegaGas mass fractions of gaseous species
		*/
		double rho_large(const double rhoGas, const Eigen::VectorXd &omegaGas) const;

		/**
		*@brief Returns the mass fraction of PAHs with 1 or 2 aromatic rings
		*/
		double omega_pah_1_2_rings(const Eigen::VectorXd &omegaGas) const;

		/**
		*@brief Returns the mass fraction of PAHs with 3 or 4 aromatic rings
		*/
		double omega_pah_3_4_rings(const Eigen::VectorXd &omegaGas) const;

		/**
		*@brief Returns the mass fraction of PAHs with more than 4 aromatic rings
		*/
		double omega_pah_large(const Eigen::VectorXd &omegaGas) const;

		/**
		*@brief Returns the mole fraction of PAHs with 1 or 2 aromatic rings
		*/
		double x_pah_1_2_rings(const Eigen::VectorXd& xGas) const;

		/**
		*@brief Returns the mole fraction of PAHs with 3 or 4 aromatic rings
		*/
		double x_pah_3_4_rings(const Eigen::VectorXd& xGas) const;

		/**
		*@brief Returns the mole fraction of PAHs with more than 4 aromatic rings
		*/
		double x_pah_large(const Eigen::VectorXd& xGas) const;

	public:		// writing on files

		/**
		*@brief Writes the summary of the soot model on a file
		*/
		void WriteBinData();

		/**
		*@brief Writes the Particle Size Distribution Function (PSDF) on file
		*@param t time [s]
		*@param x x spatial coordinate [m]
		*@param y y spatial coordinate [m]
		*@param z z spatial coordinate [m]
		*@param temperature temperature [K]
		*/
		void WriteDistribution(std::ofstream& fSootDistribution, const double t, const double x, const double y, const double z, const double temperature);

		/**
		*@brief Writes head line for file on which the soot particle size distribution function will be written
		*@param fSoot reference to the file
		*/
		void WriteDistributionLabel(std::ofstream &fSoot);

		/**
		*@brief Writes the integral quantities about soot on a file
		*@param fOutput reference to the file
		*@param width width of colums
		*/
		void WriteIntegralDataFile(std::ofstream& fOutput, const unsigned int width);

		/**
		*@brief Writes the head line for file on which integral quantities will be written
		*@param fOutput reference to the file
		*@param counter current counter
		*@param width width of columns
		*/
		void WriteLabelIntegralDataFile(std::ofstream& fOutput, unsigned int& counter, const unsigned int width);

		/**
		*@brief Writes the formation rates of soot classes
		*@param fOutput reference to the file
		*@param width width of columns
		*/
		void WriteFormationSootClassesDataFile(std::ofstream& fOutput, const unsigned int width);

		/**
		*@brief Writes the head line for file on which formation rates of soot classes will be written
		*@param fOutput reference to the file
		*@param counter current counter
		*@param width width of columns
		*/
		void WriteLabelFormationSootClassesDataFile(std::ofstream& fOutput, unsigned int& counter, const unsigned int width);
		

		/**
		*@brief Return the correction coefficients for the soot reactions to be modified
		*@param indices reaction indices (0-based)
		*@param coefficients correction coefficients
		*@return true if correction coefficients are enabled/available
		*/
		bool ClassesCorrectionCoefficients(std::vector<unsigned int>& indices, std::vector<double>& coefficients);

		void Analysis(	OpenSMOKE::KineticsMap_CHEMKIN& kinetics,
						const double T, const double P_Pa, const double rhoGas, const Eigen::VectorXd& omegaGas, const Eigen::VectorXd& xGas);

		PolimiSootClasses* soot_classes() { return soot_classes_; }

		bool formation_rates_by_classes() const { return formation_rates_by_classes_;  }

	private:

		/**
		*@brief Calculates soot diameters (dmean = d^a/d^b)
		*@param a parameter in the mean diameter definition
		*@param b parameter in the mean diameter definition
		*@param bin_indices indices of BINs to be accounted for in the diameter evaluation
		*@param threshold minimum mole fraction to activate calculation
		*@return the mean diameter (in m)
		*/
		double SootDiameters(const double a, const double b, const std::vector<unsigned int>& bin_indices, const double threshold);

	public:		// indices

		const std::vector<unsigned int>& bin_indices() const { return bin_indices_; }
		const std::vector<std::string>& bin_names() const { return bin_names_; }
		const std::vector<unsigned int>& bin_indices_large() const { return bin_indices_large_; }
		const std::vector<unsigned int>& bin_indices_small() const { return bin_indices_small_; }
		const std::vector<unsigned int>& bin_indices_large_spherical() const { return bin_indices_large_spherical_; }
		const std::vector<unsigned int>& bin_indices_large_aggregates() const { return bin_indices_large_aggregates_; }

		const std::vector<unsigned int>& bin_indices_large_global() const { return bin_indices_large_global_; }
		const std::vector<unsigned int>& bin_indices_small_global() const { return bin_indices_small_global_; }
		const std::vector<unsigned int>& bin_indices_large_spherical_global() const { return bin_indices_large_spherical_global_; }
		const std::vector<unsigned int>& bin_indices_large_aggregates_global() const { return bin_indices_large_aggregates_global_; }

		const std::vector<unsigned int>& soot_dimer_indices_global() const { return soot_dimer_indices_global_;   }
		const std::vector<unsigned int>& pah_1_2_rings_indices_global() const { return pah_1_2_rings_indices_global_;   }
		const std::vector<unsigned int>& pah_3_4_rings_indices_global() const { return pah_3_4_rings_indices_global_;   }
		const std::vector<unsigned int>& pah_large_indices_global() const { return pah_1_2_rings_indices_global_; }

		const std::vector<double>& bin_density_large() const { return bin_density_large_; }
		const std::vector<double>& bin_V_large() const { return bin_V_large_; }
		const std::vector<double>& bin_density_small() const { return bin_density_small_; }
		const std::vector<double>& bin_V_small() const { return bin_V_small_; }

		const std::vector<double>& bin_baskets_d() const {return bin_baskets_d_;}
		const std::vector<double>& dN_over_dlog10d() const {return dN_over_dlog10d_;}

		const std::vector<unsigned int>& bin_indices_thermophoresis() const { return bin_indices_thermophoresis_; }

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermo_;

		std::string bin_label_;
		std::string bin_minimum_spheres_;
		std::string bin_minimum_aggregates_;
		int bin_minimum_spheres_section_;
		int bin_minimum_aggregates_section_;

		unsigned int bin_density_index_zero_;
		unsigned int bin_density_index_final_;
		double bin_density_value_zero_;
		double bin_density_value_final_;

		double Df_;

		bool thermophoretic_effect_;
		bool thermophoretic_effect_included_in_correction_;
		ThermophoreticEffectInEnthalpyFluxes thermophoretic_effect_in_enthalpy_fluxes_;
		double thermophoretic_effect_smoothing_time_;
		double thermophoretic_effect_amplification_factor_;
		int thermophoretic_effect_minimum_bin_;
		
		bool radiative_heat_transfer_;
		bool physical_diffusion_;

		unsigned int iC_;
		unsigned int iO_;
		unsigned int iH_;
		unsigned nspecies_;

		std::vector<double> bin_physical_diffusion_correction_factors_;
		int bin_physical_diffusion_reference_species_;

		std::vector<double> bin_density_;
		std::vector<unsigned int> bin_indices_;
		std::vector<std::string> bin_names_;
		std::vector<double> bin_mw_;
		std::vector<double> bin_m_;
		std::vector<double> bin_ds_;
		std::vector<unsigned int> bin_section_;				// BIN section [-]
		std::vector<bool> bin_is_radical_;					// BIN is radical
		std::vector<std::string> bin_hydrogenation_level_;	// BIN hydrogenation level
		
		std::vector<std::string> bin_hydrogenation_level_unique_values_;	// BIN hydrogenation level (unique values)

	
		std::vector<std::vector<unsigned int>> bin_indices_sections_;
		std::vector<std::vector<unsigned int>> bin_indices_sections_global_;

		std::vector<double> bin_V_;
		std::vector<double> bin_V_large_;
		std::vector<double> bin_V_small_;

		std::vector<double> bin_c_;
		std::vector<double> bin_h_;
		std::vector<double> bin_o_;
		std::vector<double> bin_h_over_c_;
		std::vector<double> bin_o_over_c_;
		std::vector<double> bin_o_over_h_;

		std::vector<unsigned int> bin_indices_large_;
		std::vector<unsigned int> bin_indices_small_;
		std::vector<unsigned int> bin_indices_large_global_;
		std::vector<unsigned int> bin_indices_small_global_;
		std::vector<unsigned int> bin_indices_large_spherical_;
		std::vector<unsigned int> bin_indices_large_aggregates_;
		std::vector<unsigned int> bin_indices_large_spherical_global_;
		std::vector<unsigned int> bin_indices_large_aggregates_global_;
		std::vector<unsigned int> bin_indices_thermophoresis_;

		std::vector<double> bin_np_;
		std::vector<double> bin_dc_;
		std::vector<double> bin_d_;

		std::vector<double> bin_omega_;
		std::vector<double> bin_x_;
		std::vector<double> bin_fv_;
		std::vector<double> bin_rho_;
		std::vector<double> bin_N_;

		std::vector<double> bin_density_small_;
		std::vector<double> bin_density_large_;

		std::vector< std::vector<unsigned int> >bin_baskets_indices_;

		std::vector<double> bin_baskets_c_;
		std::vector<double> bin_baskets_d_;
		std::vector<double> bin_baskets_mw_;
		std::vector<double> bin_baskets_log10d_;
		std::vector<double> bin_baskets_dlog10d_;
		std::vector<double> dN_over_dlog10d_;

		std::vector<double> bin_baskets_V_;
		std::vector<double> bin_baskets_log10V_;
		std::vector<double> bin_baskets_dlog10V_;
		std::vector<double> dN_over_dlog10V_;

		std::vector<double> bin_baskets_m_;
		std::vector<double> bin_baskets_log10m_;
		std::vector<double> bin_baskets_dlog10m_;
		std::vector<double> dN_over_dlog10m_;

		std::vector<double> bin_baskets_N_;
		std::vector<double> bin_baskets_fv_;
		std::vector<double> bin_baskets_rho_;
		std::vector<double> bin_baskets_x_;
		std::vector<double> bin_baskets_omega_;

		std::vector<unsigned int> soot_dimer_indices_global_;
		std::vector<unsigned int> pah_1_2_rings_indices_global_;
		std::vector<unsigned int> pah_3_4_rings_indices_global_;
		std::vector<unsigned int> pah_large_indices_global_;

		double fv_small_;
		double rho_small_;
		double N_small_;
		double omega_small_;
		double x_small_;
		double h_over_c_small_;
		double o_over_c_small_;
		double o_over_h_small_;

		double fv_large_;
		double rho_large_;
		double N_large_;
		double omega_large_;
		double x_large_;
		double h_over_c_large_;
		double o_over_c_large_;
		double o_over_h_large_;

		double fv_large_spherical_;
		double rho_large_spherical_;
		double N_large_spherical_;
		double omega_large_spherical_;
		double x_large_spherical_;
		double h_over_c_large_spherical_;
		double o_over_c_large_spherical_;
		double o_over_h_large_spherical_;

		double fv_large_aggregates_;
		double rho_large_aggregates_;
		double N_large_aggregates_;
		double omega_large_aggregates_;
		double x_large_aggregates_;
		double h_over_c_large_aggregates_;
		double o_over_c_large_aggregates_;
		double o_over_h_large_aggregates_;

		double fv_pah_1_2_rings_;
		double fv_pah_3_4_rings_;
		double fv_pah_large_;

		double omega_pah_1_2_rings_;
		double omega_pah_3_4_rings_;
		double omega_pah_large_;

		double x_pah_1_2_rings_;
		double x_pah_3_4_rings_;
		double x_pah_large_;

		double N_pah_1_2_rings_;
		double N_pah_3_4_rings_;

		double rho_pah_1_2_rings_;
		double rho_pah_3_4_rings_;

		double d10_N_large_;
		double d32_N_large_;
		double dvariance_N_;
		double dstd_N_;

		double d10_small_;
		double d32_small_;
		double d43_small_;

		double d10_large_;
		double d32_large_;
		double d43_large_;

		double d10_large_spherical_;
		double d32_large_spherical_;
		double d43_large_spherical_;

		double d10_large_aggregates_;
		double d32_large_aggregates_;
		double d43_large_aggregates_;

		SootPlanckCoefficient soot_planck_coefficient_;
		double physical_diffusion_reduction_coefficient_;
		int physical_diffusion_bin_to_start_;
		int physical_diffusion_bin_to_cut_;
		double physical_diffusion_exp_;

		bool write_psdf_;
		double threshold_for_psdf_;

		// Reaction classes
		PolimiSootClasses* soot_classes_;
		bool formation_rates_by_classes_;
		std::vector<std::string> list_class_corrections_names_;
		std::vector<double> list_class_corrections_coefficients_;
		Eigen::VectorXd R_sum_small_;
		Eigen::VectorXd R_sum_large_;
		Eigen::VectorXd R_sum_large_spherical_;
		Eigen::VectorXd R_sum_large_aggregates_;
		Eigen::VectorXd Omega_sum_small_;
		Eigen::VectorXd Omega_sum_large_;
		Eigen::VectorXd Omega_sum_large_spherical_;
		Eigen::VectorXd Omega_sum_large_aggregates_;
	};
}

#include "OpenSMOKE_PolimiSoot_Analyzer.hpp"

#endif /* OpenSMOKE_PolimiSoot_Analyzer_H */
