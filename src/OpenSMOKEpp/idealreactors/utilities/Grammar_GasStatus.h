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
|   License                                                               |
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

#ifndef OpenSMOKE_Grammar_GasStatus_H
#define	OpenSMOKE_Grammar_GasStatus_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{

	class Grammar_GasStatus : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Temperature", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Temperature of the mixture (i.e. 500 K)", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Pressure", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Pressure of the mixture (i.e. 1 atm)", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Density", 
																OpenSMOKE::SINGLE_MEASURE, 
																"Density of the mixture (i.e. 1 g/cm3)", 
																false) );	

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MoleFractions", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Mole fractions of the mixture (i.e. CH4 0.60 H2 0.40)", 
																true,
																"@MassFractions @Moles @Masses @EquivalenceRatio",
																"none",
																"none") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MassFractions", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Mass fractions of the mixture (i.e. CH4 0.60 H2 0.40)", 
																true,
																"@MoleFractions @Moles @Masses @EquivalenceRatio",
																"none",
																"none") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Moles", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Moles (relative) of the mixture (i.e. CH4 2 H2 1)", 
																true,
																"@MoleFractions @MassFractions @Masses @EquivalenceRatio",
																"none",
																"none") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Masses", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Masses (relative) of the mixture (i.e. CH4 2 H2 1)", 
																true,
																"@MoleFractions @MassFractions @Moles @EquivalenceRatio",
																"none",
																"none") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@EquivalenceRatio", 
																OpenSMOKE::VECTOR_DOUBLE, 
																"EquivalenceRatio", 
																true,
																"@MoleFractions @MassFractions @Moles @Masses",
																"none",
																"none") );
		
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FuelMoleFractions", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Fuel mole fractions", 
																false,
																"none",
																"@EquivalenceRatio",
																"@FuelMassFractions @FuelMoles @FuelMasses") );
		
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FuelMassFractions", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Fuel mass fractions", 
																false,
																"none",
																"@EquivalenceRatio",
																"@FuelMoleFractions @FuelMoles @FuelMasses") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FuelMoles", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Fuel moles", 
																false,
																"none",
																"@EquivalenceRatio",
																"@FuelMoleFractions @FuelMassFractions @FuelMasses") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FuelMasses", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Fuel masses", 
																false,
																"none",
																"@EquivalenceRatio",
																"@FuelMoleFractions @FuelMassFractions @FuelMoles") );

		
			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OxidizerMoleFractions", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Oxidizer mole fractions", 
																false,
																"none",
																"@EquivalenceRatio",
																"@OxidizerMassFractions @OxidizerMoles @OxidizerMasses") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OxidizerMassFractions", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Oxidizer mass fractions", 
																false,
																"none",
																"@EquivalenceRatio",
																"@OxidizerMoleFractions @OxidizerMoles @OxidizerMasses") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OxidizerMoles", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Oxidizer moles", 
																false,
																"none",
																"@EquivalenceRatio",
																"@OxidizerMoleFractions @OxidizerMassFractions @OxidizerMasses") );

			AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OxidizerMasses", 
																OpenSMOKE::VECTOR_STRING_DOUBLE, 
																"Oxidizer masses", 
																false,
																"none",
																"@EquivalenceRatio",
																"@OxidizerMoleFractions @OxidizerMassFractions @OxidizerMoles") );
															
		}
	};


	template<typename Thermodynamics>
	void GetGasStatusFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, Thermodynamics& thermodynamicsMapXML,
									double& T, double& P_Pa, OpenSMOKE::OpenSMOKEVectorDouble& omega)
	{
		Grammar_GasStatus grammar_gas_status;
		dictionary.SetGrammar(grammar_gas_status);

		unsigned int state_variables = 0;
		bool temperature_assigned = false;
		bool pressure_assigned = false;
		bool density_assigned = false;
		
		// Temperature
		{
			if (dictionary.CheckOption("@Temperature") == true)
			{
				double value;
				std::string units;
				dictionary.ReadMeasure("@Temperature", value, units);
		
				if (units == "K")			T = value;
				else if (units == "C")		T = value + 273.15;
				else OpenSMOKE::FatalErrorMessage("Unknown temperature units");

				state_variables++;
				temperature_assigned = true;
			}
		}

		// Pressure
		{
			if (dictionary.CheckOption("@Pressure") == true)
			{
				double value;
				std::string units;
				dictionary.ReadMeasure("@Pressure", value, units);

				if (units == "Pa")			P_Pa = value;
				else if (units == "bar")	P_Pa = value*1.e5;
				else if (units == "atm")	P_Pa = value*101325.;
				else OpenSMOKE::FatalErrorMessage("Unknown pressure units");

				state_variables++;
				pressure_assigned = true;
			}
		}

		// Density
		double rho;
		{
			if (dictionary.CheckOption("@Density") == true)
			{
				double value;
				std::string units;
				dictionary.ReadMeasure("@Density", value, units);

				if (units == "kg/m3")		rho = value;
				else if (units == "g/cm3")	rho = value*1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown density units");

				state_variables++;
				density_assigned = true;
			}
		}

		if (state_variables != 2)
			OpenSMOKE::FatalErrorMessage("The status of a gas mixture requires any 2 (and only 2) among: @Temperature, @Pressure and @Density");

		// Composition
		{
			{
				std::vector<std::string> names;
				std::vector<double> values;
		
				if (dictionary.CheckOption("@MoleFractions") == true)
				{
					dictionary.ReadOption("@MoleFractions", names, values);
				
					const double sum =std::accumulate(values.begin(),values.end(),0.);
					if (sum<(1.-1e-6) || sum>(1.+1e-6))
						OpenSMOKE::FatalErrorMessage("The mole fractions must sum to 1.");

					OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML.NumberOfSpecies());
					for(unsigned int i=0;i<names.size();i++)
						x[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i]/sum;
					ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega, true);
					double MW;
					thermodynamicsMapXML.MassFractions_From_MoleFractions(omega.GetHandle(),MW,x.GetHandle());
				}
				else if (dictionary.CheckOption("@MassFractions") == true)
				{
					dictionary.ReadOption("@MassFractions", names, values);

					const double sum =std::accumulate(values.begin(),values.end(),0.);
					if (sum<(1.-1e-6) || sum>(1.+1e-6))
						OpenSMOKE::FatalErrorMessage("The mass fractions must sum to 1.");

					ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega, true);
					for(unsigned int i=0;i<names.size();i++)
						omega[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i]/sum;
				}
				else if (dictionary.CheckOption("@Moles") == true)
				{
					dictionary.ReadOption("@Moles", names, values);
					const double sum =std::accumulate(values.begin(),values.end(),0.);
					OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML.NumberOfSpecies());
					for(unsigned int i=0;i<names.size();i++)
						x[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i]/sum;
					ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega, true);
					double MW;
					thermodynamicsMapXML.MassFractions_From_MoleFractions(omega.GetHandle(),MW,x.GetHandle());
				}
				else if (dictionary.CheckOption("@Masses") == true)
				{
					dictionary.ReadOption("@Masses", names, values);
					const double sum =std::accumulate(values.begin(),values.end(),0.);
					ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega, true);
					for(unsigned int i=0;i<names.size();i++)
						omega[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i]/sum;
				}
			}
		
			if (dictionary.CheckOption("@EquivalenceRatio") == true)
			{
				std::vector<double> equivalence_ratios;
				dictionary.ReadOption("@EquivalenceRatio", equivalence_ratios);

				if (equivalence_ratios.size() != 1)
					OpenSMOKE::FatalErrorMessage("Only a single @EquivalenceRatio can be specified");

				double equivalence_ratio = equivalence_ratios[0];
			
				std::vector<std::string> names_fuel;
				std::vector<double> values_fuel;

				if (dictionary.CheckOption("@FuelMoleFractions") == true)
				{
					dictionary.ReadOption("@FuelMoleFractions", names_fuel, values_fuel);

					const double sum =std::accumulate(values_fuel.begin(),values_fuel.end(),0.);
					if (sum<(1.-1e-6) || sum>(1.+1e-6))
						OpenSMOKE::FatalErrorMessage("The fuel mass fractions must sum to 1.");

					for(unsigned int i=0;i<names_fuel.size();i++)
						values_fuel[i] /= sum;
				}
				else if (dictionary.CheckOption("@FuelMassFractions") == true)
				{
					dictionary.ReadOption("@FuelMassFractions", names_fuel, values_fuel);
				
					const double sum =std::accumulate(values_fuel.begin(),values_fuel.end(),0.);
					if (sum<(1.-1e-6) || sum>(1.+1e-6))
						OpenSMOKE::FatalErrorMessage("The fuel mass fractions must sum to 1.");

					OpenSMOKE::OpenSMOKEVectorDouble omega_fuel(thermodynamicsMapXML.NumberOfSpecies());
					for(unsigned int i=0;i<names_fuel.size();i++)
						omega_fuel[thermodynamicsMapXML.IndexOfSpecies(names_fuel[i])] = values_fuel[i]/sum;

					double MW_fuel;
					OpenSMOKE::OpenSMOKEVectorDouble x_fuel(thermodynamicsMapXML.NumberOfSpecies());
					thermodynamicsMapXML.MoleFractions_From_MassFractions(x_fuel.GetHandle(), MW_fuel, omega_fuel.GetHandle());
					for(unsigned int i=0;i<names_fuel.size();i++)
						values_fuel[i] = x_fuel[thermodynamicsMapXML.IndexOfSpecies(names_fuel[i])];
				}
				else if (dictionary.CheckOption("@FuelMoles") == true)
				{
					dictionary.ReadOption("@FuelMoles", names_fuel, values_fuel);
				}
				else if (dictionary.CheckOption("@FuelMasses") == true)
				{
					dictionary.ReadOption("@FuelMasses", names_fuel, values_fuel);

					OpenSMOKE::OpenSMOKEVectorDouble omega_fuel(thermodynamicsMapXML.NumberOfSpecies());
					for(unsigned int i=0;i<names_fuel.size();i++)
						omega_fuel[thermodynamicsMapXML.IndexOfSpecies(names_fuel[i])] = values_fuel[i];

					double MW_fuel;
					OpenSMOKE::OpenSMOKEVectorDouble x_fuel(thermodynamicsMapXML.NumberOfSpecies());
					thermodynamicsMapXML.MoleFractions_From_MassFractions(x_fuel.GetHandle(), MW_fuel, omega_fuel.GetHandle());
					for(unsigned int i=0;i<names_fuel.size();i++)
						values_fuel[i] = x_fuel[thermodynamicsMapXML.IndexOfSpecies(names_fuel[i])];
				}
				else
				{
					OpenSMOKE::FatalErrorMessage("The @EquivalenceRatio option requires the user specifies the fuel composition");
				}

				std::vector<std::string> names_oxidizer;
				std::vector<double> values_oxidizer;

				if (dictionary.CheckOption("@OxidizerMoleFractions") == true)
				{
					dictionary.ReadOption("@OxidizerMoleFractions", names_oxidizer, values_oxidizer);

					const double sum =std::accumulate(values_oxidizer.begin(),values_oxidizer.end(),0.);
					if (sum<(1.-1e-6) || sum>(1.+1e-6))
						OpenSMOKE::FatalErrorMessage("The oxidizer mass fractions must sum to 1.");

					for(unsigned int i=0;i<names_oxidizer.size();i++)
						values_oxidizer[i] /= sum;
				}
				else if (dictionary.CheckOption("@OxidizerMassFractions") == true)
				{
					dictionary.ReadOption("@OxidizerMassFractions", names_oxidizer, values_oxidizer);
				
					const double sum =std::accumulate(values_oxidizer.begin(),values_oxidizer.end(),0.);
					if (sum<(1.-1e-6) || sum>(1.+1e-6))
						OpenSMOKE::FatalErrorMessage("The oxidizer mass fractions must sum to 1.");

					OpenSMOKE::OpenSMOKEVectorDouble omega_oxidizer(thermodynamicsMapXML.NumberOfSpecies());
					for(unsigned int i=0;i<names_oxidizer.size();i++)
						omega_oxidizer[thermodynamicsMapXML.IndexOfSpecies(names_oxidizer[i])] = values_oxidizer[i]/sum;

					double MW_oxidizer;
					OpenSMOKE::OpenSMOKEVectorDouble x_oxidizer(thermodynamicsMapXML.NumberOfSpecies());
					thermodynamicsMapXML.MoleFractions_From_MassFractions(x_oxidizer.GetHandle(), MW_oxidizer, omega_oxidizer.GetHandle());
					for(unsigned int i=0;i<names_oxidizer.size();i++)
						values_oxidizer[i] = x_oxidizer[thermodynamicsMapXML.IndexOfSpecies(names_oxidizer[i])];
				}
				else if (dictionary.CheckOption("@OxidizerMoles") == true)
				{
					dictionary.ReadOption("@OxidizerMoles", names_oxidizer, values_oxidizer);
				}
				else if (dictionary.CheckOption("@OxidizerMasses") == true)
				{
					dictionary.ReadOption("@OxidizerMasses", names_oxidizer, values_oxidizer);

					OpenSMOKE::OpenSMOKEVectorDouble omega_oxidizer(thermodynamicsMapXML.NumberOfSpecies());
					for(unsigned int i=0;i<names_oxidizer.size();i++)
						omega_oxidizer[thermodynamicsMapXML.IndexOfSpecies(names_oxidizer[i])] = values_oxidizer[i];

					double MW_oxidizer;
					OpenSMOKE::OpenSMOKEVectorDouble x_oxidizer(thermodynamicsMapXML.NumberOfSpecies());
					thermodynamicsMapXML.MoleFractions_From_MassFractions(x_oxidizer.GetHandle(), MW_oxidizer, omega_oxidizer.GetHandle());
					for(unsigned int i=0;i<names_oxidizer.size();i++)
						values_oxidizer[i] = x_oxidizer[thermodynamicsMapXML.IndexOfSpecies(names_oxidizer[i])];
				}
				else
				{
					names_oxidizer.resize(2);	names_oxidizer[0] = "O2";	names_oxidizer[1] = "N2";
					values_oxidizer.resize(2);	values_oxidizer[0] = 0.21;	values_oxidizer[1] = 0.79;
				}

				std::vector<double> values = thermodynamicsMapXML.GetMoleFractionsFromEquivalenceRatio(equivalence_ratio, names_fuel, values_fuel, names_oxidizer, values_oxidizer );

				const double sum =std::accumulate(values.begin(),values.end(),0.);
				if (sum<(1.-1e-6) || sum>(1.+1e-6))
					OpenSMOKE::FatalErrorMessage("The mole fractions must sum to 1.");

				OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML.NumberOfSpecies());
				for(unsigned int i=0;i<thermodynamicsMapXML.NumberOfSpecies();i++)
					x[i+1] = values[i]/sum;
				ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega, true);
				double MW;
				thermodynamicsMapXML.MassFractions_From_MoleFractions(omega.GetHandle(), MW, x.GetHandle());
			}
		}

		if (density_assigned == true)
		{
			const double MW = thermodynamicsMapXML.MolecularWeight_From_MassFractions(omega.GetHandle());
			if (temperature_assigned == true)	P_Pa	= rho*PhysicalConstants::R_J_kmol*T/MW;
			if (pressure_assigned == true)		T		= P_Pa*MW/PhysicalConstants::R_J_kmol/rho;
		}
	}

	template<typename Thermodynamics>
	void GetGasStatusFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, Thermodynamics& thermodynamicsMapXML,
		std::vector<double>& T, std::vector<double>& P_Pa, std::vector<OpenSMOKE::OpenSMOKEVectorDouble>& omega, std::vector<double>& equivalence_ratios)
	{
		Grammar_GasStatus grammar_gas_status;
		dictionary.SetGrammar(grammar_gas_status);

		unsigned int state_variables = 0;
		bool temperature_assigned = false;
		bool pressure_assigned = false;
		bool density_assigned = false;

		// Temperature
		{
			T.resize(1);
			if (dictionary.CheckOption("@Temperature") == true)
			{
				double value;
				std::string units;
				dictionary.ReadMeasure("@Temperature", value, units);

				if (units == "K")			T[0] = value;
				else if (units == "C")		T[0] = value + 273.15;
				else OpenSMOKE::FatalErrorMessage("Unknown temperature units");

				state_variables++;
				temperature_assigned = true;
			}
		}

		// Pressure
		{
			P_Pa.resize(1);
			if (dictionary.CheckOption("@Pressure") == true)
			{
				double value;
				std::string units;
				dictionary.ReadMeasure("@Pressure", value, units);

				if (units == "Pa")			P_Pa[0] = value;
				else if (units == "bar")	P_Pa[0] = value*1.e5;
				else if (units == "atm")	P_Pa[0] = value*101325.;
				else OpenSMOKE::FatalErrorMessage("Unknown pressure units");

				state_variables++;
				pressure_assigned = true;
			}
		}

		// Density
		double rho;
		{
			if (dictionary.CheckOption("@Density") == true)
			{
				double value;
				std::string units;
				dictionary.ReadMeasure("@Density", value, units);

				if (units == "kg/m3")		rho = value;
				else if (units == "g/cm3")	rho = value*1.e3;
				else OpenSMOKE::FatalErrorMessage("Unknown density units");

				state_variables++;
				density_assigned = true;
			}
		}

		if (state_variables != 2)
			OpenSMOKE::FatalErrorMessage("The status of a gas mixture requires any 2 (and only 2) among: @Temperature, @Pressure and @Density");

		// Composition
		{
			{
				std::vector<std::string> names;
				std::vector<double> values;

				if (dictionary.CheckOption("@MoleFractions") == true)
				{
					omega.resize(1);

					dictionary.ReadOption("@MoleFractions", names, values);

					const double sum = std::accumulate(values.begin(), values.end(), 0.);
					if (sum<(1. - 1e-6) || sum>(1. + 1e-6))
						OpenSMOKE::FatalErrorMessage("The mole fractions must sum to 1.");

					OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML.NumberOfSpecies());
					for (unsigned int i = 0; i<names.size(); i++)
						x[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i] / sum;
					ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega[0], true);
					double MW;
					thermodynamicsMapXML.MassFractions_From_MoleFractions(omega[0].GetHandle(), MW, x.GetHandle());
				}
				else if (dictionary.CheckOption("@MassFractions") == true)
				{
					omega.resize(1);

					dictionary.ReadOption("@MassFractions", names, values);

					const double sum = std::accumulate(values.begin(), values.end(), 0.);
					if (sum<(1. - 1e-6) || sum>(1. + 1e-6))
						OpenSMOKE::FatalErrorMessage("The mass fractions must sum to 1.");

					ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega[0], true);
					for (unsigned int i = 0; i<names.size(); i++)
						omega[0][thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i] / sum;
				}
				else if (dictionary.CheckOption("@Moles") == true)
				{
					omega.resize(1);

					dictionary.ReadOption("@Moles", names, values);
					const double sum = std::accumulate(values.begin(), values.end(), 0.);
					OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML.NumberOfSpecies());
					for (unsigned int i = 0; i<names.size(); i++)
						x[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i] / sum;
					ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega[0], true);
					double MW;
					thermodynamicsMapXML.MassFractions_From_MoleFractions(omega[0].GetHandle(), MW, x.GetHandle());
				}
				else if (dictionary.CheckOption("@Masses") == true)
				{
					omega.resize(1);
					dictionary.ReadOption("@Masses", names, values);
					const double sum = std::accumulate(values.begin(), values.end(), 0.);
					ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega[0], true);
					for (unsigned int i = 0; i<names.size(); i++)
						omega[0][thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i] / sum;
				}
			}

			if (dictionary.CheckOption("@EquivalenceRatio") == true)
			{
				dictionary.ReadOption("@EquivalenceRatio", equivalence_ratios);

				omega.resize(equivalence_ratios.size());

				std::vector<std::string> names_fuel;
				std::vector<double> values_fuel;

				if (dictionary.CheckOption("@FuelMoleFractions") == true)
				{
					dictionary.ReadOption("@FuelMoleFractions", names_fuel, values_fuel);

					const double sum = std::accumulate(values_fuel.begin(), values_fuel.end(), 0.);
					if (sum<(1. - 1e-6) || sum>(1. + 1e-6))
						OpenSMOKE::FatalErrorMessage("The fuel mass fractions must sum to 1.");

					for (unsigned int i = 0; i<names_fuel.size(); i++)
						values_fuel[i] /= sum;
				}
				else if (dictionary.CheckOption("@FuelMassFractions") == true)
				{
					dictionary.ReadOption("@FuelMassFractions", names_fuel, values_fuel);

					const double sum = std::accumulate(values_fuel.begin(), values_fuel.end(), 0.);
					if (sum<(1. - 1e-6) || sum>(1. + 1e-6))
						OpenSMOKE::FatalErrorMessage("The fuel mass fractions must sum to 1.");

					OpenSMOKE::OpenSMOKEVectorDouble omega_fuel(thermodynamicsMapXML.NumberOfSpecies());
					for (unsigned int i = 0; i<names_fuel.size(); i++)
						omega_fuel[thermodynamicsMapXML.IndexOfSpecies(names_fuel[i])] = values_fuel[i] / sum;

					double MW_fuel;
					OpenSMOKE::OpenSMOKEVectorDouble x_fuel(thermodynamicsMapXML.NumberOfSpecies());
					thermodynamicsMapXML.MoleFractions_From_MassFractions(x_fuel.GetHandle(), MW_fuel, omega_fuel.GetHandle());
					for (unsigned int i = 0; i<names_fuel.size(); i++)
						values_fuel[i] = x_fuel[thermodynamicsMapXML.IndexOfSpecies(names_fuel[i])];
				}
				else if (dictionary.CheckOption("@FuelMoles") == true)
				{
					dictionary.ReadOption("@FuelMoles", names_fuel, values_fuel);
				}
				else if (dictionary.CheckOption("@FuelMasses") == true)
				{
					dictionary.ReadOption("@FuelMasses", names_fuel, values_fuel);

					OpenSMOKE::OpenSMOKEVectorDouble omega_fuel(thermodynamicsMapXML.NumberOfSpecies());
					for (unsigned int i = 0; i<names_fuel.size(); i++)
						omega_fuel[thermodynamicsMapXML.IndexOfSpecies(names_fuel[i])] = values_fuel[i];

					double MW_fuel;
					OpenSMOKE::OpenSMOKEVectorDouble x_fuel(thermodynamicsMapXML.NumberOfSpecies());
					thermodynamicsMapXML.MoleFractions_From_MassFractions(x_fuel.GetHandle(), MW_fuel, omega_fuel.GetHandle());
					for (unsigned int i = 0; i<names_fuel.size(); i++)
						values_fuel[i] = x_fuel[thermodynamicsMapXML.IndexOfSpecies(names_fuel[i])];
				}
				else
				{
					OpenSMOKE::FatalErrorMessage("The @EquivalenceRatio option requires the user specifies the fuel composition");
				}

				std::vector<std::string> names_oxidizer;
				std::vector<double> values_oxidizer;

				if (dictionary.CheckOption("@OxidizerMoleFractions") == true)
				{
					dictionary.ReadOption("@OxidizerMoleFractions", names_oxidizer, values_oxidizer);

					const double sum = std::accumulate(values_oxidizer.begin(), values_oxidizer.end(), 0.);
					if (sum<(1. - 1e-6) || sum>(1. + 1e-6))
						OpenSMOKE::FatalErrorMessage("The oxidizer mass fractions must sum to 1.");

					for (unsigned int i = 0; i<names_oxidizer.size(); i++)
						values_oxidizer[i] /= sum;
				}
				else if (dictionary.CheckOption("@OxidizerMassFractions") == true)
				{
					dictionary.ReadOption("@OxidizerMassFractions", names_oxidizer, values_oxidizer);

					const double sum = std::accumulate(values_oxidizer.begin(), values_oxidizer.end(), 0.);
					if (sum<(1. - 1e-6) || sum>(1. + 1e-6))
						OpenSMOKE::FatalErrorMessage("The oxidizer mass fractions must sum to 1.");

					OpenSMOKE::OpenSMOKEVectorDouble omega_oxidizer(thermodynamicsMapXML.NumberOfSpecies());
					for (unsigned int i = 0; i<names_oxidizer.size(); i++)
						omega_oxidizer[thermodynamicsMapXML.IndexOfSpecies(names_oxidizer[i])] = values_oxidizer[i] / sum;

					double MW_oxidizer;
					OpenSMOKE::OpenSMOKEVectorDouble x_oxidizer(thermodynamicsMapXML.NumberOfSpecies());
					thermodynamicsMapXML.MoleFractions_From_MassFractions(x_oxidizer.GetHandle(), MW_oxidizer, omega_oxidizer.GetHandle());
					for (unsigned int i = 0; i<names_oxidizer.size(); i++)
						values_oxidizer[i] = x_oxidizer[thermodynamicsMapXML.IndexOfSpecies(names_oxidizer[i])];
				}
				else if (dictionary.CheckOption("@OxidizerMoles") == true)
				{
					dictionary.ReadOption("@OxidizerMoles", names_oxidizer, values_oxidizer);
				}
				else if (dictionary.CheckOption("@OxidizerMasses") == true)
				{
					dictionary.ReadOption("@OxidizerMasses", names_oxidizer, values_oxidizer);

					OpenSMOKE::OpenSMOKEVectorDouble omega_oxidizer(thermodynamicsMapXML.NumberOfSpecies());
					for (unsigned int i = 0; i<names_oxidizer.size(); i++)
						omega_oxidizer[thermodynamicsMapXML.IndexOfSpecies(names_oxidizer[i])] = values_oxidizer[i];

					double MW_oxidizer;
					OpenSMOKE::OpenSMOKEVectorDouble x_oxidizer(thermodynamicsMapXML.NumberOfSpecies());
					thermodynamicsMapXML.MoleFractions_From_MassFractions(x_oxidizer.GetHandle(), MW_oxidizer, omega_oxidizer.GetHandle());
					for (unsigned int i = 0; i<names_oxidizer.size(); i++)
						values_oxidizer[i] = x_oxidizer[thermodynamicsMapXML.IndexOfSpecies(names_oxidizer[i])];
				}
				else
				{
					names_oxidizer.resize(2);	names_oxidizer[0] = "O2";	names_oxidizer[1] = "N2";
					values_oxidizer.resize(2);	values_oxidizer[0] = 0.21;	values_oxidizer[1] = 0.79;
				}

				for (unsigned int k = 0; k < equivalence_ratios.size(); k++)
				{
					std::vector<double> values = thermodynamicsMapXML.GetMoleFractionsFromEquivalenceRatio(equivalence_ratios[k], names_fuel, values_fuel, names_oxidizer, values_oxidizer);

					const double sum = std::accumulate(values.begin(), values.end(), 0.);
					if (sum<(1. - 1e-6) || sum>(1. + 1e-6))
						OpenSMOKE::FatalErrorMessage("The mole fractions must sum to 1.");

					OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML.NumberOfSpecies());
					for (unsigned int i = 0; i<thermodynamicsMapXML.NumberOfSpecies(); i++)
						x[i + 1] = values[i] / sum;
					ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega[k], true);
					double MW;
					thermodynamicsMapXML.MassFractions_From_MoleFractions(omega[k].GetHandle(), MW, x.GetHandle());
				}
			}
		}

		if (density_assigned == true)
		{
			if (omega.size() != 1)
				OpenSMOKE::FatalErrorMessage("The multiple @EquivalenceRatio option cannot be used together with the @Density option.");
			else
			{
				const double MW = thermodynamicsMapXML.MolecularWeight_From_MassFractions(omega[0].GetHandle());
				if (temperature_assigned == true)	P_Pa[0] = rho*PhysicalConstants::R_J_kmol*T[0] / MW;
				if (pressure_assigned == true)		T[0] = P_Pa[0] * MW / PhysicalConstants::R_J_kmol / rho;
			}
		}

		// TOADJUST
		{
			const double T_ = T[0];
			const double P_Pa_ = P_Pa[0];
			T.resize(omega.size());
			P_Pa.resize(omega.size());
			for (unsigned int k = 0; k < omega.size(); k++)
			{
				T[k] = T_;
				P_Pa[k] = P_Pa_;
			}
		}
	}

	template<typename Thermodynamics>
	void GetGasStatusFromDictionary(	OpenSMOKE::OpenSMOKE_Dictionary& dictionary, Thermodynamics& thermodynamicsMapXML,
										double& T, double& P_Pa, Eigen::VectorXd& omega)
	{
		OpenSMOKE::OpenSMOKEVectorDouble temp(thermodynamicsMapXML.NumberOfSpecies());
		GetGasStatusFromDictionary(dictionary, thermodynamicsMapXML, T, P_Pa, temp);
		omega.resize(thermodynamicsMapXML.NumberOfSpecies());
		for (unsigned int i = 0; i < omega.size(); i++)
			omega(i) = temp[i + 1];
	}

	template<typename Thermodynamics>
	void GetGasStatusFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, Thermodynamics& thermodynamicsMapXML,
					std::vector<std::string>& names_fuel, std::vector<double>& molefractions_fuel,
					std::vector<std::string>& names_oxidizer, std::vector<double>& molefractions_oxidizer)
	{
		Grammar_GasStatus grammar_gas_status;
		dictionary.SetGrammar(grammar_gas_status);

		if (dictionary.CheckOption("@EquivalenceRatio") == false)
			OpenSMOKE::FatalErrorMessage("The @EquivalenceRatio option is mandatory");

		if (dictionary.CheckOption("@FuelMoleFractions") == true)
		{
			dictionary.ReadOption("@FuelMoleFractions", names_fuel, molefractions_fuel);

			const double sum =std::accumulate(molefractions_fuel.begin(),molefractions_fuel.end(),0.);
			if (sum<(1.-1e-6) || sum>(1.+1e-6))
				OpenSMOKE::FatalErrorMessage("The fuel mass fractions must sum to 1.");

			for(unsigned int i=0;i<names_fuel.size();i++)
				molefractions_fuel[i] /= sum;
		}
		else if (dictionary.CheckOption("@FuelMassFractions") == true)
		{
			dictionary.ReadOption("@FuelMassFractions", names_fuel, molefractions_fuel);

			const double sum =std::accumulate(molefractions_fuel.begin(),molefractions_fuel.end(),0.);
			if (sum<(1.-1e-6) || sum>(1.+1e-6))
				OpenSMOKE::FatalErrorMessage("The fuel mass fractions must sum to 1.");

			OpenSMOKE::OpenSMOKEVectorDouble omega_fuel(thermodynamicsMapXML.NumberOfSpecies());
			for(unsigned int i=0;i<names_fuel.size();i++)
				omega_fuel[thermodynamicsMapXML.IndexOfSpecies(names_fuel[i])] = molefractions_fuel[i]/sum;

			double MW_fuel;
			OpenSMOKE::OpenSMOKEVectorDouble x_fuel(thermodynamicsMapXML.NumberOfSpecies());
			thermodynamicsMapXML.MoleFractions_From_MassFractions(x_fuel.GetHandle(), MW_fuel, omega_fuel.GetHandle());
			for(unsigned int i=0;i<names_fuel.size();i++)
				molefractions_fuel[i] = x_fuel[thermodynamicsMapXML.IndexOfSpecies(names_fuel[i])];
		}
		else if (dictionary.CheckOption("@FuelMoles") == true)
		{
			dictionary.ReadOption("@FuelMoles", names_fuel, molefractions_fuel);
		}
		else if (dictionary.CheckOption("@FuelMasses") == true)
		{
			dictionary.ReadOption("@FuelMasses", names_fuel, molefractions_fuel);

			OpenSMOKE::OpenSMOKEVectorDouble omega_fuel(thermodynamicsMapXML.NumberOfSpecies());
			for(unsigned int i=0;i<names_fuel.size();i++)
				omega_fuel[thermodynamicsMapXML.IndexOfSpecies(names_fuel[i])] = molefractions_fuel[i];

			double MW_fuel;
			OpenSMOKE::OpenSMOKEVectorDouble x_fuel(thermodynamicsMapXML.NumberOfSpecies());
			thermodynamicsMapXML.MoleFractions_From_MassFractions(x_fuel.GetHandle(), MW_fuel, omega_fuel.GetHandle());
			for(unsigned int i=0;i<names_fuel.size();i++)
				molefractions_fuel[i] = x_fuel[thermodynamicsMapXML.IndexOfSpecies(names_fuel[i])];
		}
		else
		{
			OpenSMOKE::FatalErrorMessage("The @EquivalenceRatio option requires the user specifies the fuel composition");
		}



		if (dictionary.CheckOption("@OxidizerMoleFractions") == true)
		{
			dictionary.ReadOption("@OxidizerMoleFractions", names_oxidizer, molefractions_oxidizer);

			const double sum =std::accumulate(molefractions_oxidizer.begin(),molefractions_oxidizer.end(),0.);
			if (sum<(1.-1e-6) || sum>(1.+1e-6))
				OpenSMOKE::FatalErrorMessage("The oxidizer mass fractions must sum to 1.");

			for(unsigned int i=0;i<names_oxidizer.size();i++)
				molefractions_oxidizer[i] /= sum;
		}
		else if (dictionary.CheckOption("@OxidizerMassFractions") == true)
		{
			dictionary.ReadOption("@OxidizerMassFractions", names_oxidizer, molefractions_oxidizer);

			const double sum =std::accumulate(molefractions_oxidizer.begin(),molefractions_oxidizer.end(),0.);
			if (sum<(1.-1e-6) || sum>(1.+1e-6))
				OpenSMOKE::FatalErrorMessage("The oxidizer mass fractions must sum to 1.");

			OpenSMOKE::OpenSMOKEVectorDouble omega_oxidizer(thermodynamicsMapXML.NumberOfSpecies());
			for(unsigned int i=0;i<names_oxidizer.size();i++)
				omega_oxidizer[thermodynamicsMapXML.IndexOfSpecies(names_oxidizer[i])] = molefractions_oxidizer[i]/sum;

			double MW_oxidizer;
			OpenSMOKE::OpenSMOKEVectorDouble x_oxidizer(thermodynamicsMapXML.NumberOfSpecies());
			thermodynamicsMapXML.MoleFractions_From_MassFractions(x_oxidizer.GetHandle(), MW_oxidizer, omega_oxidizer.GetHandle());
			for(unsigned int i=0;i<names_oxidizer.size();i++)
				molefractions_oxidizer[i] = x_oxidizer[thermodynamicsMapXML.IndexOfSpecies(names_oxidizer[i])];
		}
		else if (dictionary.CheckOption("@OxidizerMoles") == true)
		{
			dictionary.ReadOption("@OxidizerMoles", names_oxidizer, molefractions_oxidizer);
		}
		else if (dictionary.CheckOption("@OxidizerMasses") == true)
		{
			dictionary.ReadOption("@OxidizerMasses", names_oxidizer, molefractions_oxidizer);

			OpenSMOKE::OpenSMOKEVectorDouble omega_oxidizer(thermodynamicsMapXML.NumberOfSpecies());
			for(unsigned int i=0;i<names_oxidizer.size();i++)
				omega_oxidizer[thermodynamicsMapXML.IndexOfSpecies(names_oxidizer[i])] = molefractions_oxidizer[i];

			double MW_oxidizer;
			OpenSMOKE::OpenSMOKEVectorDouble x_oxidizer(thermodynamicsMapXML.NumberOfSpecies());
			thermodynamicsMapXML.MoleFractions_From_MassFractions(x_oxidizer.GetHandle(), MW_oxidizer, omega_oxidizer.GetHandle());
			for(unsigned int i=0;i<names_oxidizer.size();i++)
				molefractions_oxidizer[i] = x_oxidizer[thermodynamicsMapXML.IndexOfSpecies(names_oxidizer[i])];
		}
		else
		{
			names_oxidizer.resize(2);		names_oxidizer[0] = "O2";		names_oxidizer[1] = "N2";
			molefractions_oxidizer.resize(2);	molefractions_oxidizer[0] = 0.21;	molefractions_oxidizer[1] = 0.79;
		}
	}

}

#endif	/* OpenSMOKE_Grammar_GasStatus_Options_H */

