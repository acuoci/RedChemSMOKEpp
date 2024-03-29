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

#include "BatchReactor_OdeInterfaces.h"
#include "math/PhysicalConstants.h"

namespace OpenSMOKE
{
	#if OPENSMOKE_USE_BZZMATH == 1
	ODESystem_BzzOde_BatchReactor* ptOde_BatchReactor_NonIsothermal_UserDefinedVolume;
	void ODE_Print_BatchReactor_NonIsothermal_UserDefinedVolume(BzzVector &Y, double t)
	{
		ptOde_BatchReactor_NonIsothermal_UserDefinedVolume->MyPrint(Y,t);
	}
	#endif

	BatchReactor_NonIsothermal_UserDefinedVolume::BatchReactor_NonIsothermal_UserDefinedVolume(	
								OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
								OpenSMOKE::ODE_Parameters& ode_parameters,
								OpenSMOKE::BatchReactor_Options& batch_options,
								OpenSMOKE::OnTheFlyROPA& on_the_fly_ropa,
								OpenSMOKE::OnTheFlyCEMA& on_the_fly_cema,
								OpenSMOKE::OnTheFlyPostProcessing& on_the_fly_post_processing,
								OpenSMOKE::IgnitionDelayTimes_Analyzer& idts_analyzer,
								OpenSMOKE::PolimiSoot_Analyzer& polimi_soot_analyzer,
								const double V0, const double T0, const double P0,
								const OpenSMOKE::OpenSMOKEVectorDouble& omega0,
                                const double tStart, 
                                const double global_thermal_exchange_coefficient,
								const double exchange_area, const double T_environment):

		BatchReactor(thermodynamicsMap, kineticsMap, ode_parameters, batch_options, on_the_fly_ropa, on_the_fly_cema, on_the_fly_post_processing, idts_analyzer, polimi_soot_analyzer)
	
	{
		type_ = BATCH_REACTOR_NONISOTHERMAL_CONSTANTV;

		iteration_ = 0;
		counter_file_video_ = 0;
		counter_file_ASCII_ = 0;
		counter_file_XML_ = 0;
		counter_sensitivity_XML_ = 0;

		V0_ = V0;
		T0_ = T0;
		P0_ = P0;
		omega0_ = omega0;

		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();
		NE_ = NC_+1;

		MemoryAllocation();

		MW0_ = thermodynamicsMap_.MolecularWeight_From_MassFractions(omega0_.GetHandle());
		rho0_ = P0_ * MW0_ / (PhysicalConstants::R_J_kmol*T0_);
		mass_ = rho0_ * V0_;

		T_		= T0_;
		P_		= P0_;
		V_      = V0_;
		omega_	= omega0_;
		MW_     = MW0_;
		rho_    = rho0_;
		QR_     = 0.;
		
		VOld_ = V0_;
		POld_ = P0_;
		tOld_ = tStart;
                
        global_thermal_exchange_coefficient_ = global_thermal_exchange_coefficient;
		exchange_area_ = exchange_area;
		T_environment_ = T_environment;   

		thermodynamicsMap_.MoleFractions_From_MassFractions(x0_.GetHandle(), MW0_, omega0_.GetHandle());
		thermodynamicsMap_.SetPressure(P0_);
		thermodynamicsMap_.SetTemperature(T0_);
		H0_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x0_.GetHandle());
		U0_ = thermodynamicsMap_.uMolar_Mixture_From_MoleFractions(x0_.GetHandle());

		OpenAllFiles();
	}
	
	void BatchReactor_NonIsothermal_UserDefinedVolume::SetVolumeProfile(OpenSMOKE::BatchReactor_VolumeProfile&  profile)
	{
		volume_profile_ = &profile;
	}

	int BatchReactor_NonIsothermal_UserDefinedVolume::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
	{
		// Recover mass fractions
		#ifdef CHECK_MASSFRACTIONS
		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = std::max(y[i], 0.);
		#else
		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = y[i];
		#endif

		// Recover temperature
		T_ = y[NC_ + 1];

		// Calculates the mole fractions and the molecular weight
		thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());

		// Update the volume
		double dV_over_dt_ = 0.;
		if (volume_profile_->type() == BatchReactor_VolumeProfile::BatchReactor_VolumeProfile_PressureCoefficient)
		{
			// Updates the pressure
			P_ = P0_ + volume_profile_->pressureCoefficient()*t;

			// Calculates thermodynamic properties
			thermodynamicsMap_.SetTemperature(T_);
			thermodynamicsMap_.SetPressure(P_);
			const double CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle());
			CvMixMass_ = (CpMixMolar - PhysicalConstants::R_J_kmol) / MW_;
			CpMixMass_ = CpMixMolar / MW_;
			const double gamma = CpMixMass_ / CvMixMass_;	// isoentropic

			// Updates the volume
			V_ = V0_*std::pow(P0_ / P_, 1. / gamma);
			dV_over_dt_ = -1./gamma*V_/P_*volume_profile_->pressureCoefficient();
		}
		else if (volume_profile_->type() == BatchReactor_VolumeProfile::BatchReactor_VolumeProfile_VolumeHistory)
		{
			V_ = volume_profile_->operator()(t);
			dV_over_dt_ = (t>tOld_) ? (V_-VOld_)/(t-tOld_) : 0.;
		}
		else if (volume_profile_->type() == BatchReactor_VolumeProfile::BatchReactor_VolumeProfile_FromPressureProfile)
		{
			// Extract pressure from the user defined profile
			P_ = volume_profile_->operator()(t);

			// Calculates thermodynamic properties
			thermodynamicsMap_.SetTemperature(T_);
			thermodynamicsMap_.SetPressure(P_);
			const double CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle());
			CvMixMass_ = (CpMixMolar - PhysicalConstants::R_J_kmol) / MW_;
			CpMixMass_ = CpMixMolar / MW_;
			const double gamma = CpMixMass_/CvMixMass_;

			// Updates the volume	
			V_ = V0_*std::pow(P0_/P_, 1./gamma);
			dV_over_dt_ = (t>tOld_) ? -1./gamma*V_/P_*(P_-POld_)/(t-tOld_) : 0.;
		}

		// Updates the density
		rho_ = mass_ / V_;

		// Updates the concentrations of species
		cTot_ = rho_ / MW_;
		Product(cTot_, x_, &c_);
					
		// Update the pressure profile
		if (volume_profile_->type() == BatchReactor_VolumeProfile::BatchReactor_VolumeProfile_VolumeHistory)
			P_ = cTot_ * PhysicalConstants::R_J_kmol * T_;
			 
		// Calculates thermodynamic properties
		thermodynamicsMap_.SetTemperature(T_);
		thermodynamicsMap_.SetPressure(P_);
		const double CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle());
		CvMixMass_ = (CpMixMolar - PhysicalConstants::R_J_kmol) / MW_;
		CpMixMass_ = CpMixMolar / MW_;

		// Calculates kinetics
		kineticsMap_.SetTemperature(T_);
		kineticsMap_.SetPressure(P_);	
		kineticsMap_.ReactionRates(c_.GetHandle());
		kineticsMap_.FormationRates(R_.GetHandle());
		QR_ = kineticsMap_.HeatRelease(R_.GetHandle());

		// Number of moles
		const double sumMoleFormationRates = R_.SumElements();

		// Mass fraction equations
		for (unsigned int i=1;i<=NC_;++i)	
			dy[i] = thermodynamicsMap_.MW(i-1)*R_[i]/rho_;
						
		// Energy equation
		const double Qexchange = global_thermal_exchange_coefficient_*exchange_area_*(T_environment_-T_); // [W]
		dy[NC_+1]  = (QR_ + PhysicalConstants::R_J_kmol*T_*sumMoleFormationRates - P_*dV_over_dt_/V_ + Qexchange/V_) / 
						(rho_*CvMixMass_);

		return 0;
	}

	int BatchReactor_NonIsothermal_UserDefinedVolume::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		iteration_++;

		if (batch_options_.verbose_video() == true)
		{
			// Video output
			if (iteration_%batch_options_.n_step_video() == 1 || batch_options_.n_step_video() == 1)
			{
				counter_file_video_++;
				if (counter_file_video_ % 100 == 1)
				{
					std::cout << std::endl;
					std::cout << std::setw(6) << std::left << "#Step";
					std::cout << std::setw(16) << std::left << "Time[s]";
					std::cout << std::setw(10) << std::left << "T[K]";
					std::cout << std::setw(10) << std::left << "P[atm]";
					std::cout << std::endl;
				}
				std::cout << std::fixed << std::setw(6) << std::left << iteration_;
				std::cout << std::scientific << std::setw(16) << std::setprecision(6) << std::left << t;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(3) << T_;
				std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(3) << P_ / 101325.;
				std::cout << std::endl;
			}

			if (batch_options_.verbose_output() == true)
			{
					// ASCII file output
				if (batch_options_.verbose_ascii_file() == true)
				{
					if (iteration_%batch_options_.n_step_file() == 1 || batch_options_.n_step_file() == 1)
					{
						counter_file_ASCII_++;
						PrintFinalStatus(fASCII_, t);

						if (on_the_fly_post_processing_.is_active() == true)
							on_the_fly_post_processing_.WriteOnFile(t, 0., 0., 0., T_, P_, omega_);
					}
				}

				// XML file output
				if (batch_options_.verbose_xml_file() == true)
				{
					if (iteration_%batch_options_.n_step_file() == 1 || batch_options_.n_step_file() == 1)
					{
						counter_file_XML_++;
						fXML_ << t << " ";
						fXML_ << T_ << " ";
						fXML_ << P_ << " ";
						fXML_ << MW_ << " ";
						fXML_ << rho_ << " ";
						fXML_ << QR_ << " ";
						for (unsigned int i = 1; i <= NC_; i++)
							fXML_ << std::setprecision(12) << omega_[i] << " ";
						fXML_ << std::endl;

						// Write formation rates and reaction rates
						if (on_the_fly_post_processing_.is_active() == true)
						{
							// Write formation rates (kg/m3/s)
							for (unsigned int j = 1; j <= thermodynamicsMap_.NumberOfSpecies(); j++)
								fXML_formation_rates_ << std::scientific << std::setprecision(9) << thermodynamicsMap_.MW(j - 1) * R_[j] << " ";
							fXML_formation_rates_ << std::endl;

							// Write reaction rates (kmol/m3/s)
							std::vector<double> r = kineticsMap_.GiveMeReactionRates();
							for (unsigned int j = 0; j < r.size(); j++)
								fXML_reaction_rates_ << std::scientific << std::setprecision(9) << r[j] << " ";
							fXML_reaction_rates_ << std::endl;
						}
					}
				}

				// Rate of Production Analysis (on the fly)
				if (on_the_fly_ropa_.is_active() == true)
					on_the_fly_ropa_.Analyze(fROPA_, iteration_, t, T_, P_, c_, omega_, omega0_);

				// Polimi Soot Analyzer (on the fly)
				if (polimi_soot_analyzer_.is_active() == true)
					PolimiSootAnalysis(t);
			}
		}

		// Ignition delay times (on the fly)
		if (idts_analyzer_.is_active() == true)
			idts_analyzer_.Analyze(t, T_, P_, x_.GetHandle());

		// Sensitivity analysis (on the fly)
		if (batch_options_.sensitivity_analysis() == true)
			SensitivityAnalysis(t, y);
			
		// Saving old values
		{
			tOld_ = t;
			if (volume_profile_->type() == BatchReactor_VolumeProfile::BatchReactor_VolumeProfile_PressureCoefficient)
			{
				VOld_ = V0_ * P0_/(P0_+volume_profile_->pressureCoefficient()*t);
			}
			else if (volume_profile_->type() == BatchReactor_VolumeProfile::BatchReactor_VolumeProfile_VolumeHistory)
			{
				VOld_ = volume_profile_->operator()(t);
			}
			else if (volume_profile_->type() == BatchReactor_VolumeProfile::BatchReactor_VolumeProfile_FromPressureProfile)
			{
				POld_ = volume_profile_->operator()(t);
			}
		}

		return 0;
	}

	void BatchReactor_NonIsothermal_UserDefinedVolume::Solve(const double t0, const double tf)
	{
		if (batch_options_.verbose_video() == true)
		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
			std::cout << " Solving the batch reactor...                                                " << std::endl;
			std::cout << "-----------------------------------------------------------------------------" << std::endl;
		}

		for(unsigned int i=1;i<=NC_;++i)
			y0_[i] = omega0_[i];
		y0_[NC_+1] = T0_;

		// Print intial conditions
		{
			OpenSMOKE::OpenSMOKEVectorDouble dy0(y0_.Size());
			Equations(t0, y0_, dy0);
			Print(t0, y0_);
		}

		if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			// Min and max values
			Eigen::VectorXd yMin(NE_); for (unsigned int i = 0; i < NE_; i++) yMin(i) = 0.;  yMin(NC_) = 0.;
			Eigen::VectorXd yMax(NE_); for (unsigned int i = 0; i < NE_; i++) yMax(i) = 1.;  yMax(NC_) = 10000.;

			// Initial conditions
			Eigen::VectorXd y0_eigen(y0_.Size());
			y0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Create the solver
			typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_OpenSMOKE_BatchReactor> denseOde;
			typedef OdeSMOKE::MethodGear<denseOde> methodGear;
			OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
			ode_solver.SetReactor(this);

			// Set initial conditions
			ode_solver.SetInitialConditions(t0, y0_eigen);

			// Set linear algebra options
			ode_solver.SetLinearAlgebraSolver(ode_parameters_.linear_algebra());
			ode_solver.SetFullPivoting(ode_parameters_.full_pivoting());

			// Set relative and absolute tolerances
			ode_solver.SetAbsoluteTolerances(ode_parameters_.absolute_tolerance());
			ode_solver.SetRelativeTolerances(ode_parameters_.relative_tolerance());

			// Set minimum and maximum values
			ode_solver.SetMinimumValues(yMin);
			ode_solver.SetMaximumValues(yMax);

			// Set maximum number of steps
			if (ode_parameters_.maximum_number_of_steps() > 0)
				ode_solver.SetMaximumNumberOfSteps(ode_parameters_.maximum_number_of_steps());

			// Set maximum integration order
			if (ode_parameters_.maximum_order() > 0)
				ode_solver.SetMaximumOrder(ode_parameters_.maximum_order());

			// Set maximum step size allowed
			if (ode_parameters_.maximum_step() > 0)
				ode_solver.SetMaximumStepSize(ode_parameters_.maximum_step());

			// Set minimum step size allowed
			if (ode_parameters_.minimum_step() > 0)
				ode_solver.SetMinimumStepSize(ode_parameters_.minimum_step());

			// Set initial step size
			if (ode_parameters_.initial_step() > 0)
				ode_solver.SetFirstStepSize(ode_parameters_.initial_step());

			// Solve the system
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			OdeSMOKE::OdeStatus status = ode_solver.Solve(tf);
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			// Check the solution
			if (status > 0)
			{
				ode_solver.Solution(yf_eigen);
				yf_.CopyFrom(yf_eigen.data());
				ode_parameters_.TransferDataFromOdeSolver(ode_solver, tEnd - tStart);
			}
		}
		#if OPENSMOKE_USE_BZZMATH == 1
		else if (ode_parameters_.type() == ODE_Parameters::ODE_INTEGRATOR_BZZODE)
		{
			// Min and max values
			BzzVector yMin(NE_); yMin=0.; yMin[NC_+1] = 0.;
			BzzVector yMax(NE_); yMax=1.; yMax[NC_+1] = 10000.;
			
			// Initial conditions
			BzzVector y0_bzz(y0_.Size());
			y0_.CopyTo(y0_bzz.GetHandle());
			
			// Final solution
			BzzVector yf_bzz(y0_bzz.Size());

			ODESystem_BzzOde_BatchReactor odebatch(*this);
			BzzOdeStiffObject o(y0_bzz, t0, &odebatch);
			
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			o.SetTolAbs(ode_parameters_.absolute_tolerance());
			o.SetTolRel(ode_parameters_.relative_tolerance());

			ptOde_BatchReactor_NonIsothermal_UserDefinedVolume = &odebatch;
			o.StepPrint(ODE_Print_BatchReactor_NonIsothermal_UserDefinedVolume);

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			yf_bzz = o(tf,tf);
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			yf_.CopyFrom(yf_bzz.GetHandle());

			ode_parameters_.SetCPUTime(tEnd-tStart);
			ode_parameters_.SetNumberOfFunctionCalls(o.GetNumFunction());
			ode_parameters_.SetNumberOfJacobians(o.GetNumNumericalJacobian());
			ode_parameters_.SetNumberOfFactorizations(o.GetNumFactorization());
			ode_parameters_.SetNumberOfSteps(o.GetNumStep());
			ode_parameters_.SetLastOrderUsed(o.GetOrderUsed());
			ode_parameters_.SetLastStepUsed(o.GetHUsed());
		}
		#endif
		else 
		{
			SolveOpenSourceSolvers(t0, tf);
		}

		if (batch_options_.verbose_video() == true)
			ode_parameters_.Status(std::cout);

		FinalStatus(t0, tf);
		FinalSummary(batch_options_.output_path() / "FinalSummary.out", t0, tf);

		if (idts_analyzer_.is_active() == true)
			idts_analyzer_.PrintOnFile(batch_options_.output_path() / "IDT.out");

		CloseAllFiles();
	}

	void BatchReactor_NonIsothermal_UserDefinedVolume::SensitivityAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		if (iteration_ == 1)
		{
			// Writes the coefficients on file (only on request)
			if (iteration_%batch_options_.n_step_file() == 1 || batch_options_.n_step_file() == 1)		
			{
				counter_sensitivity_XML_++;
				for (unsigned int k=0;k<indices_of_sensitivity_species_.size();k++)
				{
					for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
						fSensitivityChildXML_[k] << 0. << " ";		
					fSensitivityChildXML_[k] << std::endl;	
				}

				for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
					fSensitivityChildXML_[indices_of_sensitivity_species_.size()] << 0. << " ";
				fSensitivityChildXML_[indices_of_sensitivity_species_.size()] << std::endl;
			}
		}
		else
		{
			// Scaling factors
			for(unsigned int j=1;j<=NC_;j++)
				scaling_Jp_[j] = thermodynamicsMap_.MW(j-1)/rho_;
			scaling_Jp_[NC_+1] = 1./(rho_*CvMixMass_);

			// Calculates the current Jacobian
			NumericalJacobian(t, y, Jnum_);

			// Recover variables
			{
				#ifdef CHECK_MASSFRACTIONS
				for(unsigned int i=1;i<=NC_;++i)
					omega_[i] = std::max(y[i], 0.);
				#else
				for(unsigned int i=1;i<=NC_;++i)
					omega_[i] = y[i];
				#endif

				// Recover temperature
				T_ = y[NC_+1];

				// Calculates the pressure and the concentrations of species
				thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
				cTot_ = rho_/MW_;
				Product(cTot_, x_, &c_);
				P_ = cTot_ * PhysicalConstants::R_J_kmol * T_;
			}

			// Calculates the current sensitivity coefficients
			sensitivityMap_->Calculate(t, T_, P_, c_, Jnum_, scaling_Jp_);

			// Writes the coefficients on file (only on request)
			if (iteration_%batch_options_.n_step_file() == 1 || batch_options_.n_step_file() == 1)		
			{
				counter_sensitivity_XML_++;
				for (unsigned int k=0;k<indices_of_sensitivity_species_.size();k++)
				{
					const unsigned int i = indices_of_sensitivity_species_[k];
					for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
					{
						double sum = 0.;
						for (unsigned int kk=1;kk<=NC_;kk++)
							sum += sensitivityMap_->sensitivity_coefficients()(kk-1,j-1)/thermodynamicsMap_.MW(kk-1);
						sum *= x_[i]*MW_;

						double coefficient = sensitivityMap_->sensitivity_coefficients()(i-1,j-1)*MW_/thermodynamicsMap_.MW(i-1) - sum;
						fSensitivityChildXML_[k] << coefficient << " ";
					}		
					fSensitivityChildXML_[k] << std::endl;	
				}

				for (unsigned int j=1;j<=sensitivityMap_->number_of_parameters();j++)
					fSensitivityChildXML_[indices_of_sensitivity_species_.size()] << sensitivityMap_->sensitivity_coefficients()(NC_,j-1) << " ";
				fSensitivityChildXML_[indices_of_sensitivity_species_.size()] << std::endl;
			}
		}
	}

	void BatchReactor_NonIsothermal_UserDefinedVolume::ChemicalExplosiveModeAnalysis(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		OpenSMOKE::ErrorMessage("BatchReactor_NonIsothermal_UserDefinedVolume", "ChemicalExplosiveModeAnalysis is not yet available for Adiabtic User Defined Batch Reactor");
	}

	void BatchReactor_NonIsothermal_UserDefinedVolume::SparseAnalyticalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::SparseMatrix<double> &J)
	{
		OpenSMOKE::ErrorMessage("BatchReactor_NonIsothermal_UserDefinedVolume", "SparseAnalyticalJacobian is not yet available for Adiabtic User Defined Batch Reactor");
	}

	void BatchReactor_NonIsothermal_UserDefinedVolume::DenseAnalyticalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::MatrixXd &J)
	{
		OpenSMOKE::ErrorMessage("BatchReactor_NonIsothermal_UserDefinedVolume::", "SparseAnalyticalJacobian is not yet available for Non Isothermal User Defined Volume Batch Reactor");
	}
}
