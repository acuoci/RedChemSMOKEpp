/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alessandro Stagni <alessandro.stagni@polimi.it>               |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                           |
|                                                                         |
|   Copyright(C) 2016-2012  Alberto Cuoci                                 |
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

#ifndef KINETICSGRAPH_H
#define KINETICSGRAPH_H

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// OpenSMOKE++ Dictionaries
#include "dictionary/OpenSMOKE_Dictionary"

//Boost
#include <boost/graph/adjacency_list.hpp>
typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::directedS,
		boost::no_property, boost::property < boost::edge_weight_t, double > > DirectedGraph;
typedef std::pair<int, int> Edge;

class KineticsGraph
{
public:
  void Setup(	OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
				OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML);
  
  //Copy cosntructor
  //KineticsGraph(const KineticsGraph& rhs);
  
  ~KineticsGraph();
  
  inline const DirectedGraph& DirectedKineticsGraph() const                             {return DirectedKineticsGraph_;};
  inline const Eigen::SparseMatrix<double>& stoichiometric_matrix_reactants() const     {return stoichiometric_matrix_reactants_;};
  inline const Eigen::SparseMatrix<double>& stoichiometric_matrix_products() const      {return stoichiometric_matrix_products_;};
  inline const Eigen::SparseMatrix<double>& stoichiometric_matrix_overall() const       {return stoichiometric_matrix_overall_;};
  
  void SetWeights(std::vector<std::vector<double> >& local_dic);
  void SetKeySpecies(const std::vector<std::string>& key_species);
  std::vector<std::vector<double> >& ShortestPaths();

private:
  
  OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;
  OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;
  
  int NS_, NR_;
  
  DirectedGraph DirectedKineticsGraph_;
  
  std::vector<DirectedGraph::vertex_descriptor> species_;
  std::vector<std::vector<DirectedGraph::vertex_descriptor> > species_map_;
  std::vector<DirectedGraph::edge_descriptor> connections_;
  std::vector<Edge>  edge_array_;
  Eigen::SparseMatrix<double> stoichiometric_matrix_reactants_;
  Eigen::SparseMatrix<double> stoichiometric_matrix_products_;
  Eigen::SparseMatrix<double> stoichiometric_matrix_overall_;
  
  boost::property_map<DirectedGraph, boost::edge_weight_t>::type weightmap_;
  boost::property_map<DirectedGraph, boost::vertex_index_t>::type indexmap_;
  
  std::vector<std::vector<double> > target_oic_;
  std::vector<std::string> key_species_;
  
  
  //Functions
  void MemoryAllocation();
  void BuildStoichiometricMatrix();
  void BuildGraph();
  
};


#include "KineticsGraph.hpp"

#endif // KINETICSGRAPH_H
