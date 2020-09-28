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

#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>

void KineticsGraph::Setup
(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
 OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML)
{
	thermodynamicsMapXML_ = thermodynamicsMapXML;
	kineticsMapXML_ = kineticsMapXML;
	MemoryAllocation();
	BuildGraph ();
}
/*
KineticsGraph::KineticsGraph (const KineticsGraph& rhs) :
thermodynamicsMapXML_ (rhs.thermodynamicsMapXML_),
kineticsMapXML_ (rhs.kineticsMapXML_)
{
  this->DirectedKineticsGraph_ = rhs.DirectedKineticsGraph_;
  this->NR_ = rhs.NR_;
  this->NS_ = rhs.NS_;
  this->species_ = rhs.species_;
  this->stoichiometric_matrix_products_ = rhs.stoichiometric_matrix_products_;
  this->stoichiometric_matrix_reactants_ = rhs.stoichiometric_matrix_reactants_;
  this->stoichiometric_matrix_overall_ = rhs.stoichiometric_matrix_overall_;
  this->connections_ = rhs.connections_;
  this->edge_array_ = rhs.edge_array_;
  this->weightmap_ = rhs.weightmap_;
  this->indexmap_ = rhs.indexmap_;
}
*/
KineticsGraph::~KineticsGraph ()
{
  species_.clear ();
  species_map_.clear ();
  connections_.clear ();
  edge_array_.clear ();
  target_oic_.clear ();
}

void KineticsGraph::MemoryAllocation ()
{
  species_.resize (thermodynamicsMapXML_->NumberOfSpecies ());
  NS_ = thermodynamicsMapXML_->NumberOfSpecies ();
  NR_ = kineticsMapXML_->NumberOfReactions ();
  edge_array_.reserve (kineticsMapXML_->NumberOfReactions ());
}

void KineticsGraph::BuildGraph ()
{

  //Adding vertices
  for (int i = 0; i < NS_; i++)
    species_[i] = add_vertex (DirectedKineticsGraph_);

  indexmap_ = get (boost::vertex_index, DirectedKineticsGraph_);

  BuildStoichiometricMatrix ();

  //Adding edges
  for (int i = 0; i < stoichiometric_matrix_overall_.outerSize (); i++)
    {
      //cout << i << endl;
      for (Eigen::SparseMatrix<double>::InnerIterator it_reactant (stoichiometric_matrix_overall_, i); it_reactant; ++it_reactant)
        for (Eigen::SparseMatrix<double>::InnerIterator it_product (stoichiometric_matrix_overall_, i); it_product; ++it_product)
          {
            Edge new_connection;
            new_connection.first = it_reactant.index ();
            new_connection.second = it_product.index ();

            edge_array_.push_back (new_connection);
          }
    }

  connections_.resize (edge_array_.size ());

  for (int i = 0; i < connections_.size (); ++i)
    {
      DirectedGraph::edge_descriptor e;
      bool inserted;
      tie (e, inserted) = add_edge (edge_array_[i].first, edge_array_[i].second, DirectedKineticsGraph_);
    }
}

void KineticsGraph::BuildStoichiometricMatrix ()
{
  OpenSMOKE::KineticsMap_CHEMKIN kineticsMapStoichiometry = *kineticsMapXML_;

  stoichiometric_matrix_products_.resize (NS_, NR_);

  stoichiometric_matrix_reactants_.resize (NS_, NR_);

  stoichiometric_matrix_overall_.resize (NS_, NR_);

  for (int z = 0; z < kineticsMapStoichiometry.stoichiometry ().stoichiometric_matrix_products ().outerSize (); z++)
    for (Eigen::SparseMatrix<double>::InnerIterator it (kineticsMapStoichiometry.stoichiometry ().stoichiometric_matrix_products (), z); it; ++it)
      {
        stoichiometric_matrix_products_.insert (z, it.index ()) = it.value ();
        stoichiometric_matrix_overall_.coeffRef (z, it.index ()) += it.value ();
      }

  for (int z = 0; z < kineticsMapStoichiometry.stoichiometry ().stoichiometric_matrix_reactants ().outerSize (); z++)
    for (Eigen::SparseMatrix<double>::InnerIterator it (kineticsMapStoichiometry.stoichiometry ().stoichiometric_matrix_reactants (), z); it; ++it)
      {
        stoichiometric_matrix_reactants_.insert (z, it.index ()) = it.value ();
        stoichiometric_matrix_overall_.coeffRef (z, it.index ()) += it.value ();
      }

  stoichiometric_matrix_reactants_.makeCompressed ();
  stoichiometric_matrix_products_.makeCompressed ();
  stoichiometric_matrix_overall_.makeCompressed ();

}

void KineticsGraph::SetKeySpecies (const std::vector<std::string>& key_species)
{
  key_species_ = key_species;

  target_oic_.resize (key_species_.size ());

  for (int i = 0; i < target_oic_.size (); i++)
    target_oic_[i].resize (thermodynamicsMapXML_->NumberOfSpecies ());

  species_map_.resize (key_species_.size ());
  for (int i = 0; i < species_map_.size (); i++)
    species_map_[i] = species_;
}

void KineticsGraph::SetWeights (std::vector<std::vector<double> >& local_dic)
{
  //Normalizing weighths, if some value > 1 + eps
  double epsilon = 1.e-6;


  for (int i = 0; i < connections_.size (); ++i)
    {
      std::pair < DirectedGraph::edge_descriptor, bool> edgePair = boost::edge (edge_array_[i].first, edge_array_[i].second, DirectedKineticsGraph_);
      DirectedGraph::edge_descriptor edge = edgePair.first;
      if (local_dic[edge_array_[i].first][edge_array_[i].second] <= 1)
        weightmap_[edge] = local_dic[edge_array_[i].first][edge_array_[i].second];
      else if (local_dic[edge_array_[i].first][edge_array_[i].second] > 1 &&
               local_dic[edge_array_[i].first][edge_array_[i].second] <= 1 + epsilon)
        weightmap_[edge] = 1.;
      else if (local_dic[edge_array_[i].first][edge_array_[i].second] > 1 + epsilon)
        OpenSMOKE::FatalErrorMessage ("Direct interaction coefficients higher than 1!");
    }
}

std::vector<std::vector<double> >& KineticsGraph::ShortestPaths ()
{

  //Calculating the shortest paths, according to Dijkstra algorithm
  //(Niemeyer K. and Sung C.J., Comb Flame 158 (8) 1439-1443 (2011))
  for (int i = 0; i < key_species_.size (); i++)
    {
	  std::vector<double> target_distance;
      target_distance.resize (thermodynamicsMapXML_->NumberOfSpecies ());

      boost::dijkstra_shortest_paths_no_color_map (DirectedKineticsGraph_,
                                                   species_[thermodynamicsMapXML_->IndexOfSpecies (key_species_[i]) - 1],
                                                   &species_map_[i][0],
                                                   &target_distance[0],
                                                   weightmap_,
                                                   indexmap_,
                                                   std::greater<double>(),
                                                   std::multiplies<double>(),
                                                   0.,
                                                   1.,
                                                   boost::default_dijkstra_visitor ());

      target_oic_[i] = target_distance;

      target_distance.clear ();
    }

  return target_oic_;
}