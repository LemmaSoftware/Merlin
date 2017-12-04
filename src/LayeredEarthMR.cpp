/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      08/28/2017 03:32:34 PM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2017, University of Utah
 * @copyright Copyright (c) 2017, Lemma Software, LLC
 */

#include "LayeredEarthMR.h"

namespace Lemma {

    // ====================  FRIEND METHODS  =====================

    std::ostream &operator << (std::ostream &stream, const LayeredEarthMR &ob) {
        stream << ob.Serialize()  << "\n"; // End of doc ---
        return stream;
    }

    // ====================  LIFECYCLE     =======================

    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  LayeredEarthMR
    // Description:  constructor (locked)
    //--------------------------------------------------------------------------------------
    LayeredEarthMR::LayeredEarthMR (const ctor_key&) : LayeredEarth( ) {

    }  // -----  end of method LayeredEarthMR::LayeredEarthMR  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  LayeredEarthMR
    // Description:  DeSerializing constructor (locked)
    //--------------------------------------------------------------------------------------
    LayeredEarthMR::LayeredEarthMR (const YAML::Node& node, const ctor_key&) : LayeredEarth(node) {
        Interfaces = node["Interfaces"].as<VectorXr>();
        T2StarBins = node["T2StarBins"].as<VectorXr>();
        T2StarBinEdges = node["T2StarBinEdges"].as<VectorXr>();
        ModelMat = node["ModelMat"].as<MatrixXr>();
    }  // -----  end of method LayeredEarthMR::LayeredEarthMR  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  NewSP()
    // Description:  public constructor returing a shared_ptr
    //--------------------------------------------------------------------------------------
    std::shared_ptr< LayeredEarthMR >  LayeredEarthMR::NewSP() {
        return std::make_shared< LayeredEarthMR >( ctor_key() );
    }

    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  ~LayeredEarthMR
    // Description:  destructor (protected)
    //--------------------------------------------------------------------------------------
    LayeredEarthMR::~LayeredEarthMR () {

    }  // -----  end of method LayeredEarthMR::~LayeredEarthMR  (destructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  Serialize
    //--------------------------------------------------------------------------------------
    YAML::Node  LayeredEarthMR::Serialize (  ) const {
        YAML::Node node = LayeredEarth::Serialize();
        node.SetTag( GetName() );
        // FILL IN CLASS SPECIFICS HERE
        node["Merlin_VERSION"] = MERLIN_VERSION;
        node["Interfaces"] = Interfaces;
        node["T2StarBins"] = T2StarBins;
        node["T2StarBinEdges"] = T2StarBinEdges;
        node["ModelMat"] = ModelMat;
        return node;
    }		// -----  end of method LayeredEarthMR::Serialize  -----

    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  DeSerialize
    //--------------------------------------------------------------------------------------
    std::shared_ptr<LayeredEarthMR> LayeredEarthMR::DeSerialize ( const YAML::Node& node  ) {
        if (node.Tag() !=  "LayeredEarthMR" ) {
            throw  DeSerializeTypeMismatch( "LayeredEarthMR", node.Tag());
        }
        return std::make_shared< LayeredEarthMR > ( node, ctor_key() );
    }		// -----  end of method LayeredEarthMR::DeSerialize  -----

    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  SetNumberOfLayers
    //--------------------------------------------------------------------------------------
    void LayeredEarthMR::SetNumberOfLayers ( const int& nlay  ) {
        NumberOfLayers = nlay;
        NumberOfInterfaces = nlay+1;
        return ;
    }		// -----  end of method LayeredEarthMR::SetNumberOfLayers  -----

    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  AlignWithKernel
    //--------------------------------------------------------------------------------------
    void LayeredEarthMR::AlignWithKernel ( std::shared_ptr<KernelV0> Kern ) {
        int nlay = Kern->GetInterfaces().size()-1;
        SetNumberOfLayers( nlay );
        Interfaces = Kern->GetInterfaces();
        LayerThickness = Kern->GetInterfaces().tail(nlay) - Kern->GetInterfaces().head(nlay) ;
        SetMagneticFieldComponents( Kern->GetSigmaModel()->GetMagneticField(), TESLA);
        return ;
    }		// -----  end of method LayeredEarthMR::AlignWithKernel  -----

    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  SetNumberOfT2StarBins
    //--------------------------------------------------------------------------------------
    void LayeredEarthMR::SetT2StarBins ( const Real& first, const Real& last, const int& nT2 ) {
        T2StarBinEdges = VectorXr::Zero(nT2+1);
        Real m = 1./(nT2);
        Real quotient = std::pow(last/first, m);
        T2StarBinEdges[0] = first;
        for (int i=1; i<nT2+1; ++i) {
            T2StarBinEdges[i] = T2StarBinEdges[i-1]*quotient;
        }
        T2StarBins = (T2StarBinEdges.head(nT2) + T2StarBinEdges.tail(nT2)) / 2;
        InitModelMat();
        return;
    }		// -----  end of method LayeredEarthMR::SetNumberOfT2StarBins  -----


    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  InitModelMat
    //--------------------------------------------------------------------------------------
    void LayeredEarthMR::InitModelMat (  ) {
        ModelMat = MatrixXr::Zero( NumberOfLayers, T2StarBins.size() );
        return ;
    }		// -----  end of method LayeredEarthMR::InitModelMat  -----


    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  GetT2StarBins
    //--------------------------------------------------------------------------------------
    VectorXr LayeredEarthMR::GetT2StarBins (  ) {
        return T2StarBins ;
    }		// -----  end of method LayeredEarthMR::GetT2StarBins  -----


    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  GetModelMatrix
    //--------------------------------------------------------------------------------------
    MatrixXr LayeredEarthMR::GetModelMatrix (  ) {
        return ModelMat ;
    }		// -----  end of method LayeredEarthMR::GetModelMatrix  -----


    //--------------------------------------------------------------------------------------
    //       Class:  LayeredEarthMR
    //      Method:  GetModelVector
    //--------------------------------------------------------------------------------------
    VectorXr LayeredEarthMR::GetModelVector (  ) {
        //VectorXr B(Eigen::Map<VectorXr>(ModelMat.data(), ModelMat.cols()*ModelMat.rows()));
        MatrixXr MT = ModelMat.transpose();
        VectorXr B(Eigen::Map<VectorXr>(MT.data(), MT.cols()*MT.rows()));
        //VectorXr B(Eigen::Map<VectorXr>(ModelMat.transpose().data(), ModelMat.cols()*ModelMat.rows()));
        return B;
    }		// -----  end of method LayeredEarthMR::GetModelVector  -----


} // ----  end of namespace Lemma  ----

/* vim: set tabstop=4 expandtab: */
/* vim: set filetype=cpp: */


