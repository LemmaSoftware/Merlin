/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      08/28/2017 09:03:03 AM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2017, University of Utah
 * @copyright Copyright (c) 2017, Lemma Software, LLC
 */

#include "ForwardFID.h"

namespace Lemma {

    // ====================  FRIEND METHODS  =====================

    std::ostream &operator << (std::ostream &stream, const ForwardFID &ob) {
        stream << ob.Serialize()  << "\n---\n"; // End of doc ---
        return stream;
    }

    // ====================  LIFECYCLE     =======================

    //--------------------------------------------------------------------------------------
    //       Class:  ForwardFID
    //      Method:  ForwardFID
    // Description:  constructor (locked)
    //--------------------------------------------------------------------------------------
    ForwardFID::ForwardFID (const ctor_key&) : MerlinObject( ) {

    }  // -----  end of method ForwardFID::ForwardFID  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  ForwardFID
    //      Method:  ForwardFID
    // Description:  DeSerializing constructor (locked)
    //--------------------------------------------------------------------------------------
    ForwardFID::ForwardFID (const YAML::Node& node, const ctor_key&) : MerlinObject(node) {

    }  // -----  end of method ForwardFID::ForwardFID  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  ForwardFID
    //      Method:  NewSP()
    // Description:  public constructor returing a shared_ptr
    //--------------------------------------------------------------------------------------
    std::shared_ptr< ForwardFID >  ForwardFID::NewSP() {
        return std::make_shared< ForwardFID >( ctor_key() );
    }

    //--------------------------------------------------------------------------------------
    //       Class:  ForwardFID
    //      Method:  ~ForwardFID
    // Description:  destructor (protected)
    //--------------------------------------------------------------------------------------
    ForwardFID::~ForwardFID () {

    }  // -----  end of method ForwardFID::~ForwardFID  (destructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  ForwardFID
    //      Method:  Serialize
    //--------------------------------------------------------------------------------------
    YAML::Node  ForwardFID::Serialize (  ) const {
        YAML::Node node = MerlinObject::Serialize();
        node.SetTag( GetName() );
        // FILL IN CLASS SPECIFICS HERE
        if (Kernel != nullptr) {
            node["Kernel"] = Kernel->Serialize();
        }
        node["WindowEdges"] = WindowEdges;
        node["WindowCentres"] = WindowCentres;
        node["RDP"] = RDP;
        return node;
    }		// -----  end of method ForwardFID::Serialize  -----

    //--------------------------------------------------------------------------------------
    //       Class:  ForwardFID
    //      Method:  DeSerialize
    //--------------------------------------------------------------------------------------
    std::shared_ptr<ForwardFID> ForwardFID::DeSerialize ( const YAML::Node& node  ) {
        if (node.Tag() !=  "ForwardFID" ) {
            throw  DeSerializeTypeMismatch( "ForwardFID", node.Tag());
        }
        return std::make_shared< ForwardFID > ( node, ctor_key() );
    }		// -----  end of method ForwardFID::DeSerialize  -----


    //--------------------------------------------------------------------------------------
    //       Class:  ForwardFID
    //      Method:  SetWindowEdges
    //--------------------------------------------------------------------------------------
    void ForwardFID::SetWindowEdges ( const VectorXr& Edges ) {
        WindowEdges = Edges;
        return ;
    }		// -----  end of method ForwardFID::SetWindowEdges  -----


    //--------------------------------------------------------------------------------------
    //       Class:  ForwardFID
    //      Method:  ForwardModel
    //--------------------------------------------------------------------------------------
    std::shared_ptr<DataFID> ForwardFID::ForwardModel (  ) {
        auto FID = DataFID::NewSP();
        return FID;
    }		// -----  end of method ForwardFID::ForwardModel  -----

    //--------------------------------------------------------------------------------------
    //       Class:  ForwardFID
    //      Method:  LogSpaced
    //--------------------------------------------------------------------------------------
    void ForwardFID::SetLogSpacedWindows ( const Real& first, const Real& last,
        const int& n ) {
        WindowEdges = VectorXr::Zero(n);
        Real m = 1./(n-1.);
        Real quotient = std::pow(last/first, m);
        WindowEdges[0] = first;
        for (int i=1; i<n; ++i) {
            WindowEdges[i] = WindowEdges[i-1]*quotient;
        }
        WindowCentres = (WindowEdges.head(n-1) + WindowEdges.tail(n-1)) / 2;
        return;
    }		// -----  end of method ForwardFID::LogSpaced  -----

    //--------------------------------------------------------------------------------------
    //       Class:  ForwardFID
    //      Method:  CalcQTMatrix
    //--------------------------------------------------------------------------------------
    void ForwardFID::CalcQTMatrix (  ) {
        return ;
    }		// -----  end of method ForwardFID::CalcQTMatrix  -----

} // ----  end of namespace Lemma  ----

/* vim: set tabstop=4 expandtab: */
/* vim: set filetype=cpp: */


