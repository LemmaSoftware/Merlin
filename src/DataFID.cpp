/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      08/28/2017 01:04:19 PM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2017, University of Utah
 * @copyright Copyright (c) 2017, Lemma Software, LLC
 */


#include "DataFID.h"

namespace Lemma {

    // ====================  FRIEND METHODS  =====================

    std::ostream &operator << (std::ostream &stream, const DataFID &ob) {
        stream << ob.Serialize()  << "\n---\n"; // End of doc ---
        return stream;
    }

    // ====================  LIFECYCLE     =======================

    //--------------------------------------------------------------------------------------
    //       Class:  DataFID
    //      Method:  DataFID
    // Description:  constructor (locked)
    //--------------------------------------------------------------------------------------
    DataFID::DataFID (const ctor_key&) : MerlinObject( ) {

    }  // -----  end of method DataFID::DataFID  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  DataFID
    //      Method:  DataFID
    // Description:  DeSerializing constructor (locked)
    //--------------------------------------------------------------------------------------
    DataFID::DataFID (const YAML::Node& node, const ctor_key&) : MerlinObject(node) {

    }  // -----  end of method DataFID::DataFID  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  DataFID
    //      Method:  NewSP()
    // Description:  public constructor returing a shared_ptr
    //--------------------------------------------------------------------------------------
    std::shared_ptr< DataFID >  DataFID::NewSP() {
        return std::make_shared< DataFID >( ctor_key() );
    }

    //--------------------------------------------------------------------------------------
    //       Class:  DataFID
    //      Method:  ~DataFID
    // Description:  destructor (protected)
    //--------------------------------------------------------------------------------------
    DataFID::~DataFID () {

    }  // -----  end of method DataFID::~DataFID  (destructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  DataFID
    //      Method:  Serialize
    //--------------------------------------------------------------------------------------
    YAML::Node  DataFID::Serialize (  ) const {
        YAML::Node node = MerlinObject::Serialize();
        node.SetTag( GetName() );
        // FILL IN CLASS SPECIFICS HERE
        return node;
    }		// -----  end of method DataFID::Serialize  -----

    //--------------------------------------------------------------------------------------
    //       Class:  DataFID
    //      Method:  DeSerialize
    //--------------------------------------------------------------------------------------
    std::shared_ptr<DataFID> DataFID::DeSerialize ( const YAML::Node& node  ) {
        if (node.Tag() !=  "DataFID" ) {
            throw  DeSerializeTypeMismatch( "DataFID", node.Tag());
        }
        return std::make_shared< DataFID > ( node, ctor_key() );
    }		// -----  end of method DataFID::DeSerialize  -----

} // ----  end of namespace Lemma  ----

/* vim: set tabstop=4 expandtab: */
/* vim: set filetype=cpp: */


