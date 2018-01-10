/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      08/28/2017 11:51:33 AM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2017, University of Utah
 * @copyright Copyright (c) 2017, Lemma Software, LLC
 */


#include "MerlinObject.h"

namespace Lemma {

    // ====================  FRIEND METHODS  =====================

    std::ostream &operator << (std::ostream &stream, const MerlinObject &ob) {
        stream << ob.Serialize()  << "\n"; // End of doc ---
        return stream;
    }

    // ====================  LIFECYCLE     =======================

    //--------------------------------------------------------------------------------------
    //       Class:  MerlinObject
    //      Method:  MerlinObject
    // Description:  constructor (locked)
    //--------------------------------------------------------------------------------------
    MerlinObject::MerlinObject ( const ctor_key& key ) : LemmaObject( key ) {

    }  // -----  end of method MerlinObject::MerlinObject  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  MerlinObject
    //      Method:  MerlinObject
    // Description:  DeSerializing constructor (locked)
    //--------------------------------------------------------------------------------------
    MerlinObject::MerlinObject (const YAML::Node& node, const ctor_key& key ) : LemmaObject(node, key) {

    }  // -----  end of method MerlinObject::MerlinObject  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  MerlinObject
    //      Method:  NewSP()
    // Description:  public constructor returing a shared_ptr
    //--------------------------------------------------------------------------------------
    std::shared_ptr< MerlinObject >  MerlinObject::NewSP() {
        return std::make_shared< MerlinObject >( ctor_key() );
    }

    //--------------------------------------------------------------------------------------
    //       Class:  MerlinObject
    //      Method:  ~MerlinObject
    // Description:  destructor (protected)
    //--------------------------------------------------------------------------------------
    MerlinObject::~MerlinObject () {

    }  // -----  end of method MerlinObject::~MerlinObject  (destructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  MerlinObject
    //      Method:  Serialize
    //--------------------------------------------------------------------------------------
    YAML::Node  MerlinObject::Serialize (  ) const {
        YAML::Node node = LemmaObject::Serialize();
        node.SetTag( GetName() );
        // FILL IN CLASS SPECIFICS HERE
        node["Merlin_VERSION"] = MERLIN_VERSION;
        return node;
    }		// -----  end of method MerlinObject::Serialize  -----

    //--------------------------------------------------------------------------------------
    //       Class:  MerlinObject
    //      Method:  DeSerialize
    //--------------------------------------------------------------------------------------
    std::shared_ptr<MerlinObject> MerlinObject::DeSerialize ( const YAML::Node& node  ) {
        if (node.Tag() !=  "MerlinObject" ) {
            throw  DeSerializeTypeMismatch( "MerlinObject", node.Tag());
        }
        return std::make_shared< MerlinObject > ( node, ctor_key() );
    }		// -----  end of method MerlinObject::DeSerialize  -----

} // ----  end of namespace Lemma  ----

/* vim: set tabstop=4 expandtab: */
/* vim: set filetype=cpp: */


