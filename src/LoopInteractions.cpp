/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      11/11/2016 01:47:25 PM
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2016, University of Utah
 * @copyright Copyright (c) 2016, Lemma Software, LLC
 * @copyright Copyright (c) 2008, Colorado School of Mines
 */
#include "LoopInteractions.h"

namespace Lemma {

    std::string enum2String(const INTERACTION& type) {
        std::string t;
        switch (type) {
        case COUPLING:
            t = std::string("COUPLING");
            break;
        case INTERFERENCE:
            t = std::string("INTERFERENCE");
            break;
        case PHASE:
            t = std::string("PHASE");
            break;
        }
        return t;
    }

    //--------------------------------------------------------------------------------------
    //       Class:  LoopInteractions
    //      Method:  f
    //--------------------------------------------------------------------------------------
    template <>
    Complex LoopInteractions<COUPLING>::f( const Vector3r& r, const Real& volume, const Vector3cr& Ht, const Vector3cr& Hr ) {
        return volume * ( Ht.dot(Hr) );                              // coupling
    }

    template <>
    Complex LoopInteractions<INTERFERENCE>::f( const Vector3r& r, const Real& volume, const Vector3cr& Ht, const Vector3cr& Hr ) {
        return volume * (1.-((Ht+Hr).norm()/(Hr.norm() + Ht.norm()))); // interference
    }

    template <>
    Complex LoopInteractions<PHASE>::f( const Vector3r& r, const Real& volume, const Vector3cr& Ht, const Vector3cr& Hr ) {
        return volume * std::acos( (Ht.dot(Hr) / (Ht.norm()*Hr.norm())) ); // angle
    }

} // ----  end of namespace Lemma  ----

/* vim: set tabstop=4 expandtab */
/* vim: set filetype=cpp */

