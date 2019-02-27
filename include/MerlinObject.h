/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      08/28/2017 11:40:40 AM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2017, University of Utah
 * @copyright Copyright (c) 2017, Lemma Software, LLC
 */

#ifndef  MERLINOBJECT_INC
#define  MERLINOBJECT_INC


#pragma once
#include "LemmaObject.h"
#include "LemmaObject.h"
#include "MerlinConfig.h"

namespace Lemma {

    /**
     * \ingroup Merlin
     * \brief  Abstract base class for all of Merlin
     * \details
     */
    class MerlinObject : public LemmaObject {

        friend std::ostream &operator<<(std::ostream &stream, const MerlinObject &ob);

        public:

        // ====================  LIFECYCLE     =======================

        /**
         * Default constructor.
         * @note This method is locked, and cannot be called directly.
         *       The reason that the method is public is to enable the use
         *       of make_shared whilst enforcing the use of shared_ptr,
         *       in c++-17, this curiosity may be resolved.
         * @see MerlinObject::NewSP
         */
        explicit MerlinObject ( const ctor_key& key );

        /**
         * DeSerializing constructor.
         * @note This method is locked, and cannot be called directly.
         *       The reason that the method is public is to enable the use
         *       of make_shared whilst enforcing the use of shared_ptr,
         *       in c++-17, this curiosity may be resolved.
         * @see MerlinObject::DeSerialize
         */
        MerlinObject ( const YAML::Node& node, const ctor_key& key );

        /**
         * Default destructor.
         * @note This method should never be called due to the mandated
         *       use of smart pointers. It is necessary to keep the method
         *       public in order to allow for the use of the more efficient
         *       make_shared constructor.
         */
        virtual ~MerlinObject ();

        /**
         *  Uses YAML to serialize this object.
         *  @return a YAML::Node
         *  @see MerlinObject::DeSerialize
         */
        virtual YAML::Node Serialize() const ;

        /*
         *  Factory method for generating concrete class.
         *  @return a std::shared_ptr of type MerlinObject
         */
        static std::shared_ptr< MerlinObject > NewSP();

        /**
         *   Constructs an MerlinObject object from a YAML::Node.
         *   @see MerlinObject::Serialize
         */
        static std::shared_ptr<MerlinObject> DeSerialize(const YAML::Node& node);

        // ====================  OPERATORS     =======================

        // ====================  OPERATIONS    =======================

        // ====================  ACCESS        =======================

        // ====================  INQUIRY       =======================
        /**
         *  Returns the name of the underlying class, similiar to Python's type
         *  @return string of class name
         */
        virtual inline std::string GetName() const {
            return CName;
        }

        protected:

        // ====================  LIFECYCLE     =======================

        /** Copy is disabled */
        MerlinObject( const MerlinObject& ) = delete;

        // ====================  DATA MEMBERS  =========================

        private:

        /** ASCII string representation of the class name */
        static constexpr auto CName = "MerlinObject";

    }; // -----  end of class  MerlinObject  -----

}  // -----  end of namespace Lemma ----

/* vim: set tabstop=4 expandtab: */
/* vim: set filetype=cpp: */


#endif   // ----- #ifndef MERLINOBJECT_INC  -----

