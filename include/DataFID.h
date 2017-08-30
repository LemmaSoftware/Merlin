/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      08/28/2017 01:04:12 PM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2017, University of Utah
 * @copyright Copyright (c) 2017, Lemma Software, LLC
 */

#ifndef  DATAFID_INC
#define  DATAFID_INC

#pragma once
#include "MerlinObject.h"

namespace Lemma {

    /**
     * \ingroup Merlin
     * \brief   Base class for representation of FID data in Merlin.
     * \details Simplest form of FID data used in Merlin. Field data derives
                from this class for specific instrumentation details. See
                the Akvo project for further details (akvo.lemmasoftware.org).
     */
    class DataFID : public MerlinObject {

        friend std::ostream &operator<<(std::ostream &stream, const DataFID &ob);

        friend class ForwardFID;

        protected:
        /*
         *  This key is used to lock the constructor. It is protected so that inhereted
         *  classes also have the key to contruct their base class.
         */
        struct ctor_key {};

        public:

        // ====================  LIFECYCLE     =======================

        /**
         * Default constructor.
         * @note This method is locked, and cannot be called directly.
         *       The reason that the method is public is to enable the use
         *       of make_shared whilst enforcing the use of shared_ptr,
         *       in c++-17, this curiosity may be resolved.
         * @see DataFID::NewSP
         */
        explicit DataFID ( const ctor_key& );

        /**
         * DeSerializing constructor.
         * @note This method is locked, and cannot be called directly.
         *       The reason that the method is public is to enable the use
         *       of make_shared whilst enforcing the use of shared_ptr,
         *       in c++-17, this curiosity may be resolved.
         * @see DataFID::DeSerialize
         */
        DataFID ( const YAML::Node& node, const ctor_key& );

        /**
         * Default destructor.
         * @note This method should never be called due to the mandated
         *       use of smart pointers. It is necessary to keep the method
         *       public in order to allow for the use of the more efficient
         *       make_shared constructor.
         */
        virtual ~DataFID ();

        /**
         *  Uses YAML to serialize this object.
         *  @return a YAML::Node
         *  @see DataFID::DeSerialize
         */
        virtual YAML::Node Serialize() const;

        /*
         *  Factory method for generating concrete class.
         *  @return a std::shared_ptr of type DataFID
         */
        static std::shared_ptr< DataFID > NewSP();

        /**
         *   Constructs an DataFID object from a YAML::Node.
         *   @see DataFID::Serialize
         */
        static std::shared_ptr<DataFID> DeSerialize(const YAML::Node& node);

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
        DataFID( const DataFID& ) = delete;

        // ====================  DATA MEMBERS  =========================

        private:

        /** ASCII string representation of the class name */
        static constexpr auto CName = "DataFID";

        /** Holds the actual data */
        MatrixXcr   FIDData;

        VectorXr    WindowCentres;

        VectorXr    WindowEdges;

        VectorXr    PulseMoment;

    }; // -----  end of class  DataFID  -----
}  // -----  end of namespace Lemma ----

/* vim: set tabstop=4 expandtab: */
/* vim: set filetype=cpp: */


#endif   // ----- #ifndef DATAFID_INC  -----
