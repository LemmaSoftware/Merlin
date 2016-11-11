/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      11/11/2016 01:47:34 PM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2016, University of Utah
 * @copyright Copyright (c) 2016, Lemma Software, LLC
 */
#ifndef  KERNELV0_INC
#define  KERNELV0_INC

#pragma once
#include "LemmaObject.h"
#include "LayeredEarthEM.h"
#include "PolygonalWireAntenna.h"

namespace Lemma {

    /**
     * \ingroup Merlin
     * \brief
     * \details
     */
    class KernelV0 : public LemmaObject {

        friend std::ostream &operator<<(std::ostream &stream, const KernelV0 &ob);

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
         * @see KernelV0::NewSP
         */
        explicit KernelV0 ( const ctor_key& );

        /**
         * DeSerializing constructor.
         * @note This method is locked, and cannot be called directly.
         *       The reason that the method is public is to enable the use
         *       of make_shared whilst enforcing the use of shared_ptr,
         *       in c++-17, this curiosity may be resolved.
         * @see KernelV0::DeSerialize
         */
        KernelV0 ( const YAML::Node& node, const ctor_key& );

        /**
         * Default destructor.
         * @note This method should never be called due to the mandated
         *       use of smart pointers. It is necessary to keep the method
         *       public in order to allow for the use of the more efficient
         *       make_shared constructor.
         */
        virtual ~KernelV0 ();

        /**
         *  Uses YAML to serialize this object.
         *  @return a YAML::Node
         *  @see KernelV0::DeSerialize
         */
        virtual YAML::Node Serialize() const;

        /*
         *  Factory method for generating concrete class.
         *  @return a std::shared_ptr of type KernelV0
         */
        static std::shared_ptr< KernelV0 > NewSP();

        /**
         *   Constructs an KernelV0 object from a YAML::Node.
         *   @see KernelV0::Serialize
         */
        static std::shared_ptr<KernelV0> DeSerialize(const YAML::Node& node);

        // ====================  OPERATORS     =======================

        // ====================  OPERATIONS    =======================

        /**
         * @return std::shared_ptr<LayeredEarthEM>
         */
        inline std::shared_ptr<LayeredEarthEM> GetSigmaModel (  ) {
            return SigmaModel;
        }		// -----  end of method KernelV0::get_SigmaModel  -----

        /**
         * @param[in] value the 1D-EM model used for calculations
         */
        inline void SetLayeredEarthEM ( std::shared_ptr< LayeredEarthEM > value ) {
            SigmaModel	= value;
            return ;
        }		// -----  end of method KernelV0::set_SigmaModel  -----

        /**
         *   Assign transmiter coils
         */
        inline void PushCoil( const std::string& label, std::shared_ptr<PolygonalWireAntenna> ant ) {
            TxRx[label] = ant;
        }

        /**
         *
         */
        void CalculateK0 (const std::vector< std::string >& tx, const std::vector< std::string >& rx );

        void CalculateK0 (const char* tx, const char* rx );

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
        KernelV0( const KernelV0& ) = delete;

        // ====================  DATA MEMBERS  =========================

        private:

        std::shared_ptr< LayeredEarthEM >        SigmaModel;

        std::map< std::string , std::shared_ptr< PolygonalWireAntenna > >  TxRx;

        /** ASCII string representation of the class name */
        static constexpr auto CName = "KernelV0";

    }; // -----  end of class  KernelV0  -----
}  // -----  end of namespace Lemma ----

/* vim: set tabstop=4 expandtab: */
/* vim: set filetype=cpp: */

#endif   // ----- #ifndef KERNELV0_INC  -----
