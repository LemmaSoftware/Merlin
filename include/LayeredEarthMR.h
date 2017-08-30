/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      08/28/2017 03:32:26 PM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2017, University of Utah
 * @copyright Copyright (c) 2017, Lemma Software, LLC
 */

#ifndef  LAYEREDEARTHMR_INC
#define  LAYEREDEARTHMR_INC


#pragma once
#include "LayeredEarth.h"
#include "MerlinConfig.h"
#include "KernelV0.h"

namespace Lemma {

    /**
     * \ingroup Merlin
     * \brief
     * \details
     */
    class LayeredEarthMR : public LayeredEarth {

        friend std::ostream &operator<<(std::ostream &stream, const LayeredEarthMR &ob);

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
         * @see LayeredEarthMR::NewSP
         */
        explicit LayeredEarthMR ( const ctor_key& );

        /**
         * DeSerializing constructor.
         * @note This method is locked, and cannot be called directly.
         *       The reason that the method is public is to enable the use
         *       of make_shared whilst enforcing the use of shared_ptr,
         *       in c++-17, this curiosity may be resolved.
         * @see LayeredEarthMR::DeSerialize
         */
        LayeredEarthMR ( const YAML::Node& node, const ctor_key& );

        /**
         * Default destructor.
         * @note This method should never be called due to the mandated
         *       use of smart pointers. It is necessary to keep the method
         *       public in order to allow for the use of the more efficient
         *       make_shared constructor.
         */
        virtual ~LayeredEarthMR ();

        /**
         *  Uses YAML to serialize this object.
         *  @return a YAML::Node
         *  @see LayeredEarthMR::DeSerialize
         */
        virtual YAML::Node Serialize() const;

        /**
         *  Factory method for generating concrete class.
         *  @return a std::shared_ptr of type LayeredEarthMR
         */
        static std::shared_ptr< LayeredEarthMR > NewSP();

        /**
         *   Constructs an LayeredEarthMR object from a YAML::Node.
         *   @see LayeredEarthMR::Serialize
         */
        static std::shared_ptr<LayeredEarthMR> DeSerialize(const YAML::Node& node);

        // ====================  OPERATORS     =======================

        // ====================  OPERATIONS    =======================

        // ====================  ACCESS        =======================



        /**
         *  Sets the T2StarBins to solve for, these are log spaced
         *  @param[in] first is the beginning of the bins
         *  @param[in] last is the end of the bins
         *  @param[in] nT2 is the number of bins
         */
        void SetT2StarBins(const Real& first, const Real& last, const int& nT2);

        /**
         *  Convenience method, that aligns model with a Kernel
         *  @param[in] Kern input kernel to align with
         */
        void AlignWithKernel( std::shared_ptr<KernelV0> Kern );

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
        LayeredEarthMR( const LayeredEarthMR& ) = delete;

        // ====================  DATA MEMBERS  =========================

        private:

        /** ASCII string representation of the class name */
        static constexpr auto CName = "LayeredEarthMR";

        /**
         * Sets the number of layers
         */
        void SetNumberOfLayers(const int& nlay);

        /** Initializes the model matrix */
        void InitModelMat();

        VectorXr Interfaces;      // Layer interfaces, for pcolor
        VectorXr T2StarBins;      // the actual T2* values
        VectorXr T2StarBinEdges;  // Convenience, for pcolor
        MatrixXr ModelMat;        // The NMR model, in matrix form

    }; // -----  end of class  LayeredEarthMR  -----
}  // -----  end of namespace Lemma ----

/* vim: set tabstop=4 expandtab: */
/* vim: set filetype=cpp: */


#endif   // ----- #ifndef LAYEREDEARTHMR_INC  -----

