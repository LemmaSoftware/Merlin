/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      08/28/2017 08:49:51 AM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2017, University of Utah
 * @copyright Copyright (c) 2017, Lemma Software, LLC
 */

#ifndef  FORWARDFID_INC
#define  FORWARDFID_INC

#pragma once
#include <random>
#include "LayeredEarthEM.h"
#include "PolygonalWireAntenna.h"
#include "EMEarth1D.h"

#include "MerlinObject.h"
#include "DataFID.h"
#include "KernelV0.h"
#include "LayeredEarthMR.h"

namespace Lemma {

    /**
     * \ingroup Merlin
     * \brief Forward modelling for FID sNMR pulses
     * \details This class performs forward modelling of sNMR
     *          FID experiments.
     */
    class ForwardFID : public MerlinObject {

        friend std::ostream &operator<<(std::ostream &stream, const ForwardFID &ob);

        public:

        // ====================  LIFECYCLE     =======================

        /**
         * Default constructor.
         * @note This method is locked, and cannot be called directly.
         *       The reason that the method is public is to enable the use
         *       of make_shared whilst enforcing the use of shared_ptr,
         *       in c++-17, this curiosity may be resolved.
         * @see ForwardFID::NewSP
         */
        explicit ForwardFID ( const ctor_key& );

        /**
         * DeSerializing constructor.
         * @note This method is locked, and cannot be called directly.
         *       The reason that the method is public is to enable the use
         *       of make_shared whilst enforcing the use of shared_ptr,
         *       in c++-17, this curiosity may be resolved.
         * @see ForwardFID::DeSerialize
         */
        ForwardFID ( const YAML::Node& node, const ctor_key& );

        /**
         * Default destructor.
         * @note This method should never be called due to the mandated
         *       use of smart pointers. It is necessary to keep the method
         *       public in order to allow for the use of the more efficient
         *       make_shared constructor.
         */
        virtual ~ForwardFID ();

        /**
         *  Uses YAML to serialize this object.
         *  @return a YAML::Node
         *  @see ForwardFID::DeSerialize
         */
        virtual YAML::Node Serialize() const;

        /*
         *  Factory method for generating concrete class.
         *  @return a std::shared_ptr of type ForwardFID
         */
        static std::shared_ptr< ForwardFID > NewSP();

        /**
         *   Constructs an ForwardFID object from a YAML::Node.
         *   @see ForwardFID::Serialize
         */
        static std::shared_ptr<ForwardFID> DeSerialize(const YAML::Node& node);

        // ====================  OPERATORS     =======================

        // ====================  OPERATIONS    =======================

        /**
         *  Performs forward model calculation based on input parameters
         *  @return Merlin class representing data
         */
        std::shared_ptr<DataFID> ForwardModel( std::shared_ptr<LayeredEarthMR> );

        // ====================  ACCESS        =======================
        /**
         *  Sets windows using the edges in a Eigen VectorXr
         *  @param[in] Edges are the edges of the time gate windows
         */
        void SetWindowEdges( const VectorXr& Edges );

        /**
         *  Sets the windows for calculation as a series of log spaced windows.
         *  This method calculates the window edges, so the centres will be n-1
         *  @param[in] first is the beginning of the vector
         *  @param[in] last is the last element in the vector
         *  @param[in] n is the number of elements
         */
        void SetLogSpacedWindows(const Real& first, const Real& last, const int& n);

        /**
         *   @param[in] K0 is the initial amplitude imaging kernel
         */
        void SetKernel( std::shared_ptr< KernelV0 > K0 );

        /**
         *  @param[in] floor is the standard deviation of the noise (zero mean),
         *             defaults to zero
         */
        void SetNoiseFloor( const Real& floor );

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
        ForwardFID( const ForwardFID& ) = delete;


        private:

        /**
         * Calculates the QT matrix
         * @param[in] T2StarBins are the T2* bins to use
         */
        void CalcQTMatrix( VectorXr T2StarBins );

        // ====================  DATA MEMBERS  =========================
        /** ASCII string representation of the class name */
        static constexpr auto CName = "ForwardFID";

        /** Imaging kernel used in calculation */
        std::shared_ptr< KernelV0 >   Kernel = nullptr;

        /** Noise floor for additive Gaussian noise */
        Real NoiseFloor = 0.;

        /** Time gate windows */
        VectorXr WindowEdges;

        /** Time gate centres */
        VectorXr WindowCentres;

        /** QT matrix */
        MatrixXcr QT;

        /** Include RDP effects? */
        bool RDP = false;

    }; // -----  end of class  ForwardFID  -----

}  // -----  end of namespace Lemma ----

/* vim: set tabstop=4 expandtab: */
/* vim: set filetype=cpp: */

#endif   // ----- #ifndef FORWARDFID_INC  -----
