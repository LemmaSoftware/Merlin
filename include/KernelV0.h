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
#include "EMEarth1D.h"

#ifdef LEMMAUSEVTK
#include "vtkHyperOctree.h"
#include "vtkHyperOctreeCursor.h"
#include "vtkXMLHyperOctreeWriter.h"
#include "vtkDoubleArray.h"
#endif

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
         *   Calculates a single imaging kernel, however, phased arrays are supported
         *   so that more than one transmitter and/or receiver can be specified.
         *   @param[in] tx is the list of transmitters to use for a kernel, use the same labels as
         *              used in PushCoil.
         *   @param[in] rx is the list of receivers to use for a kernel, use the same labels as
         *              used in PushCoil. @see PushCoil
         *   @param[in] vtkOutput generates a VTK hyperoctree file as well, useful for visualization.
         *              requires compilation of Lemma with VTK.
         */
        void CalculateK0 (const std::vector< std::string >& tx, const std::vector< std::string >& rx,
                bool vtkOutput=false );

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

        private:

        /**
         *  Returns the kernel value for an input prism
         */
        Complex f( const Vector3r& r, const Real& volume , const Vector3cr& Bt);

        void IntegrateOnOctreeGrid( const Real& tolerance , bool vtkOutput=false );

        /**
         *  Recursive call to integrate a function on an adaptive Octree Grid.
         *  For efficiency's sake the octree grid is not stored, as only the
         *  integral (sum) is of interest. The logic for grid refinement is based
         *  on an Octree representation of the domain. If an Octree representation
         *  of the kernel is desired, call alternative version @see EvaluateKids2
         *  @param[in] size gives the domain size, in  metres
         *  @param[in] level gives the current level of the octree grid, call with 0 initially
         *  @param[in] cpos is the centre position of the parent cuboid
         */
        bool EvaluateKids(  const Vector3r& size, const int& level, const Vector3r& cpos,
                            const Complex& parentVal );

        #ifdef LEMMAUSEVTK
        /**
         *  Same functionality as @see EvaluateKids, but includes generation of a VTK
         *  HyperOctree, which is useful for visualization.
         */
        bool EvaluateKids2(  const Vector3r& size, const int& level, const Vector3r& cpos,
                            const Complex& parentVal, vtkHyperOctree* octree, vtkHyperOctreeCursor* curse );

        void GetPosition( vtkHyperOctreeCursor* Cursor, Real* p );
        #endif

        // ====================  DATA MEMBERS  =========================

        Complex                                   SUM;

        Real                                      VOLSUM;

        Real                                      tol=1e-3;

        int                                       nleaves;

        Vector3r   Size;
        Vector3r   Origin;

        std::shared_ptr< LayeredEarthEM >         SigmaModel = nullptr;

        std::map< std::string , std::shared_ptr< PolygonalWireAntenna > >  TxRx;

        std::vector< std::shared_ptr<EMEarth1D> > EMEarths;

        #ifdef LEMMAUSEVTK
        std::map< int, Complex  >                 LeafDict;
        #endif

        /** ASCII string representation of the class name */
        static constexpr auto CName = "KernelV0";

    }; // -----  end of class  KernelV0  -----
}  // -----  end of namespace Lemma ----

/* vim: set tabstop=4 expandtab */
/* vim: set filetype=cpp */

#endif   // ----- #ifndef KERNELV0_INC  -----
