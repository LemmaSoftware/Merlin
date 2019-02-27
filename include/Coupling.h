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
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2016, University of Utah
 * @copyright Copyright (c) 2016, Lemma Software, LLC
 * @copyright Copyright (c) 2008, Colorado School of Mines
 */
#ifndef  COUPLING_INC
#define  COUPLING_INC

#pragma once
#include "LemmaObject.h"
#include "LayeredEarthEM.h"
#include "PolygonalWireAntenna.h"
#include "EMEarth1D.h"

#ifdef LEMMAUSEVTK
//#include "vtkHyperOctree.h"
//#include "vtkHyperOctreeCursor.h"
//#include "vtkXMLHyperOctreeWriter.h"
#include "vtkHyperTreeGrid.h"
#include "vtkHyperTreeCursor.h"
#include "vtkXMLHyperTreeGridWriter.h"
#include "vtkDoubleArray.h"
#endif

namespace Lemma {

    /**
     * \ingroup Merlin
     * \brief
     * \details
     */
    class Coupling : public LemmaObject {

        friend std::ostream &operator<<(std::ostream &stream, const Coupling &ob);

        public:

        // ====================  LIFECYCLE     =======================

        /**
         * Default constructor.
         * @note This method is locked, and cannot be called directly.
         *       The reason that the method is public is to enable the use
         *       of make_shared whilst enforcing the use of shared_ptr,
         *       in c++-17, this curiosity may be resolved.
         * @see Coupling::NewSP
         */
        explicit Coupling ( const ctor_key& );

        /**
         * DeSerializing constructor.
         * @note This method is locked, and cannot be called directly.
         *       The reason that the method is public is to enable the use
         *       of make_shared whilst enforcing the use of shared_ptr,
         *       in c++-17, this curiosity may be resolved.
         * @see Coupling::DeSerialize
         */
        Coupling ( const YAML::Node& node, const ctor_key& );

        /**
         * Default destructor.
         * @note This method should never be called due to the mandated
         *       use of smart pointers. It is necessary to keep the method
         *       public in order to allow for the use of the more efficient
         *       make_shared constructor.
         */
        virtual ~Coupling ();

        /**
         *  Uses YAML to serialize this object.
         *  @return a YAML::Node
         *  @see Coupling::DeSerialize
         */
        virtual YAML::Node Serialize() const;

        /*
         *  Factory method for generating concrete class.
         *  @return a std::shared_ptr of type Coupling
         */
        static std::shared_ptr< Coupling > NewSP();

        /**
         *   Constructs an Coupling object from a YAML::Node.
         *   @see Coupling::Serialize
         */
        static std::shared_ptr<Coupling> DeSerialize(const YAML::Node& node);

        // ====================  OPERATORS     =======================

        // ====================  OPERATIONS    =======================

        /**
         * @return std::shared_ptr<LayeredEarthEM>
         */
        inline std::shared_ptr<LayeredEarthEM> GetSigmaModel (  ) {
            return SigmaModel;
        }		// -----  end of method Coupling::get_SigmaModel  -----

        /**
         * @param[in] value the 1D-EM model used for calculations
         */
        inline void SetLayeredEarthEM ( std::shared_ptr< LayeredEarthEM > value ) {
            SigmaModel	= value;
            return ;
        }		// -----  end of method Coupling::set_SigmaModel  -----

        /**
         *
         */
        inline void SetIntegrationSize ( const Vector3r& size ) {
            Size = size;
            return ;
        }		// -----  end of method Coupling::SetIntegrationSize  -----

        /**
         *
         */
        inline void SetIntegrationOrigin ( const Vector3r& origin ) {
            Origin = origin;
            return ;
        }		// -----  end of method Coupling::SetIntegrationOrigin  -----

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
        Complex Calculate (const std::vector< std::string >& tx, const std::vector< std::string >& rx,
                bool vtkOutput=false );

        /**
         *  Sets the tolerance to use for making the adaptive mesh
         *
         */
        inline void SetTolerance(const Real& ttol) {
            tol = ttol;
        }

        inline void SetPulseDuration(const Real& taup) {
            Taup = taup;
        }

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
        Coupling( const Coupling& ) = delete;

        private:

        /**
         *  Returns the kernel value for an input prism
         */
        Complex f( const Vector3r& r, const Real& volume , const Vector3cr& Ht, const Vector3cr& Hr);

        void IntegrateOnOctreeGrid( bool vtkOutput=false );

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
        void EvaluateKids(  const Vector3r& size, const int& level, const Vector3r& cpos,
                            const Complex& parentVal );

        #ifdef LEMMAUSEVTK6
        /**
         *  Same functionality as @see EvaluateKids, but includes generation of a VTK
         *  HyperOctree, which is useful for visualization.
         */
        void EvaluateKids2(  const Vector3r& size, const int& level, const Vector3r& cpos,
                            const Complex& parentVal, vtkHyperOctree* octree, vtkHyperOctreeCursor* curse );

        void GetPosition( vtkHyperOctreeCursor* Cursor, Real* p );
        #endif

        // ====================  DATA MEMBERS  =========================

        int                                       ilay;
        int                                       nleaves;
        int                                       minLevel=4;
        int                                       maxLevel=8;

        Real                                      VOLSUM;
        Real                                      tol=1e-11;
        Real                                      Taup = .020;  // Sec

        Complex                                   SUM;

        Vector3r                                  Size;
        Vector3r                                  Origin;

        std::shared_ptr< LayeredEarthEM >         SigmaModel = nullptr;
        std::shared_ptr< FieldPoints >            cpoints;

        std::map< std::string , std::shared_ptr< PolygonalWireAntenna > >  TxRx;
        std::map< std::string , std::shared_ptr< EMEarth1D > >             EMEarths;

        #ifdef LEMMAUSEVTK
        std::map< int, Complex  >                 LeafDict;
        std::map< int, int     >                  LeafDictIdx;
        std::map< int, Real     >                 LeafDictErr;
        #endif

        /** ASCII string representation of the class name */
        static constexpr auto CName = "Coupling";

    }; // -----  end of class  Coupling  -----
}  // -----  end of namespace Lemma ----

/* vim: set tabstop=4 expandtab */
/* vim: set filetype=cpp */

#endif   // ----- #ifndef COUPLING_INC  -----
