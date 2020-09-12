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
#ifndef  KERNELV0_INC
#define  KERNELV0_INC

#pragma once
#include "MerlinObject.h"
#include "LayeredEarthEM.h"
#include "PolygonalWireAntenna.h"
#include "EMEarth1D.h"

#include "ProgressBar.h"

#ifdef LEMMAUSEVTK
//#include "vtkHyperOctree.h"
//#include "vtkHyperOctreeCursor.h"
//#include "vtkXMLHyperOctreeWriter.h"
//#include "vtkXMLHyperOctreeWriter.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkHyperTree.h"
#include "vtkHyperTree.h"
#include "vtkHyperTreeGrid.h"
#include "vtkXMLHyperTreeGridWriter.h"
#include "vtkHyperTreeCursor.h"  // not in VTK 8.90
#ifdef LEMMA_VTK9_SUPPORT
#include "vtkHyperTreeGridNonOrientedCursor.h"
#endif
//#include "vtkHyperTreeGridLevelEntry.h" VTK 9
#include "vtkDoubleArray.h"
#endif

namespace Lemma {

    // Holds the elliptic field construction of Bperp
    // commented out variables are for error checking
    struct EllipticB {
        Real      alpha;
        Real      beta;
        Real      zeta;
//        Real      err;
        Complex   eizt;
//        Complex   BperpdotB;
        Vector3r  bhat;
        Vector3r  bhatp;
//        Vector3cr Bperp;
    };

    template <typename T> int sgn(T val) {
        return (val > T(0)) - (val < T(0));
    }

    /**
     * \ingroup Merlin
     * \brief  Calculated the initial amplitude imaging kernel of a sNMR experiment
     * \details This class calculates the imaging kernel for a free induction decay
     *          pulse. The methodology follows from Weichman et al., 2000.
     */
    class KernelV0 : public MerlinObject {

        friend std::ostream &operator<<(std::ostream &stream, const KernelV0 &ob);

        protected:
        /*
         *  This key is used to lock the constructor. It is protected so that inhereted
         *  classes also have the key to contruct their base class.
         */

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

        /**
         *  Factory method for generating concrete class.
         *  @return a std::shared_ptr of type KernelV0
         */
        static std::shared_ptr< KernelV0 > NewSP();

        /**
         *   Constructs an KernelV0 object from a YAML::Node.
         *   @see KernelV0::Serialize
         */
        static std::shared_ptr<KernelV0> DeSerialize(const YAML::Node& node);

        /**
         *   Constructs an object from a string representation of a YAML::Node. This is primarily
         *   used in Python wrapping
         */
        static std::shared_ptr<KernelV0> DeSerialize( const std::string& node ) {
            return KernelV0::DeSerialize(YAML::LoadFile(node));
        }

        // ====================  OPERATORS     =======================

        // ====================  OPERATIONS    =======================

        /**
         *   Calculates a single imaging kernel, however, phased arrays are supported
         *   so that more than one transmitter and/or receiver can be specified.
         *   @param[in] tx is the list of transmitters to use for a kernel, use the same labels as
         *              used in PushCoil.
         *   @param[in] rx is the list of receivers to use for a kernel, use the same labels as
         *              used in PushCoil. @see PushCoil
         *   @param[in] vtkOutput generates a VTK hyperoctree file as well, useful for visualization.
         *              requires compilation of Lemma with VTK. The VTK files can become very large.
         */
        void CalculateK0 (const std::vector< std::string >& tx, const std::vector< std::string >& rx,
                bool vtkOutput );
                //bool vtkOutput=false );

        /**
         *  Aligns the kernel pulse settings with an Akvo Processed dataset.
         */
        void AlignWithAkvoDataset( const YAML::Node& node ) ;

        /**
         *  Aligns the kernel pulse settings with an Akvo Processed dataset.
         */
        void AlignWithAkvoDataset( const std::string& file ) {
            AlignWithAkvoDataset( YAML::LoadFile(file) );
        }

        /**
         *   Assign transmiter coils
         */
        inline void PushCoil( const std::string& label, std::shared_ptr<PolygonalWireAntenna> ant ) {
            TxRx[label] = ant;
        }

        // ====================  INQUIRY       =======================
        /**
         * @return std::shared_ptr<LayeredEarthEM>
         */
        inline std::shared_ptr<LayeredEarthEM> GetSigmaModel (  ) {
            return SigmaModel;
        }		// -----  end of method KernelV0::get_SigmaModel  -----

        /**
         * @return the kernel matrix
         */
        inline MatrixXcr GetKernel ( ) {
            return Kern;
        }

        /**
         * @return the integration tolerance
         */
        inline Real GetTolerance ( ) {
            return tol;
        }

        /**
         * @return the layer interfaces
         */
        inline VectorXr GetInterfaces ( ) {
            return Interfaces;
        }

        /**
         * @return the pulse peak current
         */
        inline VectorXr GetPulseCurrent ( ) {
            return PulseI;
        }

        /**
         * @param[in] value the 1D-EM model used for calculations
         */
        inline void SetLayeredEarthEM ( std::shared_ptr< LayeredEarthEM > value ) {
            SigmaModel	= value;
            return ;
        }		// -----  end of method KernelV0::set_SigmaModel  -----

        /**
         *  @param[in] size the size of the volume to be integrated
         */
        inline void SetIntegrationSize ( const Vector3r& size ) {
            Size = size;
            return ;
        }		// -----  end of method KernelV0::SetIntegrationSize  -----

        /**
         * @param[in] type The type of Hankel transform that will be used.
         */
        inline void SetHankelTransformType ( const HANKELTRANSFORMTYPE& type ) {
            HankelType = type;
            return ;
        }		// -----  end of method KernelV0::SetIntegrationOrigin  -----


        /**
         * @param[in] min is the minimum leaf level, defaults to 0
         */
        inline void SetMinLevel ( const int& min ) {
            minLevel = min;
            return ;
        }		// -----  end of method KernelV0::SetMinLevel  -----

        /**
         * @param[in] max is the maximum leaf level, defaults to 12
         */
        inline void SetMaxLevel ( const int& max ) {
            maxLevel = max;
            return ;
        }		// -----  end of method KernelV0::SetMaxLevel  -----

        /**
         *  @param[in] origin The origin location (corner) for the integration volume
         */
        inline void SetIntegrationOrigin ( const Vector3r& origin ) {
            Origin = origin;
            return ;
        }		// -----  end of method KernelV0::SetIntegrationOrigin  -----

        /**
         * @param[in] Amps is the current for each pulse moment
         */
        inline void SetPulseCurrent ( const VectorXr& Amps ) {
            PulseI = Amps;
            return ;
        }		// -----  end of method KernelV0::SetIntegrationOrigin  -----

        /**
         *  Sets the temperature, which has implications in calculation of \f$ M_N^{(0)}\f$. Units in
         *  Kelvin.
         */
        inline void SetTemperature(const Real& tempK) {
            Temperature = tempK;
        }

        /**
         *  Sets the tolerance to use for making the adaptive mesh
         *  @param[in] ttol is the tolerance to use
         */
        inline void SetTolerance(const Real& ttol) {
            tol = ttol;
        }

        /**
         *  @param[in] taup sets the pulse duration
         */
        inline void SetPulseDuration(const Real& taup) {
            Taup = taup;
        }

        inline Real GetPulseDuration( ) {
            return Taup;
        }

        inline void SetDepthLayerInterfaces( const VectorXr& iface ){
            Interfaces = iface;
        }

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
        VectorXcr f( const Vector3r& r, const Real& volume , const Vector3cr& Ht, const Vector3cr& Hr);

//         Complex ComputeV0Cell(const EllipticB& EBT, const EllipticB& EBR,
//                 const Real& sintheta, const Real& phase, const Real& Mn0Abs,
//                 const Real& vol);

        EllipticB EllipticFieldRep (const Vector3cr& B, const Vector3r& B0hat);

        Vector3r ComputeMn0(const Real& Porosity, const Vector3r& B0);

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
                            const VectorXcr& parentVal );

        #ifdef LEMMAUSEVTK
        /**
         *  Same functionality as @see EvaluateKids, but includes generation of a VTK
         *  HyperOctree, which is useful for visualization.
         */
        #ifdef LEMMA_VTK8_SUPPORT
        void EvaluateKids2(  const Vector3r& size, const int& level, const Vector3r& cpos,
                            const VectorXcr& parentVal, vtkHyperTreeGrid* octree, vtkHyperTreeCursor* curse );
        #elif LEMMA_VTK9_SUPPORT
        void EvaluateKids2(  const Vector3r& size, const int& level, const Vector3r& cpos,
                            const VectorXcr& parentVal, vtkHyperTree* octree, vtkHyperTreeGridNonOrientedCursor* curse );
        #endif

        void GetPosition( vtkHyperTreeCursor* Cursor, Real* p );
        #endif

        // ====================  DATA MEMBERS  =========================

        int                                       ilay;
        int                                       nleaves;
        int                                       minLevel=0;
        int                                       maxLevel=12;

        Real                                      VOLSUM;
        Real                                      tol=1e-11;
        Real                                      Temperature=283.;
        Real                                      Taup = .020;  // Sec
        Real                                      Larmor;

        Vector3r                                  Size;
        Vector3r                                  Origin;

        VectorXr                                  PulseI;
        VectorXr                                  Interfaces;

        MatrixXcr                                 Kern;

        HANKELTRANSFORMTYPE                       HankelType=ANDERSON801;

        std::shared_ptr< LayeredEarthEM >         SigmaModel = nullptr;
        std::shared_ptr< FieldPoints >            cpoints = nullptr;

        std::map< std::string , std::shared_ptr< PolygonalWireAntenna > >  TxRx;
        std::map< std::string , std::shared_ptr< EMEarth1D > >             EMEarths;

        #ifdef LEMMAUSEVTK
        std::map< int, VectorXcr >                LeafDict;      // kernel sum for each q
        std::map< int, VectorXcr >                LeafHt;        // Transmitter field
        std::map< int, VectorXcr >                LeafHr;        // Receiver field
        std::map< int, int >                      LeafDictIdx;   // index
        std::map< int, Real >                     LeafDictErr;   // error value
        #endif

        ProgressBar* disp;
        int percent_done;

        // Physical constants and conversion factors
        static constexpr Real GAMMA = 2.67518e8;                  // MKS units
        static constexpr Real INVSQRT2 = 0.70710678118654746;     // 1/sqrt(2)
        static constexpr Real HBAR = 1.05457148e-34;              // m2 kg / s
        static constexpr Real NH2O = 6.692e28;                    // [m^3]
        static constexpr Real KB = 1.3805e-23;                    // m^2 kg s-2 K-1
        static constexpr Real CHI_N = 3.29e-3;                    // MKS units

        /** ASCII string representation of the class name */
        static constexpr auto CName = "KernelV0";

    }; // -----  end of class  KernelV0  -----
}  // -----  end of namespace Lemma ----

/* vim: set tabstop=4 expandtab */
/* vim: set filetype=cpp */

#endif   // ----- #ifndef KERNELV0_INC  -----
