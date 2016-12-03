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
#ifndef  LOOPINTERACTIONS
#define  LOOPINTERACTIONS

#pragma once
#include "LemmaObject.h"
#include "LayeredEarthEM.h"
#include "PolygonalWireAntenna.h"
#include "EMEarth1D.h"
#include "FieldPoints.h"

#ifdef LEMMAUSEVTK
#include "vtkHyperOctree.h"
#include "vtkHyperOctreeCursor.h"
#include "vtkXMLHyperOctreeWriter.h"
#include "vtkDoubleArray.h"
#endif

namespace Lemma {

    enum INTERACTION {COUPLING, INTERFERENCE, PHASE};
    /// convert enums to string saves repeated code useful for YAML serializing
    std::string enum2String(const INTERACTION& type);

    /**
     * \ingroup Merlin
     * \brief
     * \details
     */
    template< INTERACTION Type  >
    class LoopInteractions : public LemmaObject {

        friend std::ostream &operator << (std::ostream &stream, const LoopInteractions &ob) {
            stream << ob.Serialize()  << "\n---\n"; // End of doc ---
            return stream;
        }

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
         * @see LoopInteractions::NewSP
         */
        explicit LoopInteractions ( const ctor_key& ) : LemmaObject () { }

        /**
         * DeSerializing constructor.
         * @note This method is locked, and cannot be called directly.
         *       The reason that the method is public is to enable the use
         *       of make_shared whilst enforcing the use of shared_ptr,
         *       in c++-17, this curiosity may be resolved.
         * @see LoopInteractions::DeSerialize
         */
        LoopInteractions ( const YAML::Node& node, const ctor_key& ) : LemmaObject(node) { }

        /**
         * Default destructor.
         * @note This method should never be called due to the mandated
         *       use of smart pointers. It is necessary to keep the method
         *       public in order to allow for the use of the more efficient
         *       make_shared constructor.
         */
        virtual ~LoopInteractions () { }

        /**
         *  Uses YAML to serialize this object.
         *  @return a YAML::Node
         *  @see LoopInteractions::DeSerialize
         */
        virtual YAML::Node Serialize() const {
            YAML::Node node = LemmaObject::Serialize();
            node.SetTag( GetName() );

            // Coils Transmitters & Receivers
            for ( auto txm : TxRx) {
                node[txm.first] = txm.second->Serialize();
            }
            // LayeredEarthEM
            node["SigmaModel"] = SigmaModel->Serialize();

            node["tol"] = tol;
            node["minLevel"] = minLevel;
            node["maxLevel"] = maxLevel;

            return node;
        }

        /*
         *  Factory method for generating concrete class.
         *  @return a std::shared_ptr of type LoopInteractions
         */
        static std::shared_ptr< LoopInteractions > NewSP() {
            return std::make_shared< LoopInteractions >( ctor_key() );
        }

        /**
         *   Constructs an LoopInteractions object from a YAML::Node.
         *   @see LoopInteractions::Serialize
         */
        static std::shared_ptr<LoopInteractions> DeSerialize(const YAML::Node& node) {
            if (node.Tag() !=  "LoopInteractions" ) {
                throw  DeSerializeTypeMismatch( "LoopInteractions", node.Tag());
            }
            return std::make_shared< LoopInteractions > ( node, ctor_key() );
        }		// -----  end of method LoopInteractions::DeSerialize  -----

        // ====================  OPERATORS     =======================

        // ====================  OPERATIONS    =======================

        /**
         * @return std::shared_ptr<LayeredEarthEM>
         */
        inline std::shared_ptr<LayeredEarthEM> GetSigmaModel (  ) {
            return SigmaModel;
        }		// -----  end of method LoopInteractions::get_SigmaModel  -----

        /**
         * @param[in] value the 1D-EM model used for calculations
         */
        inline void SetLayeredEarthEM ( std::shared_ptr< LayeredEarthEM > value ) {
            SigmaModel	= value;
            return ;
        }		// -----  end of method LoopInteractions::set_SigmaModel  -----

        /**
         *
         */
        inline void SetIntegrationSize ( const Vector3r& size ) {
            Size = size;
            return ;
        }		// -----  end of method LoopInteractions::SetIntegrationSize  -----

        /**
         *
         */
        inline void SetIntegrationOrigin ( const Vector3r& origin ) {
            Origin = origin;
            return ;
        }		// -----  end of method LoopInteractions::SetIntegrationOrigin  -----

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
        LoopInteractions( const LoopInteractions& ) = delete;

        private:

        /**
         *  Returns the kernel value for an input prism
         */
        virtual Complex f( const Vector3r& r, const Real& volume , const Vector3cr& Ht, const Vector3cr& Hr);

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

        #ifdef LEMMAUSEVTK
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
        static constexpr auto CName = "LoopInteractions";

    }; // -----  end of class  LoopInteractions  -----

    ///////////////////////////////////////////////////////////////
    // Implimentation of non specialized args -- templated class //
    ///////////////////////////////////////////////////////////////

    // forward declare specs

    template <>
    Complex LoopInteractions<COUPLING>::f( const Vector3r& r, const Real& volume, const Vector3cr& Ht, const Vector3cr& Hr );

    template <>
    Complex LoopInteractions<INTERFERENCE>::f( const Vector3r& r, const Real& volume, const Vector3cr& Ht, const Vector3cr& Hr );

    template <>
    Complex LoopInteractions<PHASE>::f( const Vector3r& r, const Real& volume, const Vector3cr& Ht, const Vector3cr& Hr );

    //--------------------------------------------------------------------------------------
    //       Class:  LoopInteractions
    //      Method:  DeSerialize
    //--------------------------------------------------------------------------------------
    template< INTERACTION Type  >
    Complex LoopInteractions<Type>::Calculate (const std::vector< std::string>& Tx, const std::vector<std::string >& Rx,
            bool vtkOutput ) {

        static bool first = false; // a little hackish
        if (!first) {
            // All EM calculations will share same field points
            cpoints = FieldPoints::NewSP();
                cpoints->SetNumberOfPoints(8);
        }
        first = true;

        for (auto tx : Tx) {
            // Set up EMEarth
            EMEarths[tx] = EMEarth1D::NewSP();
                EMEarths[tx]->AttachWireAntenna(TxRx[tx]);
                EMEarths[tx]->AttachLayeredEarthEM(SigmaModel);
                EMEarths[tx]->AttachFieldPoints( cpoints );
         		EMEarths[tx]->SetFieldsToCalculate(H);
                // TODO query for method, altough with flat antennae, this is fastest
                EMEarths[tx]->SetHankelTransformMethod(ANDERSON801);
                EMEarths[tx]->SetTxRxMode(TX);
                TxRx[tx]->SetCurrent(1.);
        }
        for (auto rx : Rx) {
            if (EMEarths.count(rx)) {
                EMEarths[rx]->SetTxRxMode(TXRX);
            } else {
                EMEarths[rx] = EMEarth1D::NewSP();
                    EMEarths[rx]->AttachWireAntenna(TxRx[rx]);
                    EMEarths[rx]->AttachLayeredEarthEM(SigmaModel);
                    EMEarths[rx]->AttachFieldPoints( cpoints );
         		    EMEarths[rx]->SetFieldsToCalculate(H);
                    // TODO query for method, altough with flat antennae, this is fastest
                    EMEarths[rx]->SetHankelTransformMethod(ANDERSON801);
                    EMEarths[rx]->SetTxRxMode(RX);
                    TxRx[rx]->SetCurrent(1.);
            }
        }
        SUM = 0;
        IntegrateOnOctreeGrid( vtkOutput );
        std::cout << "\nFinished KERNEL\n";
        EMEarths.clear();
        return SUM;
    }


    //--------------------------------------------------------------------------------------
    //       Class:  LoopInteractions
    //      Method:  IntegrateOnOctreeGrid
    //--------------------------------------------------------------------------------------
    template< INTERACTION Type  >
    void LoopInteractions<Type>::IntegrateOnOctreeGrid( bool vtkOutput ) {

        static int count = 0;

        Vector3r cpos = Origin + Size/2.;

        VOLSUM = 0;
        nleaves = 0;
        if (!vtkOutput) {
            EvaluateKids( Size, 0, cpos, Complex(100.));
        } else {
        #ifdef LEMMAUSEVTK
            vtkHyperOctree* oct = vtkHyperOctree::New();
                oct->SetDimension(3);
                oct->SetOrigin( Origin(0), Origin(1), Origin(2) );
                oct->SetSize( Size(0), Size(1), Size(2) );
            vtkHyperOctreeCursor* curse = oct->NewCellCursor();
                curse->ToRoot();
            EvaluateKids2( Size, 0, cpos, Complex(100.0), oct, curse );

            // Fill in leaf data
            vtkDoubleArray* kr = vtkDoubleArray::New();
                kr->SetNumberOfComponents(1);
                kr->SetName("Re($\\sum$)");
                kr->SetNumberOfTuples( oct->GetNumberOfLeaves() );
            vtkDoubleArray* ki = vtkDoubleArray::New();
                ki->SetNumberOfComponents(1);
                ki->SetName("Im($\\sum$)");
                ki->SetNumberOfTuples( oct->GetNumberOfLeaves() );
            vtkDoubleArray* km = vtkDoubleArray::New();
                km->SetNumberOfComponents(1);
                km->SetName("mod($\\sum$)");
                km->SetNumberOfTuples( oct->GetNumberOfLeaves() );
            vtkIntArray* kid = vtkIntArray::New();
                kid->SetNumberOfComponents(1);
                kid->SetName("ID");
                kid->SetNumberOfTuples( oct->GetNumberOfLeaves() );
            vtkIntArray* kerr = vtkIntArray::New();
                kerr->SetNumberOfComponents(1);
                kerr->SetName("nleaf");

            //Real LeafVol(0);
            for (auto leaf : LeafDict) {
                kr->InsertTuple1( leaf.first, std::real(leaf.second) );
                ki->InsertTuple1( leaf.first, std::imag(leaf.second) );
                km->InsertTuple1( leaf.first, std::abs(leaf.second) );
                kid->InsertTuple1( leaf.first, leaf.first );
                //LeafVol += std::real(leaf.second);
            }
            //std::cout << "\n\nLeafVol=" << LeafVol << std::endl;

            for (auto leaf : LeafDictIdx) {
                kerr->InsertTuple1( leaf.first, leaf.second );
            }

            auto kri = oct->GetLeafData()->AddArray(kr);
            auto kii = oct->GetLeafData()->AddArray(ki);
            auto kmi = oct->GetLeafData()->AddArray(km);
            auto kidi = oct->GetLeafData()->AddArray(kid);
            auto keri = oct->GetLeafData()->AddArray(kerr);

            auto write = vtkXMLHyperOctreeWriter::New();
                //write.SetDataModeToAscii()
                write->SetInputData(oct);
                std::string fname = std::string("octree-") + enum2String(Type) + std::string("-")
                                    + to_string(count) + std::string(".vto");
                write->SetFileName(fname.c_str());
                write->Write();
                write->Delete();

            oct->GetLeafData()->RemoveArray( kri );
            oct->GetLeafData()->RemoveArray( kii );
            oct->GetLeafData()->RemoveArray( kmi );
            oct->GetLeafData()->RemoveArray( kidi );
            oct->GetLeafData()->RemoveArray( keri );

            kerr->Delete();
            kid->Delete();
            kr->Delete();
            ki->Delete();
            km->Delete();

            curse->Delete();
            oct->Delete();
        #else
            throw std::runtime_error("IntegrateOnOctreeGrid with vtkOutput requires Lemma with VTK support");
        #endif

        }
        std::cout << "\nVOLSUM=" << VOLSUM << "\tActual=" <<  Size(0)*Size(1)*Size(2)
                  << "\tDifference=" << VOLSUM - (Size(0)*Size(1)*Size(2)) <<  std::endl;
        count += 1;
    }

    //--------------------------------------------------------------------------------------
    //       Class:  LoopInteractions
    //      Method:  EvaluateKids
    //--------------------------------------------------------------------------------------
    template< INTERACTION Type  >
    void LoopInteractions<Type>::EvaluateKids( const Vector3r& size, const int& level, const Vector3r& cpos,
        const Complex& parentVal ) {

        std::cout << "\r" << (int)(1e2*VOLSUM/(Size[0]*Size[1]*Size[2])) << "\t" << nleaves;
        std::cout.flush();

        // Next level step, interested in one level below
        // bitshift requires one extra, faster than, and equivalent to std::pow(2, level+1)
        Vector3r step  = size.array() / (Real)(1 << (level+1) );
        Real vol = (step(0)*step(1)*step(2));     // volume of each child

        Vector3r pos =  cpos - step/2.;
        Eigen::Matrix<Real, 8, 3> posadd = (Eigen::Matrix<Real, 8, 3>() <<
                        0,       0,       0,
                  step[0],       0,       0,
                        0, step[1],       0,
                  step[0], step[1],       0,
                        0,       0, step[2],
                  step[0],       0, step[2],
                        0, step[1], step[2],
                  step[0], step[1], step[2] ).finished();

        VectorXcr kvals(8);       // individual kernel vals
        cpoints->ClearFields();
        for (int ichild=0; ichild<8; ++ichild) {
            Vector3r cp = pos;    // Eigen complains about combining these
            cp += posadd.row(ichild);
            cpoints->SetLocation( ichild, cp );
        }

        Eigen::Matrix<Complex, 3, 8> Ht = Eigen::Matrix<Complex, 3, 8>::Zero();
        Eigen::Matrix<Complex, 3, 8> Hr = Eigen::Matrix<Complex, 3, 8>::Zero();
        for ( auto EMCalc : EMEarths ) {
            EMCalc.second->GetFieldPoints()->ClearFields();
            EMCalc.second->CalculateWireAntennaFields();
            switch (EMCalc.second->GetTxRxMode()) {
                case TX:
                    Ht += EMCalc.second->GetFieldPoints()->GetHfield(0);
                    break;
                case RX:
                    Hr += EMCalc.second->GetFieldPoints()->GetHfield(0);
                    break;
                case TXRX:
                    Ht += EMCalc.second->GetFieldPoints()->GetHfield(0);
                    Hr += EMCalc.second->GetFieldPoints()->GetHfield(0);
                    break;
                default:
                    break;
            }
        }

        for (int ichild=0; ichild<8; ++ichild) {
            Vector3r cp = pos;    // Eigen complains about combining these
            cp += posadd.row(ichild);
            kvals(ichild) = f(cp, vol, Ht.col(ichild), Hr.col(ichild));
        }

        Complex ksum = kvals.sum();     // Kernel sum
        // Evaluate whether or not furthur splitting is needed
        if ( (std::abs(ksum-parentVal) > tol && level < maxLevel) || level < minLevel ) {
            // Not a leaf dive further in
            for (int ichild=0; ichild<8; ++ichild) {
                Vector3r cp = pos; // Eigen complains about combining these
                cp += posadd.row(ichild);
                EvaluateKids( size, level+1, cp, kvals(ichild) );
            }
            return; // not leaf
        }
        // implicit else, is a leaf
        SUM += ksum;
        VOLSUM  += 8.*vol;
        nleaves += 1; // could say += 8 just as fairly
        return;     // is leaf
    }

    #ifdef LEMMAUSEVTK
    //--------------------------------------------------------------------------------------
    //       Class:  LoopInteractions
    //      Method:  EvaluateKids2 -- same as Evaluate Kids, but include VTK octree generation
    //--------------------------------------------------------------------------------------
    template< INTERACTION Type  >
    void LoopInteractions<Type>::EvaluateKids2( const Vector3r& size, const int& level, const Vector3r& cpos,
        const Complex& parentVal, vtkHyperOctree* oct, vtkHyperOctreeCursor* curse) {

        std::cout << "\r" << (int)(1e2*VOLSUM/(Size[0]*Size[1]*Size[2])) << "\t" << nleaves;
        std::cout.flush();

        // Next level step, interested in one level below
        // bitshift requires one extra, faster than, and equivalent to std::pow(2, level+1)
        Vector3r step  = size.array() / (Real)(1 << (level+1) );
        Real vol = (step(0)*step(1)*step(2));         // volume of each child

        Vector3r pos =  cpos - step/2.;
        Eigen::Matrix<Real, 8, 3> posadd = (Eigen::Matrix<Real, 8, 3>() <<
                        0,       0,       0,
                  step[0],       0,       0,
                        0, step[1],       0,
                  step[0], step[1],       0,
                        0,       0, step[2],
                  step[0],       0, step[2],
                        0, step[1], step[2],
                  step[0], step[1], step[2] ).finished();

        VectorXcr kvals(8);       // individual kernel vals
        cpoints->ClearFields();
        for (int ichild=0; ichild<8; ++ichild) {
            Vector3r cp = pos;    // Eigen complains about combining these
            cp += posadd.row(ichild);
            cpoints->SetLocation( ichild, cp );
        }

        Eigen::Matrix<Complex, 3, 8> Ht = Eigen::Matrix<Complex, 3, 8>::Zero();
        Eigen::Matrix<Complex, 3, 8> Hr = Eigen::Matrix<Complex, 3, 8>::Zero();
        for ( auto EMCalc : EMEarths ) {
            //EMCalc->GetFieldPoints()->ClearFields();
            EMCalc.second->CalculateWireAntennaFields();
            switch (EMCalc.second->GetTxRxMode()) {
                case TX:
                    Ht += EMCalc.second->GetFieldPoints()->GetHfield(0);
                    break;
                case RX:
                    Hr += EMCalc.second->GetFieldPoints()->GetHfield(0);
                    break;
                case TXRX:
                    Ht += EMCalc.second->GetFieldPoints()->GetHfield(0);
                    Hr += EMCalc.second->GetFieldPoints()->GetHfield(0);
                    break;
                default:
                    break;
            }
        }

        for (int ichild=0; ichild<8; ++ichild) {
            Vector3r cp = pos;    // Eigen complains about combining these
            cp += posadd.row(ichild);
            kvals(ichild) = f(cp, vol, Ht.col(ichild), Hr.col(ichild));
        }

        Complex ksum = kvals.sum();     // Kernel sum
        // Evaluate whether or not furthur splitting is needed
        if ( (std::abs(ksum-parentVal) > tol && level < maxLevel) || level < minLevel ) {
            oct->SubdivideLeaf(curse);
            for (int ichild=0; ichild<8; ++ichild) {
                curse->ToChild(ichild);
                Vector3r cp = pos; // Eigen complains about combining these
                cp += posadd.row(ichild);
                /* Test for position via alternative means */
                /*
                Real p[3];
                GetPosition(curse, p);
                if ( (Vector3r(p) - cp).norm() > 1e-8 ) {
                    std::cout << "ERROR @ nleaves" << nleaves << "\n" << cp[0] << "\t" << p[0] << "\t" << cp[1] << "\t" << p[1]
                              << "\t" << cp[2] << "\t" << p[2] << "\t" << vol<< std::endl;
                    throw std::runtime_error("doom");
                }
                */
                /* End of position test */
                EvaluateKids2( size, level+1, cp, kvals(ichild), oct, curse );
                curse->ToParent();
            }
            return;  // not a leaf
        }
        LeafDict[curse->GetLeafId()] = ksum/(8.*vol);
        LeafDictIdx[curse->GetLeafId()] = nleaves;
        SUM += ksum;
        VOLSUM += 8*vol;
        nleaves += 1;
        return;     // is a leaf
    }

    //--------------------------------------------------------------------------------------
    //       Class:  LoopInteractions
    //      Method:  GetPosition
    //--------------------------------------------------------------------------------------
    template< INTERACTION Type  >
    void LoopInteractions<Type>::GetPosition( vtkHyperOctreeCursor* Cursor, Real* p ) {
        Real ratio=1.0/(1<<(Cursor->GetCurrentLevel()));
        //step  = ((Size).array() / std::pow(2.,Cursor->GetCurrentLevel()));
        p[0]=(Cursor->GetIndex(0)+.5)*ratio*this->Size[0]+this->Origin[0] ;//+ .5*step[0];
        p[1]=(Cursor->GetIndex(1)+.5)*ratio*this->Size[1]+this->Origin[1] ;//+ .5*step[1];
        p[2]=(Cursor->GetIndex(2)+.5)*ratio*this->Size[2]+this->Origin[2] ;//+ .5*step[2];
    }
    #endif

}  // -----  end of namespace Lemma ----

/* vim: set tabstop=4 expandtab */
/* vim: set filetype=cpp */

#endif   // ----- #ifndef LOOPINTERACTIONS  -----
