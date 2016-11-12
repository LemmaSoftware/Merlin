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
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2016, University of Utah
 * @copyright Copyright (c) 2016, Lemma Software, LLC
 */


#include "KernelV0.h"
#include "EMEarth1D.h"
#include "FieldPoints.h"

namespace Lemma {

    // ====================  FRIEND METHODS  =====================

    std::ostream &operator << (std::ostream &stream, const KernelV0 &ob) {
        stream << ob.Serialize()  << "\n---\n"; // End of doc ---
        return stream;
    }

    // ====================  LIFECYCLE     =======================

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  KernelV0
    // Description:  constructor (locked)
    //--------------------------------------------------------------------------------------
    KernelV0::KernelV0 (const ctor_key&) : LemmaObject( ) {

    }  // -----  end of method KernelV0::KernelV0  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  KernelV0
    // Description:  DeSerializing constructor (locked)
    //--------------------------------------------------------------------------------------
    KernelV0::KernelV0 (const YAML::Node& node, const ctor_key&) : LemmaObject(node) {

    }  // -----  end of method KernelV0::KernelV0  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  NewSP()
    // Description:  public constructor returing a shared_ptr
    //--------------------------------------------------------------------------------------
    std::shared_ptr< KernelV0 >  KernelV0::NewSP() {
        return std::make_shared< KernelV0 >( ctor_key() );
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  ~KernelV0
    // Description:  destructor (protected)
    //--------------------------------------------------------------------------------------
    KernelV0::~KernelV0 () {

    }  // -----  end of method KernelV0::~KernelV0  (destructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  Serialize
    //--------------------------------------------------------------------------------------
    YAML::Node  KernelV0::Serialize (  ) const {
        YAML::Node node = LemmaObject::Serialize();
        node.SetTag( GetName() );

        // Coils Transmitters & Receivers
        for ( auto txm : TxRx) {
            node[txm.first] = txm.second->Serialize();
        }

        // LayeredEarthEM
        node["SigmaModel"] = SigmaModel->Serialize();

        return node;
    }		// -----  end of method KernelV0::Serialize  -----

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  DeSerialize
    //--------------------------------------------------------------------------------------
    std::shared_ptr<KernelV0> KernelV0::DeSerialize ( const YAML::Node& node  ) {
        if (node.Tag() !=  "KernelV0" ) {
            throw  DeSerializeTypeMismatch( "KernelV0", node.Tag());
        }
        return std::make_shared< KernelV0 > ( node, ctor_key() );
    }		// -----  end of method KernelV0::DeSerialize  -----

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  DeSerialize
    //--------------------------------------------------------------------------------------
    void KernelV0::CalculateK0 (const std::vector< std::string>& Tx, const std::vector<std::string >& Rx ) {

        for (auto tx : Tx) {
            // Set up EMEarth
            auto EmEarth = EMEarth1D::NewSP();
                EmEarth->AttachWireAntenna(TxRx[tx]);
                EmEarth->AttachLayeredEarthEM(SigmaModel);
         		EmEarth->SetFieldsToCalculate(H);
                // TODO query for method, altough with flat antennae, this is fastest
                EmEarth->SetHankelTransformMethod(ANDERSON801);

// 		EmEarth->AttachFieldPoints(receivers);
//         //EmEarth->SetHankelTransformMethod(FHTKEY101);
// 	    EmEarth->CalculateWireAntennaFields();
//         Vector3Xcr Rx1 = receivers->GetHfield(0);
//         //receivers->ClearFields();
//
// 		//EmEarth->AttachWireAntenna(Tx2);
// 	    EmEarth->CalculateWireAntennaFields();
//         Rx1 += receivers->GetHfield(0);

        }

    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  IntegrateOnOctreeGrid
    //--------------------------------------------------------------------------------------
    void KernelV0::IntegrateOnOctreeGrid( const Real& tolerance) {

        Vector3r                Size;
        Vector3r                Origin;
        Vector3r                step;
        Vector3r                cpos;

        int                     level;
        int                     maxlevel;
        int                     index;
        int                     counter;

        Real                    cvol;
        Real                    tvol;
        Real                    tol;
        Complex                 KernelSum;

        //this->tol = tolerance;
        Real KernelSum = 0.;
        //Cursor->ToRoot();
        //Cubes->SetNumberOfReceivers(8);
        EvaluateKids( 1e9 ); // Large initial number don't waste time actually computing
        //EvaluateKids();
        //std::cout << "Kernel Sum from Generate Mesh "
        //    << std::real(KernelSum) << "\t" << std::imag(KernelSum) << std::endl;

        // old VTK thingy
        //SetLeafDataFromGridCreation();
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  EvaluateKids
    //--------------------------------------------------------------------------------------
    void KernelV0::EvaluateKids(const Complex& kval) {

        assert("Evaluate Kids pre" && Cursor->CurrentIsLeaf());
        vtkHyperOctreeCursor *tcurse = Cursor->Clone();
        Real p[3];
        Octree->SubdivideLeaf(Cursor);
        tcurse->ToSameNode(Cursor);
        std::cout << "\rPredivide Leaf count: " << Octree->GetNumberOfLeaves();

        //std::cout.flush();
        for (int child=0; child<8; ++child) {
            Cursor->ToChild(child);
            assert(Cursor->CurrentIsLeaf());
            // Build cube
            GetPosition(p);
            cpos <<  p[0], p[1], p[2];
            step  = ((Size).array() / std::pow(2.,Cursor->GetCurrentLevel()));
            Cubes->SetLocation(child, cpos);
            Cubes->SetLength(child, step);
            //std::cout << "child " << child << " cpos\t" << cpos.transpose() << std::endl;
            //std::cout << "child " << child << " step\t" << step.transpose() << std::endl;
            Cursor->ToSameNode(tcurse);
        }

        // make calculation
        Cubes->ClearFields();
        VectorXcr f = SenseKernel->ComputeSensitivity();
        if ( std::abs(std::abs(kval) - std::abs(f.array().abs().sum())) <= tol ||
            Cursor->GetCurrentLevel() >= maxlevel ) {
    	    // stop subdividing, save result
    	    for (int child=0; child < 8; ++ child) {
    	        Cursor->ToChild(child);
    	        leafdata.push_back( std::abs(f(child)) / Cubes->GetVolume(child) );
    	        // TODO fval is just a test
    	        //leafdata.push_back( fval );
    	        leafids.push_back(Cursor->GetLeafId());
    	        KernelSum += f(child);
    	        Cursor->ToParent();
            }
    	    return;
        }

        //std::cout << std::abs(kval) << "\t" <<
        //         std::abs(f.array().abs().sum()) << "\t" << tol << std::endl;
        for (int child=0; child < 8; ++ child) {
            //std::cout << "Down the rabit hole " <<std::endl;
            Cursor->ToChild(child);
            EvaluateKids( f(child) );
            //Cursor->ToParent();
            Cursor->ToSameNode(tcurse);
        }
        tcurse->Delete();
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  EvaluateKids
    //--------------------------------------------------------------------------------------
    void OctreeGrid::GetPosition( Real* p ) {
        Real ratio=1.0/(1<<(Cursor->GetCurrentLevel()));
        //step  = ((Size).array() / std::pow(2.,Cursor->GetCurrentLevel()));
        p[0]=(Cursor->GetIndex(0)+.5)*ratio*Size[0]+Origin[0] ;//+ .5*step[0];
        p[1]=(Cursor->GetIndex(1)+.5)*ratio*Size[1]+Origin[1] ;//+ .5*step[1];
        p[2]=(Cursor->GetIndex(2)+.5)*ratio*Size[2]+Origin[2] ;//+ .5*step[2];
    }

} // ----  end of namespace Lemma  ----

/* vim: set tabstop=4 expandtab */
/* vim: set filetype=cpp */

