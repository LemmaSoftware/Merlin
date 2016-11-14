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

        IntegrateOnOctreeGrid( 1e-2 );

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

        this->tol = tolerance;
        Vector3r                Size;
            Size << 100,100,100;
        //Vector3r                Origin;
        Vector3r                cpos;
            cpos << 50,50,50;
        int                     maxlevel;

        SUM = 0;
        nleaves = 0;
        EvaluateKids( Size, 0, cpos, -1e2 );
        std::cout << "SUM\t" << SUM << "\t" << 100*100*100 << "\t" << SUM - Complex(100.*100.*100.) <<  std::endl;
        std::cout << "nleaves\t" << nleaves << std::endl;

    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  f
    //--------------------------------------------------------------------------------------
    Complex KernelV0::f( const Vector3r& r, const Real& volume ) {
        return Complex(volume);
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  EvaluateKids
    //--------------------------------------------------------------------------------------
    bool KernelV0::EvaluateKids( const Vector3r& size, const int& level, const Vector3r& cpos,
        const Complex& parentVal ) {

        // Next level step, interested in one level below
        // bitshift requires one extra, faster than, and equivalent to std::pow(2, level+1)
        Vector3r step = size.array() / (Real)(1 << (level+2) );

        Real vol = step(0)*step(1)*step(2);     // volume of each child

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

        VectorXcr kvals(8);                     // individual kernel vals
        for (int ichild=0; ichild<8; ++ichild) {
            Vector3r cp = pos; // Eigen complains about combining these
            cp += posadd.row(ichild);
            kvals(ichild) = f(cp, vol);
        }
        Complex ksum = kvals.sum();     // Kernel sum

        // Evaluate whether or not furthur splitting is needed
        if ( std::abs(ksum - parentVal) > tol || level < 5 ) {
            for (int ichild=0; ichild<8; ++ichild) {
                Vector3r cp = pos; // Eigen complains about combining these
                cp += posadd.row(ichild);
                bool isleaf = EvaluateKids( size, level+1, pos, kvals(ichild) );
                if (isleaf) {  // Include result in final integral
//                  Id = curse.GetLeafId()     // VTK
//                  LeafDict[Id] = vals[child] // VTK
                    SUM += ksum;
                    nleaves += 1;
                }
            }
            return false;  // not leaf
        }
        return true;       // leaf
    }

} // ----  end of namespace Lemma  ----

/* vim: set tabstop=4 expandtab */
/* vim: set filetype=cpp */

