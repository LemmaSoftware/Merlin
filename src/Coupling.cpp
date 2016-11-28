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
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2016, University of Utah
 * @copyright Copyright (c) 2016, Lemma Software, LLC
 * @copyright Copyright (c) 2008, Colorado School of Mines
 */


#include "Coupling.h"
#include "FieldPoints.h"

namespace Lemma {

    // ====================  FRIEND METHODS  =====================

    std::ostream &operator << (std::ostream &stream, const Coupling &ob) {
        stream << ob.Serialize()  << "\n---\n"; // End of doc ---
        return stream;
    }

    // ====================  LIFECYCLE     =======================

    //--------------------------------------------------------------------------------------
    //       Class:  Coupling
    //      Method:  Coupling
    // Description:  constructor (locked)
    //--------------------------------------------------------------------------------------
    Coupling::Coupling (const ctor_key&) : LemmaObject( ) {

    }  // -----  end of method Coupling::Coupling  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  Coupling
    //      Method:  Coupling
    // Description:  DeSerializing constructor (locked)
    //--------------------------------------------------------------------------------------
    Coupling::Coupling (const YAML::Node& node, const ctor_key&) : LemmaObject(node) {

    }  // -----  end of method Coupling::Coupling  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  Coupling
    //      Method:  NewSP()
    // Description:  public constructor returing a shared_ptr
    //--------------------------------------------------------------------------------------
    std::shared_ptr< Coupling >  Coupling::NewSP() {
        return std::make_shared< Coupling >( ctor_key() );
    }

    //--------------------------------------------------------------------------------------
    //       Class:  Coupling
    //      Method:  ~Coupling
    // Description:  destructor (protected)
    //--------------------------------------------------------------------------------------
    Coupling::~Coupling () {

    }  // -----  end of method Coupling::~Coupling  (destructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  Coupling
    //      Method:  Serialize
    //--------------------------------------------------------------------------------------
    YAML::Node  Coupling::Serialize (  ) const {
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
    }		// -----  end of method Coupling::Serialize  -----

    //--------------------------------------------------------------------------------------
    //       Class:  Coupling
    //      Method:  DeSerialize
    //--------------------------------------------------------------------------------------
    std::shared_ptr<Coupling> Coupling::DeSerialize ( const YAML::Node& node  ) {
        if (node.Tag() !=  "Coupling" ) {
            throw  DeSerializeTypeMismatch( "Coupling", node.Tag());
        }
        return std::make_shared< Coupling > ( node, ctor_key() );
    }		// -----  end of method Coupling::DeSerialize  -----

    //--------------------------------------------------------------------------------------
    //       Class:  Coupling
    //      Method:  DeSerialize
    //--------------------------------------------------------------------------------------
    Complex Coupling::Calculate (const std::vector< std::string>& Tx, const std::vector<std::string >& Rx,
            bool vtkOutput ) {
        // All EM calculations will share same field points
        cpoints = FieldPoints::NewSP();
            cpoints->SetNumberOfPoints(8);
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
        return SUM;
    }

    //--------------------------------------------------------------------------------------
    //       Class:  Coupling
    //      Method:  IntegrateOnOctreeGrid
    //--------------------------------------------------------------------------------------
    void Coupling::IntegrateOnOctreeGrid( bool vtkOutput) {

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
                kr->SetName("Re($K_0$)");
                kr->SetNumberOfTuples( oct->GetNumberOfLeaves() );
            vtkDoubleArray* ki = vtkDoubleArray::New();
                ki->SetNumberOfComponents(1);
                ki->SetName("Im($K_0$)");
                ki->SetNumberOfTuples( oct->GetNumberOfLeaves() );
            vtkDoubleArray* km = vtkDoubleArray::New();
                km->SetNumberOfComponents(1);
                km->SetName("mod($K_0$)");
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
                std::string fname = std::string("octree-couple") + std::string(".vto");
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
    }

    //--------------------------------------------------------------------------------------
    //       Class:  Coupling
    //      Method:  f
    //--------------------------------------------------------------------------------------
    Complex Coupling::f( const Vector3r& r, const Real& volume, const Vector3cr& Ht, const Vector3cr& Hr ) {
        return volume*Ht.dot(Hr);
        //return Ht.dot(Hr);
    }

    //--------------------------------------------------------------------------------------
    //       Class:  Coupling
    //      Method:  EvaluateKids
    //--------------------------------------------------------------------------------------
    void Coupling::EvaluateKids( const Vector3r& size, const int& level, const Vector3r& cpos,
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
        if ( std::abs(ksum-parentVal) > tol || level < minLevel && level < maxLevel ) {
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
    //       Class:  Coupling
    //      Method:  EvaluateKids2 -- same as Evaluate Kids, but include VTK octree generation
    //--------------------------------------------------------------------------------------
    void Coupling::EvaluateKids2( const Vector3r& size, const int& level, const Vector3r& cpos,
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
        if ( std::abs(ksum-parentVal) > tol || level < minLevel && level < maxLevel ) {
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
    //       Class:  Coupling
    //      Method:  GetPosition
    //--------------------------------------------------------------------------------------
    void Coupling::GetPosition( vtkHyperOctreeCursor* Cursor, Real* p ) {
        Real ratio=1.0/(1<<(Cursor->GetCurrentLevel()));
        //step  = ((Size).array() / std::pow(2.,Cursor->GetCurrentLevel()));
        p[0]=(Cursor->GetIndex(0)+.5)*ratio*this->Size[0]+this->Origin[0] ;//+ .5*step[0];
        p[1]=(Cursor->GetIndex(1)+.5)*ratio*this->Size[1]+this->Origin[1] ;//+ .5*step[1];
        p[2]=(Cursor->GetIndex(2)+.5)*ratio*this->Size[2]+this->Origin[2] ;//+ .5*step[2];
    }

    #endif

} // ----  end of namespace Lemma  ----

/* vim: set tabstop=4 expandtab */
/* vim: set filetype=cpp */

