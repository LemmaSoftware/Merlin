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
 * @email     Trevor.Irons@lemmasoftware.org
 * @copyright Copyright (c) 2016, University of Utah
 * @copyright Copyright (c) 2019, Trevor P. Irons
 * @copyright Copyright (c) 2016, Lemma Software, LLC
 * @copyright Copyright (c) 2008, Colorado School of Mines
 */

#include "MerlinConfig.h"
#include "KernelV0.h"
#include "FieldPoints.h"

namespace Lemma {

    // ====================  FRIEND METHODS  =====================

    std::ostream &operator << (std::ostream &stream, const KernelV0 &ob) {
        stream << ob.Serialize()  << "\n"; // End of doc ---
        return stream;
    }

    // ====================  LIFECYCLE     =======================

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  KernelV0
    // Description:  constructor (locked)
    //--------------------------------------------------------------------------------------
    KernelV0::KernelV0 (const ctor_key& key) : MerlinObject( key ) {

    }  // -----  end of method KernelV0::KernelV0  (constructor)  -----

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  KernelV0
    // Description:  DeSerializing constructor (locked)
    //--------------------------------------------------------------------------------------
    KernelV0::KernelV0 (const YAML::Node& node, const ctor_key& key) : MerlinObject(node, key) {

        //node["PulseType"] = "FID";
        Larmor = node["Larmor"].as<Real>();
        Temperature = node["Temperature"].as<Real>();
        tol = node["tol"].as<Real>();
        minLevel = node["minLevel"].as<int>();
        maxLevel = node["maxLevel"].as<int>();
        Interfaces = node["Interfaces"].as<VectorXr>();
        Size = node["IntegrationSize"].as<Vector3r>();
        Origin = node["IntegrationOrigin"].as<Vector3r>();

        if (node["AlignWithAkvoData"]) {
            // Match pulse info with dataset
            AlignWithAkvoDataset( YAML::LoadFile( node["AlignWithAkvoData"].as<std::string>()));
        } else {
            // Read Pulse info direct from Kernel file
            PulseI = node["PulseI"].as<VectorXr>();
            Taup = node["Taup"].as<Real>();
        }

        if (node["SigmaModel"]) {
            if (node["SigmaModel"].Tag() == "LayeredEarthEM") {
                SigmaModel = LayeredEarthEM::DeSerialize(node["SigmaModel"]);
            } else {
                SigmaModel = LayeredEarthEM::DeSerialize( YAML::LoadFile( node["SigmaModel"].as<std::string>() ));
            }
        }

        if (node["Coils"]) {
            for ( auto coil : node["Coils"] ) {
                if ( coil.second.Tag() == "PolygonalWireAntenna" ) {
                    TxRx[ coil.first.as<std::string>() ] = PolygonalWireAntenna::DeSerialize( coil.second );
                } else {
                    TxRx[ coil.first.as<std::string>() ] =
                        PolygonalWireAntenna::DeSerialize( YAML::LoadFile(coil.second.as<std::string>()) );
                }
            }
        }

        if (node["K0"]) {
            Kern = MatrixXcr::Zero( Interfaces.size()-1, PulseI.size()  ).array() + 1.;
            for ( int ilay=0; ilay<Interfaces.size()-1; ++ilay ) {
                Kern.row(ilay) = node["K0"]["layer-" + to_string(ilay) ].as<VectorXcr>();
            }
        }
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

        YAML::Node node = MerlinObject::Serialize();
        node.SetTag( GetName() );

        // Coils Transmitters & Receivers
        if (!TxRx.empty()) {
            for ( auto txm : TxRx) {
                node["Coils"][txm.first] = txm.second->Serialize();
            }
        }

        // LayeredEarthEM
        if (SigmaModel != nullptr) {
            node["SigmaModel"] = SigmaModel->Serialize();
        }

        node["PulseType"] = "FID";
        node["Larmor"] = Larmor;
        node["Temperature"] = Temperature;
        node["tol"] = tol;
        node["minLevel"] = minLevel;
        node["maxLevel"] = maxLevel;
        node["Taup"] = Taup;
        node["PulseI"] = PulseI;
        node["Interfaces"] = Interfaces;
        node["IntegrationSize"] = Size;
        node["IntegrationOrigin"] = Origin;

        // TODO, use better matrix encapulation
        if (Kern.array().abs().any() > 1e-16) {
            for ( int ilay=0; ilay<Interfaces.size()-1; ++ilay ) {
                node["K0"]["layer-" + to_string(ilay) ] = static_cast<VectorXcr>(Kern.row(ilay));
            }
        }
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
    //      Method:  AlignWithAkvoDataset
    //--------------------------------------------------------------------------------------
    void KernelV0::AlignWithAkvoDataset( const YAML::Node& node ) {
        if (node["processed"].as<std::string>().substr(0,4) == "Akvo") {
            std::cout << "Akvo file read\n";
            std::cout << node["processed"] << std::endl;
        }
        if (node["pulseType"].as<std::string>() == "FID") {
            std::cout << "FID pulse detected" << std::endl;
            PulseI  = node["Pulses"]["Pulse 1"]["current"].as<VectorXr>();
            this->SetPulseDuration( node["pulseLength"][0].as<double>() );
        } else {
            std::cerr << "Pulse Type " << node["PulseType"] << "is not supported\n";
        }
        std::cout << "Finished with Akvo file read" << std::endl;
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  DeSerialize
    //--------------------------------------------------------------------------------------
    void KernelV0::CalculateK0 (const std::vector< std::string>& Tx,
                                const std::vector<std::string >& Rx, bool vtkOutput ) {
        // Set up
        Larmor = SigmaModel->GetMagneticFieldMagnitude()*GAMMA; // in rad  2246.*2.*PI;

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
                //EMEarths[tx]->SetHankelTransformMethod(FHTKEY201);
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
                    //EMEarths[rx]->SetHankelTransformMethod(FHTKEY201);
                    EMEarths[rx]->SetHankelTransformMethod(ANDERSON801);
                    EMEarths[rx]->SetTxRxMode(RX);
                    TxRx[rx]->SetCurrent(1.);
            }
        }

        std::cout << "Calculating K0 kernel\n";
        Kern = MatrixXcr::Zero( Interfaces.size()-1, PulseI.size() );
        for (ilay=0; ilay<Interfaces.size()-1; ++ilay) {
            std::cout << "Layer " << ilay << "\tfrom " << Interfaces(ilay) <<" to "
                      << Interfaces(ilay+1) << std::endl;
            Size(2) = Interfaces(ilay+1) - Interfaces(ilay);
            Origin(2) = Interfaces(ilay);
            IntegrateOnOctreeGrid( vtkOutput );
        }
        std::cout << "\nFinished KERNEL\n";
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  IntegrateOnOctreeGrid
    //--------------------------------------------------------------------------------------
    void KernelV0::IntegrateOnOctreeGrid( bool vtkOutput) {

        Vector3r cpos = Origin + Size/2.;

        VOLSUM = 0;
        nleaves = 0;
        if (!vtkOutput) {
            EvaluateKids( Size, 0, cpos, VectorXcr::Ones(PulseI.size()) );
        } else {
        #ifdef LEMMAUSEVTK
            vtkHyperTreeGrid* oct = vtkHyperTreeGrid::New();
                oct->SetGridSize( 1, 1, 1 ); // Important!
                oct->SetDimension(3);
                vtkNew<vtkDoubleArray> xcoords;
                    xcoords->SetNumberOfComponents(1);
                    xcoords->SetNumberOfTuples(2);
                    xcoords->SetTuple1( 0, Origin(0) );
                    xcoords->SetTuple1( 1, Origin(0) + Size(0) );
                    xcoords->SetName("northing (m)");
                    oct->SetXCoordinates(xcoords);

                vtkNew<vtkDoubleArray> ycoords;
                    ycoords->SetNumberOfComponents(1);
                    ycoords->SetNumberOfTuples(2);
                    ycoords->SetTuple1( 0, Origin(1) );
                    ycoords->SetTuple1( 1, Origin(1) + Size(1) );
                    ycoords->SetName("easting (m)");
                    oct->SetYCoordinates(ycoords);

                vtkNew<vtkDoubleArray> zcoords;
                    zcoords->SetNumberOfComponents(1);
                    zcoords->SetNumberOfTuples(2);
                    zcoords->SetTuple1( 0, Origin(2) );
                    zcoords->SetTuple1( 1, Origin(2) + Size(2) );
                    zcoords->SetName("depth (m)");
                    oct->SetZCoordinates(zcoords);

            //vtkHyperTreeGridLevelEntry* curse2 =  vtkHyperTreeGridLevelEntry::New(); // VTK 9
            // I belive the index in NewCursor maybe points to which cell in the Rectilinear Grid?
            vtkHyperTreeCursor* curse = oct->NewCursor(0, true); // true creates the cursor
                curse->ToRoot();
            EvaluateKids2( Size, 0, cpos, VectorXcr::Ones(PulseI.size()), oct, curse );

            for (int iq=0; iq<PulseI.size(); ++iq) {

            // Fill in leaf data
            vtkDoubleArray* kr = vtkDoubleArray::New();
                kr->SetNumberOfComponents(1);
                kr->SetName("Re($\\mathcal{K}_0$)");
                kr->SetNumberOfTuples( oct->GetNumberOfLeaves() );

            vtkDoubleArray* ki = vtkDoubleArray::New();
                ki->SetNumberOfComponents(1);
                ki->SetName("Im($\\mathcal{K}_0$)");
                ki->SetNumberOfTuples( oct->GetNumberOfLeaves() );

            vtkDoubleArray* km = vtkDoubleArray::New();
                km->SetNumberOfComponents(1);
                km->SetName("mod($\\mathcal{K}_0$)");
                km->SetNumberOfTuples( oct->GetNumberOfLeaves() );

            vtkIntArray* kid = vtkIntArray::New();
                kid->SetNumberOfComponents(1);
                kid->SetName("ID");
                kid->SetNumberOfTuples( oct->GetNumberOfLeaves() );

            vtkIntArray* kerr = vtkIntArray::New();
                kerr->SetNumberOfComponents(1);
                kerr->SetName("err");
                kerr->SetNumberOfTuples( oct->GetNumberOfLeaves() );

            // Ht field
            vtkDoubleArray* htr = vtkDoubleArray::New();
                htr->SetNumberOfComponents(3);
                htr->SetName("Re($\\mathbf{\\mathcal{H}}_T$)");
                htr->SetNumberOfTuples( oct->GetNumberOfLeaves() );

            vtkDoubleArray* hti = vtkDoubleArray::New();
                hti->SetNumberOfComponents(3);
                hti->SetName("Im($\\mathbf{\\mathcal{H}}_T$)");
                hti->SetNumberOfTuples( oct->GetNumberOfLeaves() );

            // Hr field
            vtkDoubleArray* hrr = vtkDoubleArray::New();
                hrr->SetNumberOfComponents(3);
                hrr->SetName("Re($\\mathbf{\\mathcal{H}}_R$)");
                hrr->SetNumberOfTuples( oct->GetNumberOfLeaves() );

            vtkDoubleArray* hri = vtkDoubleArray::New();
                hri->SetNumberOfComponents(3);
                hri->SetName("Im($\\mathbf{\\mathcal{H}}_R$)");
                hri->SetNumberOfTuples( oct->GetNumberOfLeaves() );
            //Real LeafVol(0);

            //kr->Fill(0);
            int icc(0);
            for (auto leaf : LeafDict) {
                kr->InsertTuple1( leaf.first, std::real(leaf.second(iq)) );
                ki->InsertTuple1( leaf.first, std::imag(leaf.second(iq)) );
                km->InsertTuple1( leaf.first, std::abs(leaf.second(iq)) );
                kid->InsertTuple1( leaf.first, leaf.first );
                //LeafVol += std::real(leaf.second);
                ++icc;
            }

            for (auto leaf : LeafHt ) {
                htr->InsertTuple( leaf.first, leaf.second.real().data() );
                hti->InsertTuple( leaf.first, leaf.second.imag().data() );
            }

            for (auto leaf : LeafHr ) {
                hrr->InsertTuple( leaf.first, leaf.second.real().data() );
                hri->InsertTuple( leaf.first, leaf.second.imag().data() );
            }

            for (auto leaf : LeafDictIdx) {
                kerr->InsertTuple1( leaf.first, leaf.second );
            }

            // In VTK 8, vtkHyperTreeGrid does not support CellData.
            // the previous class vtkHyperOctreeGrid used "LeafData",
            // but for the new classes, PointData seems to function as LeafData.
            // this could change in VTK 9
            auto kri = oct->GetPointData()->AddArray(kr);
            auto kii = oct->GetPointData()->AddArray(ki);
            auto kmi = oct->GetPointData()->AddArray(km);
            auto kidi = oct->GetPointData()->AddArray(kid);
            auto keri = oct->GetPointData()->AddArray(kerr);
            auto khtr = oct->GetPointData()->AddArray(htr);
            auto khti = oct->GetPointData()->AddArray(hti);
            auto khrr = oct->GetPointData()->AddArray(hrr);
            auto khri = oct->GetPointData()->AddArray(hri);

            //std::cout << *oct << std::endl;
            auto write = vtkXMLHyperTreeGridWriter::New();
                std::string fname = std::string("octree-") + to_string(ilay)
                                  + std::string("-") + to_string(iq) + std::string(".htg");
                write->SetFileName(fname.c_str());
                write->SetInputData(oct);
                write->SetDataModeToBinary();
                //write->SetDataModeToAscii();
                write->Update();
                write->Write();
                write->Delete();

            oct->GetPointData()->RemoveArray( kri );
            oct->GetPointData()->RemoveArray( kii );
            oct->GetPointData()->RemoveArray( kmi );
            oct->GetPointData()->RemoveArray( kidi );
            oct->GetPointData()->RemoveArray( keri );
            oct->GetPointData()->RemoveArray( khtr );
            oct->GetPointData()->RemoveArray( khti );
            oct->GetPointData()->RemoveArray( khrr );
            oct->GetPointData()->RemoveArray( khri );

            kerr->Delete();
            kid->Delete();
            kr->Delete();
            ki->Delete();
            km->Delete();
            htr->Delete();
            hti->Delete();
            hrr->Delete();
            hri->Delete();
            }

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
    //       Class:  KernelV0
    //      Method:  f
    //--------------------------------------------------------------------------------------
    VectorXcr KernelV0::f( const Vector3r& r, const Real& volume, const Vector3cr& Ht, const Vector3cr& Hr ) {

        // Compute the elliptic fields
        Vector3r B0hat = SigmaModel->GetMagneticFieldUnitVector();
        Vector3r B0 = SigmaModel->GetMagneticField();

        // Elliptic representation
        EllipticB EBT = EllipticFieldRep(MU0*Ht, B0hat);
        EllipticB EBR = EllipticFieldRep(MU0*Hr, B0hat);

        // Compute Mn0
        Vector3r Mn0 = ComputeMn0(1.0, B0);
        Real Mn0Abs = Mn0.norm();
        //std::cout << "Mn0\t" << Mn0.transpose() << std::endl;

        // Compute phase delay
        // TODO add transmiiter current phase and delay induced apparent time phase!
        Complex PhaseTerm = EBR.bhat.dot(EBT.bhat) + Complex(0, (B0hat.dot(EBR.bhat.cross(EBT.bhat))));
        Complex ejztr = std::exp(Complex(0, EBR.zeta + EBT.zeta));

        // Calcuate vector of all responses
        VectorXcr F = VectorXcr::Zero( PulseI.size() );
        for (int iq=0; iq<PulseI.size(); ++iq) {
            // Compute the tipping angle
            Real sintheta = std::sin(0.5*GAMMA*PulseI(iq)*Taup*(EBT.alpha-EBT.beta));
            F(iq) = -volume*Complex(0,Larmor)*Mn0Abs*(EBR.alpha+EBR.beta)*ejztr*sintheta*PhaseTerm;
        }
        return F;
    }

//     //--------------------------------------------------------------------------------------
//     //       Class:  KernelV0
//     //      Method:  ComputeV0Cell
//     //--------------------------------------------------------------------------------------
//     Complex KernelV0::ComputeV0Cell(const EllipticB& EBT, const EllipticB& EBR,
//                 const Real& sintheta, const Real& phase, const Real& Mn0Abs,
//                 const Real& vol) {
//         // earth response of receiver adjoint field
//         Vector3r B0hat = SigmaModel->GetMagneticFieldUnitVector();
//         Complex ejztr = std::exp(Complex(0, EBR.zeta + EBT.zeta));
//         Complex PhaseTerm = EBR.bhat.dot(EBT.bhat) + (B0hat.dot(EBR.bhat.cross(EBT.bhat) ));
//         return -vol*Complex(0,Larmor)*Mn0Abs*(EBR.alpha+EBR.beta)*ejztr*sintheta*PhaseTerm;
//     }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  ComputeV0Cell
    //--------------------------------------------------------------------------------------
    Vector3r KernelV0::ComputeMn0(const Real& Porosity, const Vector3r& B0) {
        Real chi_n = NH2O*((GAMMA*GAMMA*HBAR*HBAR)/(4.*KB*Temperature));
        return chi_n*Porosity*B0;
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  ComputeV0Cell
    //--------------------------------------------------------------------------------------
    EllipticB KernelV0::EllipticFieldRep (const Vector3cr& B, const Vector3r& B0hat) {
        // This all follows Weichman et al., 2000.
        // There are some numerical stability issues that arise when the two terms in the beta
        // formulation are nearly equivalent. The current formulation will result in a null-valued
        // beta, or can underflow. However, this does not entirely recreate the true value of B perp.
        // Error is checked to be below 1%, but reformulating for numeric stability may be welcome
        EllipticB ElipB = EllipticB();
        Vector3cr Bperp = B - B0hat.dot(B)*B0hat;
        Real BperpNorm  = Bperp.norm();
        // These two are equivalent
        //Complex Bp2 = Bperp.transpose() * Bperp;
        Complex Bp2 = Bperp.conjugate().dot(Bperp);
        VectorXcr iB0 = Complex(0,1)*B0hat.cast<Complex>().array();
        ElipB.eizt = std::sqrt(Bp2 / std::abs(Bp2));
        ElipB.alpha = INVSQRT2*std::sqrt(BperpNorm*BperpNorm + std::abs(Bp2));
        //ElipB.beta = std::copysign(1, std::real(iB0.dot( Bperp.cross(Bperp.conjugate())) )) *
        ElipB.beta = sgn( std::real(iB0.dot( Bperp.cross(Bperp.conjugate())) )) *
                     (INVSQRT2*std::sqrt(BperpNorm*BperpNorm - std::abs(Bp2)));
        // Correct underflow in beta calculation
        // could use cerrno instead...
        // http://en.cppreference.com/w/cpp/numeric/math/sqrt
        if (ElipB.beta != ElipB.beta) ElipB.beta = 0;
        ElipB.bhat = ((Real)1./ElipB.alpha)*(((Real)1./ElipB.eizt)*Bperp.array()).real().array();
        ElipB.bhatp = B0hat.cross(ElipB.bhat);
        ElipB.zeta = std::real(std::log(ElipB.eizt)/Complex(0,1));
        /* as an error check decomposed field - computed actual */
//         Vector3cr Bperp2 = ElipB.eizt * (ElipB.alpha * ElipB.bhat
//                        + (Complex(0,1) * ElipB.beta * ElipB.bhatp) );
//         ElipB.err = (Bperp-Bperp2).norm();
//         if (ElipB.err > .01*Bperp.norm() ) { // 1% error
//             std::cout << "Elip error\n";
//             Real Beta2 = sgn( std::real(iB0.dot( Bperp.cross(Bperp.conjugate())) )) *
//                      (INVSQRT2*std::sqrt(BperpNorm*BperpNorm - std::abs(Bp2)));
//             Vector3cr Bperp3 = ElipB.eizt * (ElipB.alpha * ElipB.bhat
//                            + (Complex(0,1) * Beta2 * ElipB.bhatp) );
//             std::cout << "Beta term0\t" << (INVSQRT2*std::sqrt(BperpNorm*BperpNorm - std::abs(Bp2))) << std::endl;
//             std::cout << "Beta term1\t" << BperpNorm*BperpNorm << "\t" << std::abs(Bp2) << std::endl;
//             std::cout << "Beta  \t" << ElipB.beta << std::endl;
//             std::cout << "Beta2 \t" << Beta2 << std::endl;
//             std::cout << "Bperp \t" << Bperp.transpose() << std::endl;
//             std::cout << "Bperp2\t" << Bperp2.transpose() << std::endl;
//             std::cout << "Bperp3\t" << Bperp3.transpose() << std::endl;
//             std::cout << "err   \t" << ElipB.err << std::endl;
//         }
        return ElipB;
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  EvaluateKids
    //--------------------------------------------------------------------------------------
    void KernelV0::EvaluateKids( const Vector3r& size, const int& level, const Vector3r& cpos,
        const VectorXcr& parentVal ) {

        std::cout << "\r" << (int)(1e2*VOLSUM/(Size[0]*Size[1]*Size[2])) << "\t" << nleaves;
        //std::cout.flush();

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

        MatrixXcr kvals(8, PulseI.size());       // individual kernel vals
        for (int ichild=0; ichild<8; ++ichild) {
            Vector3r cp = pos;    // Eigen complains about combining these
            cp += posadd.row(ichild);
            kvals.row(ichild) = f(cp, vol, Ht.col(ichild), Hr.col(ichild));
        }

        VectorXcr ksum = kvals.colwise().sum();     // Kernel sum
        // Evaluate whether or not furthur splitting is needed
        if ( (((ksum - parentVal).array().abs() > tol).any() && level<maxLevel) || level < minLevel ) {
            // Not a leaf dive further in
            for (int ichild=0; ichild<8; ++ichild) {
                Vector3r cp = pos; // Eigen complains about combining these
                cp += posadd.row(ichild);
                EvaluateKids( size, level+1, cp, kvals.row(ichild) );
            }
            return; // not leaf
        }
        // implicit else, is a leaf
        Kern.row(ilay) += ksum;
        VOLSUM += 8.*vol;
        nleaves += 8; // reflects the number of kernel evaluations
        return;     // is leaf
    }

    #ifdef LEMMAUSEVTK
    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  EvaluateKids2 -- same as Evaluate Kids, but include VTK octree generation
    //--------------------------------------------------------------------------------------
    void KernelV0::EvaluateKids2( const Vector3r& size, const int& level, const Vector3r& cpos,
        const VectorXcr& parentVal, vtkHyperTreeGrid* oct, vtkHyperTreeCursor* curse) {

        std::cout << "\r" << (int)(1e2*VOLSUM/(Size[0]*Size[1]*Size[2])) << "\t" << nleaves;
        //std::cout.flush();

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

        MatrixXcr kvals(8, PulseI.size());       // individual kernel vals
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
            kvals.row(ichild) = f(cp, vol, Ht.col(ichild), Hr.col(ichild));
        }

        VectorXcr ksum = kvals.colwise().sum();     // Kernel sum
        // Evaluate whether or not furthur splitting is needed
        if ( (((ksum - parentVal).array().abs() > tol).any() && level<maxLevel) || level < minLevel ) {
            // 0 after curse is vtkIdType?
            oct->SubdivideLeaf(curse, 0);
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
                EvaluateKids2( size, level+1, cp, kvals.row(ichild), oct, curse );
                curse->ToParent();
            }
            return;  // not a leaf
        }
        /* just stuff with sum of the kids and don't subdivide */
        /*
        LeafDict[curse->GetLeafId()] = ksum/(8.*vol);
        LeafDictIdx[curse->GetLeafId()] = nleaves;
        */
        /* Alternatively, subdivide the VTK octree here and stuff the children to make better
         * visuals, but also 8x the storage...
         */

        // 0 after curse is vtkIdType?
        oct->SubdivideLeaf(curse, 0);
        for (int ichild=0; ichild<8; ++ichild) {
            curse->ToChild(ichild);
            LeafDict[curse->GetVertexId()] = ksum/(8.*vol);
            LeafHt[curse->GetVertexId()] = Ht.col(ichild);
            LeafHr[curse->GetVertexId()] = Hr.col(ichild);
            LeafDictIdx[curse->GetVertexId()] = nleaves;
            curse->ToParent();
        }

        Kern.row(ilay) += ksum;
        VOLSUM += 8*vol;
        nleaves += 8; // good reason to say 1 or 8 here...8 sounds better and reflects kernel evaluations
        return;     // is a leaf
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  GetPosition
    //--------------------------------------------------------------------------------------
    void KernelV0::GetPosition( vtkHyperTreeCursor* Cursor, Real* p ) {
        // TODO fix
        /*
        Real ratio=1.0/(1<<(Cursor->GetCurrentLevel()));
        //step  = ((Size).array() / std::pow(2.,Cursor->GetCurrentLevel()));
        p[0]=(Cursor->GetIndex(0)+.5)*ratio*this->Size[0]+this->Origin[0] ;//+ .5*step[0];
        p[1]=(Cursor->GetIndex(1)+.5)*ratio*this->Size[1]+this->Origin[1] ;//+ .5*step[1];
        p[2]=(Cursor->GetIndex(2)+.5)*ratio*this->Size[2]+this->Origin[2] ;//+ .5*step[2];
        */
    }

    #endif

} // ----  end of namespace Lemma  ----

/* vim: set tabstop=4 expandtab */
/* vim: set filetype=cpp */

