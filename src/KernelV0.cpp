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
    void KernelV0::CalculateK0 (const std::vector< std::string>& Tx, const std::vector<std::string >& Rx,
            bool vtkOutput ) {

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
                EMEarths[tx]->SetHankelTransformMethod(ANDERSON801);
                EMEarths[tx]->SetTxRxMode(TX);
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
            }
        }

        std::cout << "Calculating K0 kernel\n";
        MatrixXcr Kern = MatrixXcr::Zero( Interfaces.size() - 1, PulseI.size() );
        for (int ilay=0; ilay< Interfaces.size()-1; ++ilay) {
            for (int iq=0; iq< PulseI.size()-1; ++iq) {
                std::cout << "Layer " << ilay << " q " << iq << std::endl;
                Size(2) = Interfaces(ilay+1) - Interfaces(ilay);
                Origin(2) = Interfaces(ilay);
                Ip = PulseI(iq);
                Kern(ilay, iq) = IntegrateOnOctreeGrid( ilay, iq, vtkOutput );
            }
        }
        std::cout << "\rFinished KERNEL\n";
        std::cout << "real\n";
        std::cout << Kern.real() << std::endl;
        std::cout << "imag\n";
        std::cout << Kern.imag() << std::endl;
        //IntegrateOnOctreeGrid( vtkOutput );
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  IntegrateOnOctreeGrid
    //--------------------------------------------------------------------------------------
    Complex KernelV0::IntegrateOnOctreeGrid( const int& ilay, const int& iq, bool vtkOutput) {

        Vector3r cpos = (Size-Origin).array() / 2.;

        SUM = 0;
        VOLSUM = 0;
        nleaves = 0;
        if (!vtkOutput) {
            EvaluateKids( Size, 0, cpos, 1e6 );
        } else {
        #ifdef LEMMAUSEVTK
            vtkHyperOctree* oct = vtkHyperOctree::New();
                oct->SetDimension(3);
                oct->SetOrigin( Origin(0), Origin(1), Origin(2) );
                oct->SetSize( Size(0), Size(1), Size(2) );
            vtkHyperOctreeCursor* curse = oct->NewCellCursor();
                curse->ToRoot();
            EvaluateKids2( Size, 0, cpos, 1e6, oct, curse );

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

            for (auto leaf : LeafDict) {
                kr->InsertTuple1( leaf.first, std::real(leaf.second) );
                ki->InsertTuple1( leaf.first, std::imag(leaf.second) );
                km->InsertTuple1( leaf.first, std::abs(leaf.second) );
                kid->InsertTuple1( leaf.first, leaf.first );
            }

            for (auto leaf : LeafDictIdx) {
                kerr->InsertTuple1( leaf.first, leaf.second );
            }

            oct->GetLeafData()->AddArray(kr);
            oct->GetLeafData()->AddArray(ki);
            oct->GetLeafData()->AddArray(km);
            oct->GetLeafData()->AddArray(kid);
            oct->GetLeafData()->AddArray(kerr);

            auto write = vtkXMLHyperOctreeWriter::New();
                //write.SetDataModeToAscii()
                write->SetInputData(oct);
                std::string fname = std::string("octree-") + to_string(ilay)
                                  + std::string("-") + to_string(iq) + std::string(".vto");
                write->SetFileName(fname.c_str());
                write->Write();
                write->Delete();

            //kerr->Delete();
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
        std::cout << "nleaves\t" << nleaves << std::endl;
        std::cout << "KSUM\t" << SUM << std::endl;

        return SUM;
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  f
    //--------------------------------------------------------------------------------------
    Complex KernelV0::f( const Vector3r& r, const Real& volume, const Vector3cr& Ht, const Vector3cr& Hr ) {
        //return Complex(volume*Ht.dot(Hr));
        return ComputeV0Cell(MU0*Ht, MU0*Hr, volume, 1.0);
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  ComputeV0Cell
    //--------------------------------------------------------------------------------------
    Complex KernelV0::ComputeV0Cell(const Vector3cr& Bt,
                const Vector3cr& Br, const Real& vol, const Real& phi) {

        // Compute the elliptic fields
        Vector3r B0hat = SigmaModel->GetMagneticFieldUnitVector();
        Vector3r B0 = SigmaModel->GetMagneticField();

        // Elliptic representation
        EllipticB EBT = EllipticFieldRep(Bt, B0hat);
        EllipticB EBR = EllipticFieldRep(Br, B0hat);

        // Compute Mn0
        Vector3r Mn0 = ComputeMn0(phi, B0);
        Real Mn0Abs = Mn0.norm();

        // Compute the tipping angle
        Real sintheta = std::sin(0.5*GAMMA*Ip*Taup*std::abs(EBT.alpha-EBT.beta));

        // Compute phase delay, TODO add transmiiter phase and delay time phase!
        Real phase = EBR.zeta+EBT.zeta;

        return ComputeV0Cell(EBT, EBR, sintheta, phase, Mn0Abs, vol);
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  ComputeV0Cell
    //--------------------------------------------------------------------------------------
    Complex KernelV0::ComputeV0Cell(const EllipticB& EBT, const EllipticB& EBR,
                const Real& sintheta, const Real& phase, const Real& Mn0Abs,
                const Real& vol) {

        Vector3r B0hat = {1,0,0};

        // earth response of receiver adjoint field
        Complex ejztr = std::exp(Complex(0, EBR.zeta + EBT.zeta));

        Complex PhaseTerm = EBR.bhat.dot(EBT.bhat) +
               (B0hat.dot(EBR.bhat.cross(EBT.bhat) ));
        return -vol*Complex(0,Larmor)*Mn0Abs*(EBR.alpha+EBR.beta)*ejztr*sintheta*PhaseTerm;
    }

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
        EllipticB ElipB = EllipticB();
        Vector3cr Bperp = B.array() - B0hat.dot(B)*B0hat.array();
        Real BperpNorm = Bperp.norm();
        Complex Bp2 = Bperp.transpose() * Bperp;
        VectorXcr iB0 = Complex(0,1)*B0hat.cast<Complex>().array();
        ElipB.eizt = std::sqrt(Bp2 / std::abs(Bp2));
        ElipB.alpha = INVSQRT2*std::sqrt(BperpNorm*BperpNorm + std::abs(Bp2));
        ElipB.beta = sgn(std::real(iB0.dot(Bperp.cross(Bperp.conjugate())))) *
                (INVSQRT2)*std::sqrt(std::abs(BperpNorm*BperpNorm-std::abs(Bp2)));
        ElipB.bhat = ((Real)1./ElipB.alpha)*(((Real)1./ElipB.eizt)*Bperp.array()).real().array();
        ElipB.bhatp = B0hat.cross(ElipB.bhat);
        ElipB.zeta = std::real(std::log(ElipB.eizt)/Complex(0,1));
        return ElipB;
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  EvaluateKids
    //--------------------------------------------------------------------------------------
    bool KernelV0::EvaluateKids( const Vector3r& size, const int& level, const Vector3r& cpos,
        const Complex& parentVal ) {

        std::cout << "\r" << (int)(1e2*VOLSUM/(Size[0]*Size[1]*Size[2])) << "\t" << nleaves;
        std::cout.flush();

        // Next level step, interested in one level below
        // bitshift requires one extra, faster than, and equivalent to std::pow(2, level+1)
        Vector3r step  = size.array() / (Real)(1 << (level+1) );
        Vector3r step2 = size.array() / (Real)(1 << (level+2) );

        Real vol = (step2(0)*step2(1)*step2(2));     // volume of each child

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
        //Eigen::Matrix< Complex, 8, 3 > Bt;
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
        if ( std::abs(ksum - parentVal) > tol || level < minLevel && level < maxLevel ) {
            for (int ichild=0; ichild<8; ++ichild) {
                Vector3r cp = pos; // Eigen complains about combining these
                cp += posadd.row(ichild);
                bool isleaf = EvaluateKids( size, level+1, cp, kvals(ichild) );
                if (isleaf) {      // Include result in final integral
                    SUM += ksum;
                    VOLSUM += 8.*vol;
                    nleaves += 1;
                }
            }
            return false;  // not leaf
        }
        // Save here instead?
        return true;       // leaf
    }

    #ifdef LEMMAUSEVTK
    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  EvaluateKids2 -- same as Evaluate Kids, but include VTK octree generation
    //--------------------------------------------------------------------------------------
    bool KernelV0::EvaluateKids2( const Vector3r& size, const int& level, const Vector3r& cpos,
        const Complex& parentVal, vtkHyperOctree* oct, vtkHyperOctreeCursor* curse) {

        std::cout << "\r" << (int)(1e2*VOLSUM/(Size[0]*Size[1]*Size[2])) << "\t" << nleaves;
        std::cout.flush();

        // Next level step, interested in one level below
        // bitshift requires one extra, faster than, and equivalent to std::pow(2, level+1)
        Vector3r step  = size.array() / (Real)(1 << (level+1) );
        Vector3r step2 = size.array() / (Real)(1 << (level+2) );

        Real vol = (step2(0)*step2(1)*step2(2));     // volume of each child

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
            Vector3r cp = pos; // Eigen complains about combining these
            cp += posadd.row(ichild);
            kvals(ichild) = f(cp, vol, Ht.col(ichild), Hr.col(ichild));
        }

        Complex ksum = kvals.sum();     // Kernel sum
        // Evaluate whether or not furthur splitting is needed
        Real err = std::abs(ksum - parentVal);
        if ( std::abs(ksum - parentVal) > tol || level < minLevel && level < maxLevel ) {
            oct->SubdivideLeaf(curse);
            for (int ichild=0; ichild<8; ++ichild) {
                curse->ToChild(ichild);
                Vector3r cp = pos; // Eigen complains about combining these
                cp += posadd.row(ichild);
                /* Test for position via alternative means
                Real p[3];
                GetPosition(curse, p);
                std::cout << cp[0] << "\t" << p[0] << "\t" << cp[1] << "\t" << p[1]
                          << "\t" << cp[2] << "\t" << p[2] << "\t" <<  vol<< std::endl;
                */
                bool isleaf = EvaluateKids2( size, level+1, cp, kvals(ichild), oct, curse );
                if (isleaf) {  // Include result in final integral
                    LeafDict[curse->GetLeafId()] = kvals(ichild);       // VTK
                    LeafDictIdx[curse->GetLeafId()] = nleaves;          // VTK
                    SUM += ksum;
                    VOLSUM += 8*vol;
                    nleaves += 1;
                }
                curse->ToParent();
            }
            return false;  // not leaf
        }
        return true;       // leaf
    }

    //--------------------------------------------------------------------------------------
    //       Class:  KernelV0
    //      Method:  GetPosition
    //--------------------------------------------------------------------------------------
    void KernelV0::GetPosition( vtkHyperOctreeCursor* Cursor, Real* p ) {
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

