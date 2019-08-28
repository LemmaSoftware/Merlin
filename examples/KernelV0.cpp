/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      11/11/2016 02:44:37 PM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2016, University of Utah
 * @copyright Copyright (c) 2016, Lemma Software, LLC
 */

#include <Merlin>
using namespace Lemma;

std::shared_ptr<PolygonalWireAntenna> CircularLoop ( int nd, Real radius, Real Offsetx, Real Offsety ) ;

int main(int argc, char** argv) {

    if (argc < 6) {             // 1          2           3                    4                5
        std::cout << "./KernelV0 TxCoil.yaml RxCoil.yaml  EMEarthModel.yaml  AkvoDataSet.yaml Output.yaml \n";
        exit(EXIT_SUCCESS);
    }

	auto earth = LayeredEarthEM::NewSP();
		earth->SetNumberOfLayers(3);
		earth->SetLayerConductivity( (VectorXcr(3) << Complex(0.,0), Complex(1./50.,0), Complex(1./100.)).finished() );
		earth->SetLayerThickness( (VectorXr(1) << 10).finished() );
        // Set mag field info
        // From NOAA, Laramie WY, June 9 2016, aligned with mag. north
        earth->SetMagneticFieldIncDecMag( 67, 9, 52750, NANOTESLA );

        auto Tx1 = PolygonalWireAntenna::DeSerialize( YAML::LoadFile(argv[1]) );
        auto Tx2 = PolygonalWireAntenna::DeSerialize( YAML::LoadFile(argv[2]) );

        //////////////////////////////////////               _
        Tx1->SetCurrent(1.);                //               |
        Tx1->SetNumberOfTurns(1);           //    Set these from Akvo input!
        Tx1->SetNumberOfFrequencies(1);     //               |
        Tx1->SetFrequency(0,2246);          //               |
        Tx2->SetCurrent(1.);                //               |
        Tx2->SetNumberOfTurns(1);           //           \   |   /
        Tx2->SetNumberOfFrequencies(1);     //            \  |  /
        Tx2->SetFrequency(0,2246);          //             \ | /
        //////////////////////////////////////               _

    auto Kern = KernelV0::NewSP();
        Kern->PushCoil( "Coil 1", Tx1 );
        Kern->PushCoil( "Coil 2", Tx2 );
        Kern->SetLayeredEarthEM( earth );

        Kern->SetIntegrationSize( (Vector3r() << 200,200,200).finished() );
        Kern->SetIntegrationOrigin( (Vector3r() << 0,0,0).finished() );
        Kern->SetTolerance( 1e-12 ); // 1e-12

        auto AkvoDataNode = YAML::LoadFile(argv[4]);
        Kern->AlignWithAkvoDataset( AkvoDataNode );

        // These should to into AlignWithAkvoDataSet...
        Kern->SetPulseDuration( AkvoDataNode["pulseLength"][0].as<Real>() );
        Kern->SetPulseCurrent( AkvoDataNode["Pulses"]["Pulse 1"]["current"].as<VectorXr>() ); // nbins, low, high


        //VectorXr interfaces = VectorXr::LinSpaced( 41, .5, 45.5 ); // nlay, low, high
        VectorXr interfaces = VectorXr::LinSpaced( 51, .5, 45.5 ); // nlay, low, high
        Real thick = .5;
        for (int ilay=1; ilay<interfaces.size(); ++ilay) {
            interfaces(ilay) = interfaces(ilay-1) + thick;
            thick *= 1.05;
        }
        Kern->SetDepthLayerInterfaces( interfaces ); // nlay, low, high

    // We could, I suppose, take the earth model in here? For non-linear that
    // may be more natural to work with?
    //std::vector<std::string> tx = {std::string("Coil 1"), std::string("Coil 2") };
    std::vector<std::string> tx = {std::string("Coil 1")};  //, std::string("Coil 2") };
    std::vector<std::string> rx = {std::string("Coil 2")};
    Kern->CalculateK0( tx, rx, true );

/*
    std::ofstream dout = std::ofstream(std::string("k-Tx2coil-Rx1coil-offset-")+ std::string(argv[1])+ std::string(".dat"));
    //std::ofstream dout = std::ofstream(std::string("k-coincident.dat"));
        dout << interfaces.transpose() << std::endl;
        dout << I.transpose() << std::endl;
        dout << "#real\n";
        dout << Kern->GetKernel().real() << std::endl;
        dout << "#imag\n";
        dout << Kern->GetKernel().imag() << std::endl;
        dout.close();
*/

    //std::ofstream out = std::ofstream(std::string("k-Tx2coil-Rx1coil-offset-")+std::string(argv[1])+std::string(".yaml"));
    std::ofstream out = std::ofstream(std::string(argv[5]));
    //std::ofstream out = std::ofstream(std::string("k-coincident.yaml"));
    out << *Kern;
    out.close();
}

std::shared_ptr<Lemma::PolygonalWireAntenna> CircularLoop ( int nd, Real Radius, Real Offsetx, Real Offsety ) {

    auto Tx1 = Lemma::PolygonalWireAntenna::NewSP();
         Tx1->SetNumberOfPoints(nd);

    VectorXr range = VectorXr::LinSpaced(nd, 0, 2*PI);
    int ii;
    for (ii=0; ii<nd; ++ii) {
        Tx1->SetPoint(ii, Vector3r(Offsetx+Radius*std::cos(range(ii)), Offsety+Radius*std::sin(range(ii)),  -1e-3));
    }
    //Tx1->SetPoint(ii, Vector3r(Offsetx+Radius*1, Offsety,  -1e-3));

    Tx1->SetCurrent(1.);
    Tx1->SetNumberOfTurns(1);
    Tx1->SetNumberOfFrequencies(1);
    Tx1->SetFrequency(0,2246);

    return Tx1;
}
