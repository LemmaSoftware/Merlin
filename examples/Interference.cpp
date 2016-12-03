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

static constexpr Real GAMMA = 2.67518e8;                  // MKS units

std::shared_ptr<PolygonalWireAntenna> CircularLoop ( int nd, Real radius, Real Offsetx, Real Offsety, Real wL ) ;
void MoveLoop( std::shared_ptr<PolygonalWireAntenna> Loop, int nd, Real Radius, Real Offsetx, Real Offsety, Real wL );

int main(int argc, char** argv) {
    /*
    if ( argc < 2 ) {
        std::cout << "Calculates the coupling between two sNMR loops at the Larmor frequency. Usage\n"
                  << "\t./Coupling   EarthModel.yaml" << std::endl;
        exit(0);
    }
    //Real offset = atof(argv[1]);
    auto earth = LayeredEarthEM::DeSerialize( YAML::LoadFile(argv[1]) );
    */
    // RedButtes model, also how you can generate your own files
	auto earth = LayeredEarthEM::NewSP();
		earth->SetNumberOfLayers(3);
		earth->SetLayerConductivity( (VectorXcr(3) << Complex(0.,0), Complex(1./50.,0), Complex(1./100.)).finished() );
		earth->SetLayerThickness( (VectorXr(1) << 10).finished() );
        // Set mag field info
        // From NOAA, Laramie WY, June 9 2016, aligned with mag. north
        earth->SetMagneticFieldIncDecMag( 67, 0, 52750, NANOTESLA );
//     auto sig = std::ofstream("SigmaModel.yaml");
//         sig << *earth << std::endl;
//         sig.close();

    Real Larmor = earth->GetMagneticFieldMagnitude()*GAMMA/(2*PI);

    // Transmitter loops
    auto Tx1 = CircularLoop(21, 15, 50, 50, Larmor);
    auto Tx2 = CircularLoop(21, 15, 50, 50, Larmor); // initially coincident

    auto Kern = LoopInteractions<INTERFERENCE>::NewSP();
        Kern->PushCoil( "Coil 1", Tx1 );
        Kern->PushCoil( "Coil 2", Tx2 );
        Kern->SetLayeredEarthEM( earth );

        Kern->SetIntegrationSize( (Vector3r()   << 50,200,50).finished() );
        Kern->SetIntegrationOrigin( (Vector3r() << 0,0,5.0).finished() );
        Kern->SetTolerance( 1e-3 ); // 1e-12

    std::vector<std::string> tx = {std::string("Coil 1")};//,std::string("Coil 2")};
    std::vector<std::string> rx = {std::string("Coil 2")};
    VectorXr Offsets = VectorXr::LinSpaced(61, 0.00, 60.0); // nbins, low, high

    auto outfile = std::ofstream("interference.dat");
    for (int ii=0; ii< Offsets.size(); ++ii) {
        MoveLoop(Tx2, 21, 15, 50, 50 + Offsets(ii), Larmor);
        #ifdef LEMMAUSEVTK
        Complex coupling = Kern->Calculate( tx, rx, true );
        #else
        Complex coupling = Kern->Calculate( tx, rx, false );
        #endif
        std::cout << "coupling " << coupling << std::endl;
        outfile << Offsets(ii) << "\t" <<  std::real(coupling) << "\t" << std::imag(coupling) << std::endl;
    }
    outfile.close();
}

std::shared_ptr<Lemma::PolygonalWireAntenna> CircularLoop ( int nd, Real Radius, Real Offsetx, Real Offsety, Real wL ) {

    auto Tx1 = Lemma::PolygonalWireAntenna::NewSP();
         Tx1->SetNumberOfPoints(nd);

    VectorXr range = VectorXr::LinSpaced(nd, 0, 2*PI);
    int ii;
    for (ii=0; ii<nd; ++ii) {
        Tx1->SetPoint(ii, Vector3r(Offsetx+Radius*std::cos(range(ii)), Offsety+Radius*std::sin(range(ii)),  -1e-3));
    }

    Tx1->SetCurrent(1.);
    Tx1->SetNumberOfTurns(1);
    Tx1->SetNumberOfFrequencies(1);
    Tx1->SetFrequency(0,wL);

    return Tx1;
}

void MoveLoop( std::shared_ptr<Lemma::PolygonalWireAntenna> Tx1, int nd, Real Radius, Real Offsetx, Real Offsety, Real wL ) {

    Tx1->SetNumberOfPoints(nd);

    VectorXr range = VectorXr::LinSpaced(nd, 0, 2*PI);
    int ii;
    for (ii=0; ii<nd; ++ii) {
        Tx1->SetPoint(ii, Vector3r(Offsetx+Radius*std::cos(range(ii)), Offsety+Radius*std::sin(range(ii)),  -1e-3));
    }

    Tx1->SetCurrent(1.);
    Tx1->SetNumberOfTurns(1);
    Tx1->SetNumberOfFrequencies(1);
    Tx1->SetFrequency(0,wL);

}
