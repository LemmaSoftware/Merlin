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

int main() {

	auto earth = LayeredEarthEM::NewSP();
		earth->SetNumberOfLayers(3);
		earth->SetLayerConductivity( (VectorXcr(3) << Complex(0.,0), Complex(1./50.,0), Complex(1./100.)).finished() );
		earth->SetLayerThickness( (VectorXr(1) << 10).finished() );
        // Set mag field info
        // From NOAA, Laramie WY, June 9 2016, aligned with mag. north
        earth->SetMagneticFieldIncDecMag( 67, 0, 52750, NANOTESLA );

    // Transmitter loops
    auto Tx1 = CircularLoop(31, 15, 100, 100);
    auto Tx2 = CircularLoop(31, 15, 100, 120);
    //auto Tx1 = CircularLoop(60, 15, 0, 0); // was 60

    auto Kern = KernelV0::NewSP();
        Kern->PushCoil( "Coil 1", Tx1 );
        Kern->PushCoil( "Coil 2", Tx2 );
        Kern->SetLayeredEarthEM( earth );
        // std::cout << *Kern << std::endl;

        Kern->SetIntegrationSize( (Vector3r() << 200,200,2).finished() );
        Kern->SetIntegrationOrigin( (Vector3r() << 0,0,15).finished() );
        Kern->SetTolerance( 1e-13 );

        Kern->SetPulseDuration(0.020);
        Kern->SetPulseCurrent( VectorXr::LinSpaced( 20, .01, 200 )  ); // nbins, low, high
        Kern->SetDepthLayerInterfaces( VectorXr::LinSpaced( 20, .5, 50 ) );

    // We could, I suppose, take the earth model in here? For non-linear that
    // may be more natural to work with?
    std::vector<std::string> tx = {std::string("Coil 1")};
    std::vector<std::string> rx = {std::string("Coil 1")};
    Kern->CalculateK0( tx, rx, true );

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
