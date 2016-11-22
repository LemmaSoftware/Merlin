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
    auto Tx1 = CircularLoop(21, 15, 100, 100);
    auto Tx2 = CircularLoop(21, 15, 100, 124.8);
    //auto Tx1 = CircularLoop(60, 15, 0, 0); // was 60

    auto Kern = KernelV0::NewSP();
        Kern->PushCoil( "Coil 1", Tx1 );
        Kern->PushCoil( "Coil 2", Tx2 );
        Kern->SetLayeredEarthEM( earth );

        Kern->SetIntegrationSize( (Vector3r() << 200,200,200).finished() );
        Kern->SetIntegrationOrigin( (Vector3r() << 0,0,0).finished() );
        Kern->SetTolerance( 1e-10 );

        Kern->SetPulseDuration(0.020);
        VectorXr I(36);
        I << 397.4208916184016, 352.364477036168, 313.0112765842783, 278.37896394065376, 247.81424224324982,
             220.77925043190442, 196.76493264105017, 175.31662279234038, 156.0044839325404, 138.73983004230124,
             123.42064612625474, 109.82713394836259, 97.76534468972267, 87.06061858367781, 77.56000002944572, 69.1280780096311,
             61.64250263640252, 54.99473044877554, 49.091182970515476, 43.84634004556388, 39.184136917167976, 35.03619319797924,
             31.347205894128976, 28.06346770557137, 25.139117042424758, 22.53420773366429, 20.214205433283347,
             18.144318026099942, 16.299965972298878, 14.652633628829891, 13.184271405688083, 11.870540177313893,
             10.697057141915716, 9.64778948429609, 8.709338689612677, 7.871268012862094;
        //Kern->SetPulseCurrent( VectorXr::LinSpaced( 1, 10, 200 )  ); // nbins, low, high
        Kern->SetPulseCurrent( I ); // nbins, low, high

        //Kern->SetDepthLayerInterfaces( VectorXr::LinSpaced( 30, 3, 45.5 ) ); // nlay, low, high
        //10**np.linspace(np.log10(10),np.log10(19),10)
        VectorXr interfaces = VectorXr::LinSpaced(21, std::log10(2), std::log10(50)); // 30 log spaced
        for (int i=0; i<interfaces.size(); ++i) {
            interfaces(i) = std::pow(10, interfaces(i));
        }
        Kern->SetDepthLayerInterfaces( interfaces ); // nlay, low, high

    // We could, I suppose, take the earth model in here? For non-linear that
    // may be more natural to work with?
    std::vector<std::string> tx = {std::string("Coil 1"), std::string("Coil 2") };
    std::vector<std::string> rx = {std::string("Coil 1")};
    Kern->CalculateK0( tx, rx, true );

    ofstream out = ofstream("k.yaml");
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
