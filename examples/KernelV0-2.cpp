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


int main(int argc, char** argv) {

    if (argc<4) {
        std::cout << "./KernelV0-2 earth.yaml tx.yaml  rx.yaml \n";
        return(EXIT_SUCCESS);
    }
    std::cout << "Using earth model: " << argv[1] << std::endl;
    auto earth = LayeredEarthEM::DeSerialize( YAML::LoadFile(argv[1]) );

    std::cout << "Using transmitter: " << argv[2] << std::endl;
    auto Tx = PolygonalWireAntenna::DeSerialize( YAML::LoadFile(argv[2]) );

    std::cout << "Using receivers: " << argv[3] << std::endl;
    auto Rx1 = PolygonalWireAntenna::DeSerialize( YAML::LoadFile(argv[3]) );

    auto Kern = KernelV0::NewSP();
        Kern->PushCoil( "Coil 1", Tx );
        Kern->PushCoil( "Coil 2", Rx1 );
        Kern->SetLayeredEarthEM( earth );

        Kern->SetIntegrationSize( (Vector3r() << 200, 200., 100).finished() );
        Kern->SetIntegrationOrigin( (Vector3r() << -100, -100, .5).finished() );
        Real tol(1e-13); // 13
        Kern->SetTolerance( tol ); // 1e-12

//         Kern->AlignWithAkvoDataset( YAML::LoadFile(argv[2]) );

        Kern->SetPulseDuration(0.020);
        VectorXr I(36);
        // off from VC by 1.075926340216996
        // Pulses from Wyoming Red Buttes exp 0
        I << 397.4208916184016, 352.364477036168, 313.0112765842783, 278.37896394065376, 247.81424224324982,
             220.77925043190442, 196.76493264105017, 175.31662279234038, 156.0044839325404, 138.73983004230124,
             123.42064612625474, 109.82713394836259, 97.76534468972267, 87.06061858367781, 77.56000002944572, 69.1280780096311,
             61.64250263640252, 54.99473044877554, 49.091182970515476, 43.84634004556388, 39.184136917167976, 35.03619319797924,
             31.347205894128976, 28.06346770557137, 25.139117042424758, 22.53420773366429, 20.214205433283347,
             18.144318026099942, 16.299965972298878, 14.652633628829891, 13.184271405688083, 11.870540177313893,
             10.697057141915716, 9.64778948429609, 8.709338689612677, 7.871268012862094;
        //Kern->SetPulseCurrent( VectorXr::LinSpaced( 1, 10, 200 )  ); // nbins, low, high
        Kern->SetPulseCurrent( I ); // nbins, low, high

        //VectorXr interfaces = VectorXr::LinSpaced( 41, .5, 45.5 ); // nlay, low, high
        //VectorXr interfaces = VectorXr::LinSpaced( 61, .5, 45.5 ); // nlay, low, high
        VectorXr interfaces = VectorXr::LinSpaced( 2, .5, 45.5 ); // nlay, low, high
        Real thick = .1;
        for (int ilay=1; ilay<interfaces.size(); ++ilay) {
            interfaces(ilay) = interfaces(ilay-1) + thick;
            thick *= 1.05;
        }
        Kern->SetDepthLayerInterfaces( interfaces ); // nlay, low, high

    // We could, I suppose, take the earth model in here? For non-linear that
    // may be more natural to work with?
    std::vector<std::string> tx = {std::string("Coil 1")};
    std::vector<std::string> rx = {std::string("Coil 2")};

    //std::cout << "KERNEL.yaml" << std::endl;
    //std::cout << *Kern << std::endl;

    Kern->CalculateK0( tx, rx, true ); // 3rd argument is vtk output

    std::ofstream dout = std::ofstream(std::string("Rx-")+std::string(argv[3])+std::string(".dat"));
    dout << "# Transmitters: ";
    for (auto lp : tx) {
        dout << lp << "\t";
    }
    dout << "\n# Receivers: ";
    for (auto lp : rx) {
        dout << lp << "\t";
    }
    dout << "\n# Tolerance: " << tol << std::endl;

        dout << interfaces.transpose() << std::endl;
        dout << I.transpose() << std::endl;
        dout << "#real\n";
        dout << Kern->GetKernel().real() << std::endl;
        dout << "#imag\n";
        dout << Kern->GetKernel().imag() << std::endl;
        dout.close();

    std::ofstream out = std::ofstream(std::string("Rx-")+std::string(argv[2])+std::string(".yaml"));
    //std::ofstream out = std::ofstream(std::string("k-coincident.yaml"));
    out << *Kern;
    out.close();
}


