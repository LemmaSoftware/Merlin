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

    if (argc<5) {
        std::cout << "./KernelV0-2 Kernel.yaml TxString RxString  vtkoutput<true/false> \n";
        return(EXIT_FAILURE);
    }

    std::cout << "Using kernel paramaters: " << argv[1] << std::endl;
    auto Kern = KernelV0::DeSerialize( YAML::LoadFile(argv[1]) );
    std::cout << "Kernel DeSerialized successful" << std::endl;

    std::vector<std::string> tx = {std::string(argv[2])};
    std::vector<std::string> rx = {std::string(argv[3])};

    std::cout << "argv[4]\t" << argv[4] << std::endl;
    if( std::string(argv[4]) == "true" || std::string(argv[4]) == "True") {
        std::cout << "Using VTK, output files may be very large" << std::endl;
        Kern->CalculateK0( tx, rx, true ); // 3rd argument is vtk output
    } else {
        std::cout << "not using VTK" << std::endl;
        Kern->CalculateK0( tx, rx, false ); // 3rd argument is vtk output
    }

    // TODO fix python post-processing so this is not necessary
    // Save in simplified format for easy python plotting
//     std::ofstream dout = std::ofstream(std::string("Tx")+std::string(argv[2])+std::string("Rx-")+std::string(argv[3])+std::string(".dat"));
//     dout << "# Transmitters: ";
//     for (auto lp : tx) {
//         dout << lp << "\t";
//     }
//     dout << "\n# Receivers: ";
//     for (auto lp : rx) {
//         dout << lp << "\t";
//     }
//     dout << "\n# Tolerance: " << Kern->GetTolerance() << std::endl;
//         dout << Kern->GetInterfaces().transpose() << std::endl;
//         dout << Kern->GetPulseDuration()*Kern->GetPulseCurrent().transpose() << std::endl;
//         dout << "#real\n";
//         dout << Kern->GetKernel().real() << std::endl;
//         dout << "#imag\n";
//         dout << Kern->GetKernel().imag() << std::endl;
//         dout.close();

    // Save YAML kernel
    std::ofstream out = std::ofstream(std::string("Tx-")+std::string(argv[2])+std::string("_Rx-")+std::string(argv[3])+std::string(".yaml"));
    out << *Kern;
    out.close();

    return EXIT_SUCCESS;
}


