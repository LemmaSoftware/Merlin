/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      08/28/2017 09:07:04 AM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2017, University of Utah
 * @copyright Copyright (c) 2017, Lemma Software, LLC
 */

#include <Merlin>
using namespace Lemma;

int main(int argc, char** argv) {

    if (argc<5) {
        std::cout << "./ForwardFID Kernel.yaml TxString RxString  vtkoutput<true/false> \n";
    //    return(EXIT_FAILURE);
    }

    auto Forward = ForwardFID::NewSP();
        Forward->SetWindowEdges( VectorXr::LinSpaced(10,0,1) );

    std::cout << *Forward << std::endl;

    //auto FID = Forward->ForwardModel();
    //std::cout << *FID << std::endl;

    return EXIT_SUCCESS;
}
