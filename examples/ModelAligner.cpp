/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      08/30/2017 04:08:53 AM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2017, University of Utah
 * @copyright Copyright (c) 2017, Trevor Irons & Lemma Software, LLC
 */

#include <Merlin>
using namespace Lemma;

int main(int argc, char** argv) {

    if (argc<5) {
        std::cout << "ModelAligner aligns a dummy model with a pre-calculated"
                  << "imaging kernel.\n\n"
                  << "./ModelAligner Kernel.yaml T2Low T2High nT2  \n";

        return(EXIT_FAILURE);
    }

    auto Kernel = KernelV0::DeSerialize(YAML::LoadFile(argv[1]));

    auto Model = LayeredEarthMR::NewSP();
        Model->AlignWithKernel(Kernel);
        Model->SetT2StarBins(10, 500, 20);

    std::cout << *Model << std::endl;
}


