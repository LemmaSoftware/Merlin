/* This file is part of Lemma, a geophysical modelling and inversion API.
 * More information is available at http://lemmasoftware.org
 */

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/**
 * @file
 * @date      12/01/2017 10:29:26 PM
 * @version   $Id$
 * @author    Trevor Irons (ti)
 * @email     tirons@egi.utah.edu
 * @copyright Copyright (c) 2017, University of Utah
 * @copyright Copyright (c) 2017, Trevor Irons & Lemma Software, LLC
 */

#include "Merlin"

using namespace Lemma;

int main(int argc, char** argv) {
    if (argc<2) {
        std::cout << "useage\n"
                  << "./KernelAligner  <AkvoData> \n";
        exit(EXIT_SUCCESS);
    }
    auto K0 = KernelV0::NewSP();

    std::cout << *K0 << std::endl;
}
