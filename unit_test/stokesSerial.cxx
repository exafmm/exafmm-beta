/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include "serialfmm.h"

int main(int ac, char** av)
{
    int numBodies = 10000;
    int numTarget = 100;
    IMAGES = 0;
    THETA = 1 / sqrtf(4);
    Bodies bodies, jbodies;
    Cells cells, jcells;
    SerialFMM<Stokes> FMM;

    FMM.initialize();
    FMM.setDelta(.01);

    for ( int it = 0; it != 5; ++it )
    {
        numBodies = int(pow(10, (it + 24) / 8.0));
        std::cout << "N             : " << numBodies << std::endl;
        bodies.resize(numBodies);
        FMM.cube(bodies, 1, 1);
        FMM.startTimer("FMM          ");
        FMM.setDomain(bodies);
        cells.clear();
#ifdef TOPDOWN
        FMM.topdown(bodies, cells);
#else
        FMM.bottomup(bodies, cells);
#endif
        jcells = cells;
        FMM.downward(cells, jcells);
        FMM.stopTimer("FMM          ", true);
        FMM.eraseTimer("FMM          ");
        FMM.startTimer("Direct sum   ");
        FMM.buffer = bodies;
#if 1
        FMM.initTarget(FMM.buffer);
        if ( IMAGES != 0 )
        {
            jbodies = FMM.periodicBodies(FMM.buffer);
        }
        else
        {
            jbodies = FMM.buffer;
        }
        FMM.buffer.resize(numTarget);
        FMM.evalP2P(FMM.buffer, jbodies, 1);
        FMM.writeTarget(FMM.buffer);
#else
        FMM.readTarget(FMM.buffer);
#endif
        FMM.stopTimer("Direct sum   ", true);
        FMM.eraseTimer("Direct sum   ");
        FMM.resetTimer();

        real diff = 0, norm = 0;
        bodies.resize(numTarget);
        FMM.evalError(bodies, FMM.buffer, diff, norm);
        FMM.printError(diff, norm);
    }
    FMM.finalize();
    std::cout << "v1 = [";
    for (B_iter B = bodies.begin(), end = bodies.end(); B != end; ++B)
        std::cout << B->TRG[0] << "," << B->TRG[1] << "," << B->TRG[2] << ";";
    std::cout << "];" << std::endl;
    std::cout << "v2 = [";
    for (B_iter B = FMM.buffer.begin(), end = FMM.buffer.end(); B != end; ++B)
        std::cout << B->TRG[0] << "," << B->TRG[1] << "," << B->TRG[2] << ";";
    std::cout << "];" << std::endl;
    return 0;
}
