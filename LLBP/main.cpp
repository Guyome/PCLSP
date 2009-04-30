//       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
//       version 0.1

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "FreeSolver.hpp"

using namespace Ipopt;

int main(int argv, char* argc[]) {
    int T = 4;
    float a[]={100.,100.,100.,100.};//coeff demand function
    float b[]={1.,1.,1.,1.};//coeff demand function
    float v[]={20.,20.,20.,20.};//production cost
    float h[]={2.,2.,2.,2.};//storage cost
    float r[]={18., 20., 32., 12.};

    // Create an instance of your nlp...
    SmartPtr<TNLP> mynlp = new FreeSolver(a,b,v,h,r,T);

    // Create an instance of the IpoptApplication
    SmartPtr<IpoptApplication> app = new IpoptApplication();

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
        printf("\n\n*** Error during initialization!\n");
        return (int) status;
    }

    status = app->OptimizeTNLP(mynlp);

    if (status == Solve_Succeeded) {
        // Retrieve some statistics about the solve
        Index iter_count = app->Statistics()->IterationCount();
        printf("\n\n*** The problem solved in %d iterations!\n", iter_count);
        Number final_obj = app->Statistics()->FinalObjective();
        printf("\n\n*** The final value of the objective function is %e.\n", -final_obj);
    }

    return (int) status;
}
