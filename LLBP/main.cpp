//       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
//       version 0.1

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "FreeSolver.hpp"

using namespace Ipopt;

int main(int argv, char* argc[]) {
    int T = 4;//number of time T
    int J = 1;//number of product
    float r[]={18., 20., 32., 12.};//constraint
    float** al = new float*[T];//coeff demand function
    float** b = new float*[T] ;//coeff demand function
    float** v = new float*[T];//production cost
    float** h = new float*[T];//storage cost
    float** a = new float*[T]; //consumption
    for (int i = 0; i < T; i ++)
    {
        al[i] = new float [J];
        b[i] = new float [J] ;
        v[i] = new float [J];
        h[i] = new float [J];
        a[i] = new float [J];
        for (int j = 0; j < J; j ++)
        {
            al[i][j] = 100.;
            b[i][j] = 1.;
            v[i][j] = 20.;
            h[i][j] = 2.;
            a[i][j] = 1;
        }
        
    }

    // Create an instance of your nlp...
    SmartPtr<TNLP> mynlp = new FreeSolver(al,b,v,h,a,(float*)r,T,J);

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
