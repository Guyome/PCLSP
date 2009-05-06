//       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
//       version 0.1

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "LbIpopt.hpp"
#include <blitz/array.h>

using namespace Ipopt;

int ipopt(float** alpha,float** beta, float** prod,
            float** stor, float** consumption,
            float* constraint, blitz::Array<double,1> results,
            int period, int product,int verbose)
{

    // Create an instance of your nlp...
    LbIpopt* problem = new LbIpopt(alpha,beta,prod,
        stor,consumption,constraint,period,product);
    SmartPtr<TNLP> mynlp = problem;
    // Create an instance of the IpoptApplication
    SmartPtr<IpoptApplication> app = new IpoptApplication();
    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;

    // Run IpOpt
    //app->Options()->SetIntegerValue("print_level",verbose);
    //app->Options()->SetStringValue("derivative_test","None");
    status = app->Initialize();
    if (status != Solve_Succeeded) 
    {
        printf("\n\n*** Error during initialization!\n");
        return (int) 1;
    }
    
    status = app->OptimizeTNLP(mynlp);

    if (status == Solve_Succeeded) 
    {
        // Retrieve some statistics about the solve
        Number final_obj = app->Statistics()->FinalObjective();
        results(0) = -final_obj;
        printf("OK 0\n");
        for (int i = 0; i < (product+1)*period; i++)
        {
            results(i+1)=(double)problem->get_coef()[i];
        }
        printf("OK: END");
        if (verbose >0)
        {
            printf("\n\n*** The final value of the objective function is %e.\n", -final_obj);
        }
    }
    return (int) 0;
}
