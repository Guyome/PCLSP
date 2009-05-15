//       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
//       version 0.1

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "LbIpopt.hpp"
#include <blitz/array.h>

using namespace Ipopt;
using namespace blitz;

int ipopt(Array<double,2> alpha,Array<double,2> beta, Array<double,2> prod,
            Array<double,2> stor, Array<double,2> consumption, Array<double,2> setup,
            Array<double,1> constraint, Array<double,1> results,
            int period, int product,int verbose)
{

    // Create an instance of your nlp...
    LbIpopt* problem = new LbIpopt(alpha,beta,prod,stor,
        consumption,setup,constraint,period,product);
    SmartPtr<TNLP> mynlp = problem;
    // Create an instance of the IpoptApplication
    SmartPtr<IpoptApplication> app = new IpoptApplication();
    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;

    // Run IpOpt
    app->Options()->SetIntegerValue("print_level", verbose);
    app->Options()->SetStringValue("derivative_test","None");
    status = app->Initialize();
    if (status != Solve_Succeeded) 
    {
        printf("\n\n*** Error during initialization!\n");
        return status;
    }
    
    status = app->OptimizeTNLP(mynlp);

    if (status == Solve_Succeeded) 
    {
        // Retrieve some statistics about the solve
        Number final_obj = app->Statistics()->FinalObjective();
        int idx = 0;
        results(idx) = -final_obj;
        idx++;
        for (int i = product*period; i < (product+1)*period; i++)
        {
            results(idx)=problem->get_coef()(i);
            idx++;
        }
        for (int i = 0; i < 3*period*product; i++)
        {
            results(idx)=problem->get_final_value()(i);
            idx++;
        }
        if (verbose > 1)
        {
            printf("\n\n*** The final value of the objective function is %e.\n", -final_obj);
        }
    }
    return status;
}
