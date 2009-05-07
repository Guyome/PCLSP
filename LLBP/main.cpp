//       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
//       version 0.1

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "LbIpopt.hpp"
#include <blitz/array.h>

using namespace Ipopt;
using namespace blitz;

int ipopt(Array<double,2> alpha,Array<double,2> beta, Array<double,2> prod,
            Array<double,2> stor, Array<double,2> consumption,
            Array<double,1> constraint, Array<double,1> results,
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
    app->Options()->SetIntegerValue("print_level",5*verbose);
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
        results(0) = -final_obj;
        for (int i = 0; i < (product+1)*period; i++)
        {
            results(i+1)=problem->get_coef()(i);
        }
        if (verbose > 0)
        {
            printf("\n\n*** The final value of the objective function is %e.\n", -final_obj);
        }
    }
    return status;
}

int main ()
{
    int hor = 3;
    int obj = 1;
    int status;
    Array<double,2> alpha(hor,obj);
    Array<double,2> beta(hor,obj);
    Array<double,2> prod(hor,obj);
    Array<double,2> stor(hor,obj);
    Array<double,2> consumption(hor,obj);
    Array<double,1> constraint(hor);
    Array<double,1> results((obj+1)*hor+1);

    alpha = 70, 50, 80;
    beta = 1, 3, 2;
    prod = 10, 10, 10;
    stor = 1, 1, 1;
    consumption = 2, 2, 2;
    results = 0, 0,0;
    constraint = 11, 15, 11;

    
    status = ipopt(alpha, beta, prod,
            stor, consumption, constraint,
            results, hor, obj, 1);
    for (int i = 0; i < hor; i += 1)
    {
            printf("%f\n",results(0,i));
    }
    return 0;
}
