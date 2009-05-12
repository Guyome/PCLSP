#include <stdio.h>
#include <blitz/array.h>

using namespace blitz;

int thomas(Array<double,2> alpha,Array<double,2> beta, Array<double,2> prod,
            Array<double,2> stor, Array<double,2> cons, Array<double,2> setup,
             Array<double,2> results, Array<double,1> lambda, Array<long,2> ind,
            int hor, int obj, int verbose)
{
    //problem variable
    Array<double,1> c(hor);
    Array<double,1> f(hor+1);
    Array<double,2> price(hor,hor);
    Array<double,2> demand(hor,hor);
    // algorithm variable
    float sumc;
    int t;
    if (verbose >2)
    {
        printf("j\tt\tF(t)\t\tInd\n");
        printf("-------------------------------------\n");
    }
    for(int j = 0; j < obj;  j++)
    {
        //initiate price,demand since there is no storage
        t = 0;
        price(t,t) = (alpha(j,t) + (prod(j,t) + cons(j,t)*lambda(j,t))* beta(j,t)) / (2 * beta(j,t) );
        demand(t,t) = alpha(j,t) - beta(j,t) * price(t,t);
        c(t) = demand(t,t)*( price(t,t) - prod(j,t))-setup(j,t);
        f(t) = 0;
        f(t+1) = -c(t);
        // result for the first period
        ind(j,0) = 0;
        results(j,t)=price(t,(int)ind(j,t));
        for (t = 1; t < hor;  t++)
        {
            for (int t0 = 0; t0 <= t; t0++)
            {
                // compute price
                price(t,t0) = (alpha(j,t) + (prod(j,t0) + sum(stor(j,Range(t0,t-1)))
                    + cons(j,t0)*lambda(j,t0))* beta(j,t)) / (2 * beta(j,t));
                // compute demand
                demand(t,t0) = alpha(j,t) - beta(j,t) * price(t,t0);
                // compute cost function
                sumc = 0;
                for (int i = t0; i <= t; i++)
                {
                    sumc += (prod(j,t0) + cons(j,t0)*lambda(j,t0) + sum(stor(j,Range(t0,i-1))) - price(i,t0))*demand(i,t0);
                }
                c(t0) = sumc + setup(j,t0);
            }
            // find minimal criterium
            f(t+1) = min(c(Range(0,t)) + f(Range(0,t)));
            // result for period t
            ind(j,t) = min(minIndex(c(Range(0,t))+ f(Range(0,t))));
            results(j,t)=price(t,(int)ind(j,t));
        }
        
        if (verbose >2)
        {
            for (t = 0; t < hor; t++)
            {
                printf("%d\t%d\t%f\t%d\n",j+1,t+1,-f(t+1),(int)ind(j,t)+1);
            }
        }
    }
    if (verbose>1)
    {
        printf("-------------------------------------\n");
    }
    return 0;
}
