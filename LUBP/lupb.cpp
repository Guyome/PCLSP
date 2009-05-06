#include <stdio.h>
#include <blitz/array.h>

using namespace blitz;

int ipopt(Array<double,2> alpha,Array<double,2> beta, Array<double,2> prod,
            Array<double,2> stor, Array<double,2> consumption,
            Array<double,1> lambda, blitz::Array<double,1> results,
            int period, int product,int verbose)
{
    //problem variable
    Array<double,1> c(hor);
    Array<double,1> f(hor);
    Array<double,2> price(hor,hor);
    Array<double,2> demand(hor,hor);
    // algorithm variable
    Array<int,1> ind;
    int k,i,t0,t = 0,;
    float sumstor,sum,summin;
    
    //initiate price,demand since there is no storage
    for(int j = 0; j < obj;  j++)
    {
        price(t,t) = (alpha(j,t) + (prod(j,t) + cons(j,t)*lamdba(j,t))* beta(j,t)) / (2 * beta(j,t) );
        demand(t,t) = alpha(j,t) - beta(j,t) * price(t,t);
        c(t)= demand(t,t)*( price(t,t) - prod(j,t))-setup(j,t);
        f(t+1) = -c(t);
        ind(0) = 0;
        for ( t = 1; t < hor;  t++)
        {
            for (t0 = 0; t0 <= t; t0++)
            {
                // sum of first storage
                sumstor = 0;
                for (i = t0; i < t; i++)
                {
                    sumstor += stor(j,i);
                }
                // compute price and demand
                price(t,t0) = (alpha(j,t) + 
                    (prod(j,t0) + sumstor + cons(j,t0)*lamdba(j,t0))* beta(j,t))
                     / (2 * beta(j,t));
                demand(t,t0) = alpha(j,t) - beta(j,t) * price(t,t0);
            }
            for (t0 = 0; t0 <= t; t0++)
            {
                sum = 0;
                for (i = t0; i <= t; i++)
                {
                    sumstor = 0;
                    for (k = t0; k < i; k++)
                    {
                        sumstor += stor(j,i);
                    }
                    sum += (prod(j,t0) + sumstor - price(j,t0))*demand(j,t0);
                }
                c(t0) = sum + setup(j,t0);
            }
            // compute minimal criterium
            summin = 0;
            for (t0 = 0; t0 <= t; t0++)
            {
                if (c(t0) + f(t0) < summin)
                {
                    summin = c(t0) + f(t0);
                    ind(t) = t0;
                }
            }
            f(t+1) = summin;
        }
        for( t = 0; t < hor;  t++)
        {
            optiprice(j,t)=price(t,ind(t));
        }
        
        printf("t \t F(t) \t Ind\n");
        printf("--------------------------------\n");
        for (t = 0; t < hor; t++)
        {
            printf("%d\t%f\t%d\n",t+1,-f(t+1),ind(t)+1);
        }
        printf("--------------------------------\n");
    }
    return 0;
}

