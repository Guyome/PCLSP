#include <stdio.h>

int main(int argc, char *argv[])
{
    //problem parameter
    int hor = 3;
    float alpha[]={50, 50, 100};
    float beta[]={2, 2, 2};
    float setup[]={11, 15, 11};
    float prod[]={10, 10, 10};
    float stor[]={1, 1, 1};
    //problem variable
    float c[]={0, 0, 0};
    float f[]={0, 0, 0, 0};
    float price[hor][hor];
    float demand[hor][hor];
    // algorithm variable
    int k,j,t0,t = 0,ind[hor];
    float sumstor,sum,summin;


    //initiate price,demand since there is no storage
    price[t][t] = (alpha[t] + prod[t] * beta[t]) / (2 * beta[t] );
    demand[t][t] = alpha[t] - beta[t] * price[t][t];
    c[t]= demand[t][t]*( price[t][t] - prod[t])-setup[t];
    f[t+1] = -c[t];
    ind[0] = 0;
    for ( t = 1; t < hor;  t++)
    {
        for (t0 = 0; t0 <= t; t0++)
        {
            // sum of first storage
            sumstor = 0;
            for (j = t0; j < t; j++)
            {
                sumstor += stor[j];
            }
            // compute price and demand
            price[t][t0] = (alpha[t] + (prod[t0] + sumstor)* beta[t]) / (2 * beta[t] );
            demand[t][t0] = alpha[t] - beta[t] * price[t][t0];
        }
        for (t0 = 0; t0 <= t; t0++)
        {
            sum = 0;
            for (j = t0; j <= t; j++)
            {
                sumstor = 0;
                for (k = t0; k < j; k++)
                {
                    sumstor += stor[j];
                }
                sum += (prod[t0] + sumstor - price[j][t0])*demand[j][t0];
            }
            c[t0] = sum + setup[t0];
        }
        // compute minimal criterium
        summin = 0;
        for (t0 = 0; t0 <= t; t0++)
        {
            if (c[t0] + f[t0] < summin)
            {
                summin = c[t0] + f[t0];
                ind[t] = t0;
            }
        }
        f[t+1] = summin;
    }
    printf("t \t F[t] \t Ind\n");
    printf("--------------------------------\n");
    for (t = 0; t < hor; t++)
    {
        printf("%d\t%f\t%d\n",t+1,-f[t+1],ind[t]+1);
    }
    printf("--------------------------------\n");
    return 0;
}


