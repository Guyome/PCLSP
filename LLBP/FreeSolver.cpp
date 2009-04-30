//       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
//       version 0.1

#include "FreeSolver.hpp"
#include <iostream>

using namespace Ipopt;

// constructor
FreeSolver::FreeSolver(float* _alpha,float* _beta, float* _prod,float* _stor,
                            float* _constraint, float* consumption, int _period, int _product);
{
    period = _period;
    product = _product;
    alpha = new float[period][product];
    beta = new float[period][product];
    prod = new float[period][product];
    stor = new float[period][product];
    constraint = new float[period];
    consumption = new float[period][product];
    for (int i = 0; i < period; i++)
    {
        for (int j = 0; j < product; j ++)
        {
            alpha[i][j]=_alpha[i][j];
            std::cout << "alpha["<< i<<"]["<< j<<"]="<<alpha[i][j]<<"\t";
            beta[i][j]=_beta[i][j];
            std::cout << "beta["<< i<<"]["<< j<<"]="<<beta[i][j]<<"\t";
            prod[i][j]=_prod[i];
            std::cout << "prod["<< i<<"]["<< j<<"]="<<prod[i][j]<<"\t";
            stor[i][j]=_stor[i][j];
            std::cout << "stor["<< i<<"]["<< j<<"]="<<stor[i][j]<<std::endl;
            consumption[i]=_consumption[i];
            std::cout << "consumption["<< i<<"]="<<consumption[i] << std::endl;
        }
        constraint[i]=_constraint[i];
        std::cout << "constraint["<< i<<"]="<<constraint[i]<<"\t";
    }
}

//destructor
FreeSolver::~FreeSolver(){}

// returns the size of the problem
bool FreeSolver::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style) 
{
    n = 3*period*product; //number of variable
    m = (product+1)period; //number of constraint
    nnz_jac_g = 5*period*product-1;//number of non zero value in constraint jacobian
    nnz_h_lag = product*period;//number of non zero value in lagrangian hessian
    index_style = TNLP::C_STYLE;// use the C style indexing (0-based)
    return true;
}

// returns the variable bounds
bool FreeSolver::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u) 
{
    // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
    // If desired, we could assert to make sure they are what we think they are.
    assert(n == 3*period*product);
    assert(m == (product+1)*period);
    
    for (Index i=0; i<n; i++) 
    {
        x_l[i] = 0.0; // the variables are positives
        x_u[i] = 2e19;  // have no upper bounds
    }
    for (Index i=0; i<m; i++) 
    {
        g_l[i] = g_u[i] = alpha[i]; // product*priod equality contraint
        if (i >= period*product)
        {
            g_l[i] = 0.0;
            g_u[i] = constraint[i-period*product]
        }
    }
    return true;
}

// returns the initial point for the problem
bool FreeSolver::get_starting_point(Index n, bool init_x, Number* x,
                            bool init_z, Number* z_L, Number* z_U,
                            Index m, bool init_lambda,
                            Number* lambda)
    {
    // Here, we assume we only have starting values for x, if you code
    // your own NLP, you can provide starting values for the dual variables
    // if you wish
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);
    // initialize to the given starting point
    for (Index i=0; i<period; i++) {
        x[i] = constraint[i]; //production is equal to constraint
    }
    for(Index i=period;i<2*period;i++){
        x[i]=(alpha[i-period]-constraint[i-period])/beta[i-period]; //price to stay in feasible region
    }
    for(Index i=2*period;i<3*period;i++){
        x[i]=0.;//no storage
    }
    std::cout <<"*** INITIAL POINT\n";
    std::cout <<"Prod:\t"<<x[0]<<"\\"<<x[1]<<"\\"<<x[2]<<"\\"<<x[3]<<"\n";
    std::cout <<"Stor:\t"<<x[8]<<"\\"<<x[9]<<"\\"<<x[10]<<"\\"<<x[11]<<"\n";
    std::cout <<"Price:\t"<<x[4]<<"\\"<<x[5]<<"\\"<<x[6]<<"\\"<<x[7]<<"\n";
    return true;
}

// returns the value of the objective function
bool FreeSolver::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) 
{
    assert(n == 3*period*product);
    obj_value = 0.0;
    for(Index i=0; i< period; i++) 
    {
        for (Index j = 0; j < product; j++)
        {
            obj_value -= (alpha[i][j]-beta[i][j]*x[period+i])*x[period+i+j]-prod[i][j]*x[i+j]-stor[i][j]*x[2*period+i+j];
        }
    }
    return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool FreeSolver::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) 
{
    assert(n == 3*period*product);
    for(Index i=0; i<period;i++)
    {
        for (Index j = 0; j < product; j += 1)
        {
        grad_f[i] = prod[i][j];
        grad_f[period+i] = 2*beta[i][j]*x[period+i+j]-alpha[i][j];
        grad_f[2*period+i] = stor[i][j];
        }
    }
    return true;
}

// return the value of the constraints: g(x)
bool FreeSolver::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    assert(n == 3*period*product);
    assert(m == (product+1)*period);
    int idx = 0;
    //equatlity constraint
    g[idx] = x[0] - x[2*period] + beta[0]*x[period]; //no storage form time period 0
    for(Index i=1; i<period;i++)
    {
        for (Index j = 0; j < object; j ++)
        {
            g[idx] = x[i+j] + x[2*period+i-1+j] - x[2*period+i+j] + beta[i][j]*x[period+i+j];
            idx++;
        }
        idx++;
    }
    // inequality constraint
    for (Index i = 0; i < period; i++)
    {
        g[idx] = 0;
        for (Index j = 0; j < product; j ++)
        {
            g[idx] += consumption[i][j]*x[i+j];
        }
        idx++;
    }
    assert(idx == m);
    return true;
}

// return the structure or values of the jacobian
bool FreeSolver::eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow, Index *jCol,
                            Number* values) 
{
    if (values == NULL) 
    {
    // return the structure of the jacobian
        //element on 3 diagonals
        Index idx=0;
        for(Index j=0;j<3;j++)
        {
            for(Index i=0; i<period*product;i++) 
            {
                iRow[idx]=i;
                jCol[idx]=i+j*period*product;
                idx++;
            }
        }
                // element due to the inequality constraint
        for(Index i=0; i<period*product;i++)
        {
            iRow[idx]=i+period*product;
            jCol[idx]=i;
            idx++;
        }
        // elements on the sub-diagonale (i_{t-1})
        for(Index i=1;i<period*product;i++)
        {
            iRow[idx]=i;
            jCol[idx]=i-1+2*period*product;
            idx++;
        }
        assert(idx == nele_jac);
    }
    else 
    {
    // return the values of the jacobian of the constraints
        for(Index i=0;i<period;i++)
        {
            for (Index j = 0; j < product; j++)
            {
                values[i+j]=1;
                values[i+j+period=beta[i][j];
                values[i+j+2*period]=-1;
                values[i+j+3*period]=alpha[i][j];
                if (i<period-1)
                {
                    values[i+j+4*period]=1;
                }
            }
        }
    }
    return true;
}

//return the structure or values of the hessian
bool FreeSolver::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
    if (values == NULL)
    {
        // the hessian for this problem is actually dense
        Index idx=0;
        for (Index i = period*product; i <2*period*product; i++) {
            iRow[idx] = i;
            jCol[idx] = i;
            idx++;
        }
        assert(idx == nele_hess);
    }
    else
    {
        Index idx=0;
        for (Index i = 0; i <period; i++)
        {
            for (Index j = 0; j < product; j++)
            {
                values[idx] = 2*beta[i][j]*obj_factor;
                idx++;
            }
            idx++;
        }
    }
    return true;
}

void FreeSolver::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
              const IpoptData* ip_data,
              IpoptCalculatedQuantities* ip_cq) {

    std::cout <<"*** RESULTS:\n";
    std::cout <<"Dem:\t"<<alpha[0]-beta[0]*x[4]<<"\\"<<alpha[1]-beta[1]*x[5]<<"\\"<<alpha[2]-beta[2]*x[6]<<"\\"<<alpha[3]-beta[3]*x[7]<<"\n";
    std::cout <<"Prod:\t"<<x[0]<<"\\"<<x[1]<<"\\"<<x[2]<<"\\"<<x[3]<<"\n";
    std::cout <<"Stor:\t"<<x[8]<<"\\"<<x[9]<<"\\"<<x[10]<<"\\"<<x[11]<<"\n";
    std::cout <<"Price:\t"<<x[4]<<"\\"<<x[5]<<"\\"<<x[6]<<"\\"<<x[7]<<"\n";
}
