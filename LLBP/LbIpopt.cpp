//       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
//       version 0.1

#include "LbIpopt.hpp"
#include <blitz/array.h>

using namespace Ipopt;
using namespace blitz;

// constructor
LbIpopt::LbIpopt(Array<double,2> _alpha, Array<double,2> _beta, Array<double,2> _prod,
    Array<double,2> _stor, Array<double,2> _consumption, Array<double,2> _setup,
    Array<double,1> _constraint, int _period, int _product)
{
    period = _period;
    product = _product;
    alpha = new Array<double,2>(_alpha);
    beta = new Array<double,2>(_beta);
    prod = new Array<double,2>(_prod);
    stor = new Array<double,2>(_stor);
    consumption = new Array<double,2>(_consumption) ; 
    setup =  new Array<double, 2>(_setup);
    constraint = new Array<double,1>(_constraint);
    coef = new Array<double,1>((product+1)*period) ;
}

//destructor
LbIpopt::~LbIpopt(){}

// returns the size of the problem
bool LbIpopt::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style) 
{
    n = 3*period*product; //number of variable
    m = (product+1)*period; //number of constraint
    nnz_jac_g = 5*period*product-1;//number of non zero value in constraint jacobian
    nnz_h_lag = product*period;//number of non zero value in lagrangian hessian
    index_style = TNLP::C_STYLE;// use the C style indexing (0-based)
    return true;
}

// returns the variable bounds
bool LbIpopt::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u) 
{
    // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
    // If desired, we could assert to make sure they are what we think they are.
    assert(n == 3*period*product);
    assert(m == (product+1)*period);
    for (Index i = 0; i < period; i ++)
    {
        for (Index j = 0; j < product; j ++)
        {
            x_l[j+i*product] = 0.0; // the variables are positives
            x_u[j+i*product] = 2e19*(*setup)(j,i);  // null in function of setup structure
        }
    }
    for (Index i=period*product; i<n; i++) 
    {
        x_l[i] = 0.0; // the variables are positives
        x_u[i] = 2e19;  // have no upper bounds
    }
    // price upper bound
    Index idx=period*product;
    for (Index i = 0; i < period; i ++)
    {
        for ( Index j = 0; j < product; j ++)
        {
            x_u[idx] = (*alpha)(j,i)/(*beta)(j,i);
            idx++;
        }
    }
    idx=0;
    // product*period equality contraint
    for (Index i=0; i<period; i++) 
    {
        for (Index j = 0; j < product; j ++)
        {
            g_l[idx] = g_u[idx] = (*alpha)(j,i); 
            idx++;
        }
    }
    //inequality constraints
    for (Index i=0; i<period; i++) 
    {
        g_l[idx] = 0.0;
        g_u[idx] = (*constraint)(i);
        idx++;
    }
	assert(m == idx );
    return true;
}

// returns the initial point for the problem
bool LbIpopt::get_starting_point(Index n, bool init_x, Number* x,
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
    //find the minimal of consumption per time period
    Index idx = 0;
    for ( Index i = 0; i < period; i++)
    {
        for (Index j = 0; j < product; j++)
        {
            // initialize to the given starting point
            x[idx] = (*constraint)(i)/(max((*consumption)(Range::all(),i))*product); //production minimize the constraint
            //price to stay in feasible region
            x[idx+period*product]=((*alpha)(j,i)-x[idx])/(*beta)(j,i);
            x[idx+2*period*product]=0.;//no storage
            idx++;
        }
    }
    
    assert(idx == n);
    return true;
}

// returns the value of the productective function
bool LbIpopt::eval_f(Index n, const Number* x, bool new_x, Number& product_value) 
{
    assert(n == 3*period*product);
    product_value = 0.0;
    Index idx = 0;
    for(Index i=0; i< period; i++) 
    {
        for (Index j = 0; j < product; j++)
        {
            product_value -= ((*alpha)(j,i)-(*beta)(j,i)*x[period*product+idx])*x[period*product+idx]-(*prod)(j,i)*x[idx]-(*stor)(j,i)*x[2*period*product+idx];
            idx++;
        }
    }
    return true;
}

// return the gradient of the productective function grad_{x} f(x)
bool LbIpopt::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) 
{
    assert(n == 3*period*product);
    Index idx = 0;
    for(Index i=0; i<period;i++)
    {
        for (Index j = 0; j < product; j += 1)
        {
        grad_f[idx] = (*prod)(j,i);
        grad_f[idx+period*product] = 2*(*beta)(j,i)*x[period*product+idx]-(*alpha)(j,i);
        grad_f[idx+2*period*product] = (*stor)(j,i);
        idx++;
        }
    }
    return true;
}

// return the value of the constraints: g(x)
bool LbIpopt::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    assert(n == 3*period*product);
    assert(m == (product+1)*period);
    Index idx = 0;
    //equatlity constraint
    for(Index i=0; i<period;i++)
    {
        for (Index j = 0; j < product; j ++)
        {
            if (i > 0)
            {
                 g[idx] = x[idx]  + x[2*period*product+idx-1] - x[2*period*product+idx] + (*beta)(j,i)*x[period*product+idx];
            }
            else
            {
                g[idx] = x[idx] - x[2*period*product+idx] + (*beta)(j,i)*x[period*product+idx];
            }
            idx++;
        }
    }
    // inequality constraint
    for (Index i = 0; i < period; i++)
    {
        g[idx] = 0;
        for (Index j = 0; j < product; j ++)
        {
            g[idx] += (*consumption)(j,i)*x[j+i*product];
        }
        idx++;
    }
    assert(idx == m);
    return true;
}

// return the structure or values of the jacobian
bool LbIpopt::eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow, Index *jCol,
                            Number* values) 
{
    if (values == NULL) 
    {
    // return the structure of the jacobian
        //element on 3 diagonals
        Index idx=0;
        for(Index k=0;k<3;k++)
        {
            for(Index i=0; i<period*product;i++) 
            {
                iRow[idx]=i;
                jCol[idx]=i+k*period*product;
                idx++;
            }
        }
        // element due to the inequality constraint
        for(Index i=0; i<period;i++)
        {
            for (Index j = 0; j <product; j ++)
            {
                iRow[idx]=i+period*product;
                jCol[idx]=j+i*product;
                idx++;
            }
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
        Index idx =0;
        for(Index i=0;i<period;i++)
        {
            for (Index j = 0; j < product; j++)
            {
                //element on 3 diagonals
                values[idx]=1.;
                values[idx+period*product]=(*beta)(j,i);
                values[idx+2*period*product]=-1.;
                // element due to the inequality constraint
                values[idx+3*period*product]=(*consumption)(j,i);
                // elements on the sub-diagonale (i_{t-1})
                if (i<period-1)
                {
                    values[idx+4*period*product]=1.;
                }
                idx++;
            }
        }
    }
    return true;
}

//return the structure or values of the hessian
bool LbIpopt::eval_h(Index n, const Number* x, bool new_x,
                   Number product_factor, Index m, const Number* lambda,
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
        for (Index i = 0; i < period; i++)
        {
            for (Index j = 0; j < product; j++)
            {
                values[idx] = 2*(*beta)(j,i)*product_factor;
                idx++;
            }
        }
        assert(idx == nele_hess);
    }
    return true;
}

void LbIpopt::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number product_value,
              const IpoptData* ip_data,
              IpoptCalculatedQuantities* ip_cq) 
{
    for (int i = 0; i < m; i ++)
    {
        (*coef)(i) = lambda[i];
    }
}

Array<double,1> LbIpopt::get_coef(){
    return *coef;
}
