#include <hxalignmicrotubules/IPOPTForRigidAlign.h>

#include <cassert>
#include <iostream>

#include <hxalignmicrotubules/IPOPTForCPD.h>

using namespace Ipopt;

namespace {

class IPOPTForRigidAlign : public IPOPTForCPD {
  public:
    /** default constructor */
    IPOPTForRigidAlign();

    /** default destructor */
    virtual ~IPOPTForRigidAlign();

    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the nlp */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, IndexStyleEnum& index_style);

    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m,
                                 Number* g_l, Number* g_u);

    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda, Number* lambda);

    /** Method to return the objective value */
    virtual bool eval_f(Index n, const Number* x, bool new_x,
                        Number& obj_value);

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
                             Number* grad_f);

    /** Method to return the constraint residuals */
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m,
                        Number* g);

    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x, Index m,
                            Index nele_jac, Index* iRow, Index* jCol,
                            Number* values);

    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not
     * NULL)
     */
    virtual bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor,
                        Index m, const Number* lambda, bool new_lambda,
                        Index nele_hess, Index* iRow, Index* jCol,
                        Number* values);

    //@}

    /** @name Solution Methods */
    //@{
    /** This method is called when the algorithm is complete so the TNLP can
     * store/write the solution */
    virtual void finalize_solution(SolverReturn status, Index n,
                                   const Number* x, const Number* z_L,
                                   const Number* z_U, Index m, const Number* g,
                                   const Number* lambda, Number obj_value,
                                   const IpoptData* ip_data,
                                   IpoptCalculatedQuantities* ip_cq);
    //@}

  private:
    /**@name Methods to block default compiler methods.
     * The compiler automatically generates the following three methods.
     *  Since the default compiler implementation is generally not what
     *  you want (for all but the most simple classes), we usually
     *  put the declarations of these methods in the private section
     *  and never implement them. This prevents the compiler from
     *  implementing an incorrect "default" behavior without us
     *  knowing. (See Scott Meyers book, "Effective C++")
     *
     */
    //@{
    //
    IPOPTForRigidAlign(const IPOPTForRigidAlign&);
    IPOPTForRigidAlign& operator=(const IPOPTForRigidAlign&);
    //@}
};

}  // namespace

IPOPTForCPD* createIPOPTForRigidAlign() { return new IPOPTForRigidAlign; }

// constructor
IPOPTForRigidAlign::IPOPTForRigidAlign() {}

// destructor
IPOPTForRigidAlign::~IPOPTForRigidAlign() {}

// returns the size of the problem
bool IPOPTForRigidAlign::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                      Index& nnz_h_lag,
                                      IndexStyleEnum& index_style) {
    n = 4;
    m = 0;

    nnz_jac_g = 0;

    nnz_h_lag = 10;
    index_style = TNLP::C_STYLE;

    return true;
}

// returns the variable bounds
bool IPOPTForRigidAlign::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                         Index m, Number* g_l, Number* g_u) {

    assert(n == 4);
    assert(m == 0);

    x_l[0] = 1.e-6;   // scale must be larger 0
    x_l[1] = -2e19;   // rho can take any value
    x_l[2] = 1.e-12;  // kappa must be larger 0
    x_l[3] = 1.e-12;  // sigma must be larger 0

    // the variables are unconstrained
    for (Index i = 0; i < 4; i++) {
        x_u[i] = 2e19;
    }
    x_u[2] = 500;

    return true;
}

// returns the initial point for the problem
bool IPOPTForRigidAlign::get_starting_point(Index n, bool init_x, Number* x,
                                            bool init_z, Number* z_L,
                                            Number* z_U, Index m,
                                            bool init_lambda, Number* lambda) {

    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    // initialize to the given starting point
    x[0] = s;
    x[1] = rho;
    x[2] = kappa;
    x[3] = sigmaSquare;

    return true;
}

// returns the value of the objective function
bool IPOPTForRigidAlign::eval_f(Index n, const Number* x, bool new_x,
                                Number& obj_value) {
    assert(n == 4);

    McDVector<double> xVec;
    xVec.insert(0, 4, x);
    obj_value = gradAndHessAndFunc.getFuncValue(xVec);

    return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IPOPTForRigidAlign::eval_grad_f(Index n, const Number* x, bool new_x,
                                     Number* grad_f) {
    assert(n == 4);
    McDVector<double> xVec;
    McDVector<double> gradient;
    xVec.insert(0, 4, x);
    gradAndHessAndFunc.gradient(xVec, gradient);
    memcpy(grad_f, gradient.dataPtr(), 4 * sizeof(double));
    return true;
}

// return the value of the constraints: g(x)
bool IPOPTForRigidAlign::eval_g(Index n, const Number* x, bool new_x, Index m,
                                Number* g) {
    assert(n == 4);
    assert(m == 0);

    return true;
}

// return the structure or values of the jacobian
bool IPOPTForRigidAlign::eval_jac_g(Index n, const Number* x, bool new_x,
                                    Index m, Index nele_jac, Index* iRow,
                                    Index* jCol, Number* values) {

    return true;
}

// return the structure or values of the hessian
bool IPOPTForRigidAlign::eval_h(Index n, const Number* x, bool new_x,
                                Number obj_factor, Index m,
                                const Number* lambda, bool new_lambda,
                                Index nele_hess, Index* iRow, Index* jCol,
                                Number* values) {
    if (values == NULL) {
        Index idx = 0;
        for (Index row = 0; row < 4; row++) {
            for (Index col = 0; col <= row; col++) {
                iRow[idx] = row;
                jCol[idx] = col;
                idx++;
            }
        }

        assert(idx == nele_hess);
    } else {
        // return the values. This is a symmetric matrix, fill the lower left
        // triangle only
        assert(n == 4);
        McDVector<double> xVec;
        McDMatrix<double> hessian;
        xVec.insert(0, 4, x);

        gradAndHessAndFunc.hessian(xVec, hessian);
        // fill the objective portion
        values[0] = hessian[0][0];  // 0,0

        values[1] = hessian[1][0];  // 1,0
        values[2] = hessian[1][1];  // 1,1

        values[3] = hessian[2][0];  // 2,0
        values[4] = hessian[2][1];  // 2,1
        values[5] = hessian[2][2];  // 2,2

        values[6] = hessian[3][0];  // 3,0
        values[7] = hessian[3][1];  // 3,1
        values[8] = hessian[3][2];  // 3,2
        values[9] = hessian[3][3];  // 3,3
    }

    return true;
}

void IPOPTForRigidAlign::finalize_solution(
    SolverReturn status, Index n, const Number* x, const Number* z_L,
    const Number* z_U, Index m, const Number* g, const Number* lambda,
    Number obj_value, const IpoptData* ip_data,
    IpoptCalculatedQuantities* ip_cq) {
    memcpy(resultValues.dataPtr(), x, sizeof(double) * 4);
}
