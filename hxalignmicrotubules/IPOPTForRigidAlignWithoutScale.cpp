#include <hxalignmicrotubules/IPOPTForRigidAlignWithoutScale.h>

#include <hxalignmicrotubules/IPOPTForCPD.h>

using namespace Ipopt;

namespace {

class IPOPTForRigidAlignWithoutScale : public IPOPTForCPD {
  public:
    /** default constructor */
    IPOPTForRigidAlignWithoutScale();

    /** default destructor */
    virtual ~IPOPTForRigidAlignWithoutScale();

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
    IPOPTForRigidAlignWithoutScale(const IPOPTForRigidAlignWithoutScale&);
    IPOPTForRigidAlignWithoutScale&
    operator=(const IPOPTForRigidAlignWithoutScale&);
    //@}
};

}  // namespace

IPOPTForCPD* createIPOPTForRigidAlignWithoutScale() {
    return new IPOPTForRigidAlignWithoutScale;
}

IPOPTForRigidAlignWithoutScale::IPOPTForRigidAlignWithoutScale() { s = 1.0; }

IPOPTForRigidAlignWithoutScale::~IPOPTForRigidAlignWithoutScale() {}

// returns the size of the problem
bool IPOPTForRigidAlignWithoutScale::get_nlp_info(Index& n, Index& m,
                                                  Index& nnz_jac_g,
                                                  Index& nnz_h_lag,
                                                  IndexStyleEnum& index_style) {
    n = 3;
    m = 0;
    nnz_jac_g = 0;
    nnz_h_lag = 6;
    index_style = TNLP::C_STYLE;

    return true;
}

// returns the variable bounds
bool IPOPTForRigidAlignWithoutScale::get_bounds_info(Index n, Number* x_l,
                                                     Number* x_u, Index m,
                                                     Number* g_l, Number* g_u) {

    assert(n == 3);
    assert(m == 0);

    x_l[0] = -2e19;   // rho can take any value
    x_l[1] = 1.e-12;  // kappa must be larger 0
    x_l[2] = 1.e-12;  // sigma must be larger 0

    // the variables are unconstrained
    for (Index i = 0; i < 3; i++) {
        x_u[i] = 2e19;
    }
    x_u[1] = 500;

    return true;
}

// returns the initial point for the problem
bool IPOPTForRigidAlignWithoutScale::get_starting_point(
    Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U,
    Index m, bool init_lambda, Number* lambda) {

    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    x[0] = rho;
    x[1] = kappa;
    x[2] = sigmaSquare;

    return true;
}

// returns the value of the objective function
bool IPOPTForRigidAlignWithoutScale::eval_f(Index n, const Number* x,
                                            bool new_x, Number& obj_value) {
    assert(n == 3);

    McDVector<double> xVec;
    xVec.append(1.0);
    xVec.insert(1, 3, x);
    obj_value = gradAndHessAndFunc.getFuncValue(xVec);

    return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IPOPTForRigidAlignWithoutScale::eval_grad_f(Index n, const Number* x,
                                                 bool new_x, Number* grad_f) {
    assert(n == 3);
    McDVector<double> xVec;
    McDVector<double> gradient;
    xVec.append(1.0);
    xVec.insert(1, 3, x);
    gradAndHessAndFunc.gradient(xVec, gradient);
    memcpy(grad_f, gradient.dataPtr() + 1, 3 * sizeof(double));
    return true;
}

// return the value of the constraints: g(x)
bool IPOPTForRigidAlignWithoutScale::eval_g(Index n, const Number* x,
                                            bool new_x, Index m, Number* g) {
    assert(n == 3);
    assert(m == 0);

    return true;
}

// return the structure or values of the jacobian
bool IPOPTForRigidAlignWithoutScale::eval_jac_g(Index n, const Number* x,
                                                bool new_x, Index m,
                                                Index nele_jac, Index* iRow,
                                                Index* jCol, Number* values) {

    return true;
}

// return the structure or values of the hessian
bool IPOPTForRigidAlignWithoutScale::eval_h(Index n, const Number* x,
                                            bool new_x, Number obj_factor,
                                            Index m, const Number* lambda,
                                            bool new_lambda, Index nele_hess,
                                            Index* iRow, Index* jCol,
                                            Number* values) {
    if (values == NULL) {
        // return the structure. This is a symmetric matrix, fill the lower left
        // triangle only.
        Index idx = 0;
        for (Index row = 0; row < 3; row++) {
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
        assert(n == 3);
        McDVector<double> xVec;
        McDMatrix<double> hessian;
        xVec.append(1.0);
        xVec.insert(1, 3, x);

        gradAndHessAndFunc.hessian(xVec, hessian);
        // fill the objective portion
        values[0] = hessian[1][1];  // 0,0

        values[1] = hessian[2][1];  // 1,0
        values[2] = hessian[2][2];  // 1,1

        values[3] = hessian[3][1];  // 2,0
        values[4] = hessian[3][2];  // 2,1
        values[5] = hessian[3][3];  // 2,2
    }

    return true;
}

void IPOPTForRigidAlignWithoutScale::finalize_solution(
    SolverReturn status, Index n, const Number* x, const Number* z_L,
    const Number* z_U, Index m, const Number* g, const Number* lambda,
    Number obj_value, const IpoptData* ip_data,
    IpoptCalculatedQuantities* ip_cq) {

    memcpy(resultValues.dataPtr() + 1, x, sizeof(double) * 3);
    resultValues[0] = 1.0;
}
