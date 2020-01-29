#ifndef __MMA__
#define __MMA__

#include <petsc.h>

/*
Copyright (C) 2013-2019, Niels Aage
*/
class MMA {
  public:
    // Construct using defaults subproblem penalization
    MMA(PetscInt n, PetscInt m, Vec x);
    // User defined subproblem penalization
    MMA(PetscInt n, PetscInt m, Vec x, PetscScalar* a, PetscScalar* c, PetscScalar* d);
    // Initialize with restart from itr
    MMA(PetscInt n, PetscInt m, PetscInt itr, Vec xo1, Vec xo2, Vec U, Vec L);
    // Initialize with restart and specify subproblem parameters
    MMA(PetscInt n, PetscInt m, PetscInt itr, Vec xo1, Vec xo2, Vec U, Vec L, PetscScalar* a, PetscScalar* c,
        PetscScalar* d);
    // Destructor
    ~MMA();

    // Set and solve a subproblem: return new xval
    PetscErrorCode Update(Vec xval, Vec dfdx, PetscScalar* gx, Vec* dgdx, Vec xmin, Vec xmax);

    // Return necessary data for possible restart
    PetscErrorCode Restart(Vec xo1, Vec xo2, Vec U, Vec L);

    // Set the aggresivity of the moving asymptotes
    PetscErrorCode SetAsymptotes(PetscScalar init, PetscScalar decrease, PetscScalar increase);

    // do/don't add convexity approx to constraints: default=false
    PetscErrorCode ConstraintModification(PetscBool conMod) {
        constraintModification = conMod;
        return 0;
    };

    // val=0: default, val=1: increase robustness, i.e
    // control the spacing between L < alp < x < beta < U,
    PetscErrorCode SetRobustAsymptotesType(PetscInt val);

    // Sets outer movelimits on all primal design variables
    // This is often requires to prevent the solver from oscilating
    PetscErrorCode SetOuterMovelimit(PetscScalar Xmin, PetscScalar Xmax, PetscScalar movelim, Vec x, Vec xmin,
                                     Vec xmax);

    // Return KKT residual norms (norm2 and normInf)
    PetscErrorCode KKTresidual(Vec xval, Vec dfdx, PetscScalar* gx, Vec* dgdx, Vec xmin, Vec xmax, PetscScalar* norm2,
                               PetscScalar* normInf);

    // Inf norm on diff between two vectors: SHOULD NOT BE HERE - USE BASIC
    // PETSc!!!!!
    PetscScalar DesignChange(Vec x, Vec xold);

  private:
    // Set up the MMA subproblem based on old x's and xval
    PetscErrorCode GenSub(Vec xval, Vec dfdx, PetscScalar* gx, Vec* dgdx, Vec xmin, Vec xmax);

    // Interior point solver for the subproblem
    PetscErrorCode SolveDIP(Vec xval);

    // Compute primal vars based on dual solution
    PetscErrorCode XYZofLAMBDA(Vec x);

    // Dual gradient
    PetscErrorCode DualGrad(Vec x);

    // Dual Hessian
    PetscErrorCode DualHess(Vec x);

    // Dual line search
    PetscErrorCode DualLineSearch();

    // Dual residual
    PetscScalar DualResidual(Vec x, PetscScalar epsi);

    // Problem size and iteration counter
    PetscInt n, m, k;

    // "speed-control" for the asymptotes
    PetscScalar asyminit, asymdec, asyminc;

    // do/don't add convexity constraint approximation in subproblem
    PetscBool constraintModification; // default = FALSE

    // Bool specifying if non lin constraints are included or not
    PetscBool NonLinConstraints;

    // 0: (default) span between alp L x U beta,
    // 1: increase the span for further robustness
    PetscInt RobustAsymptotesType;

    // Local vectors: penalty numbers for subproblem
    PetscScalar *a, *c, *d;

    // Local vectors: elastic variables
    PetscScalar* y;
    PetscScalar  z;

    // Local vectors: Lagrange multipliers:
    PetscScalar *lam, *mu, *s;

    // Global: Asymptotes, bounds, objective approx., constraint approx.
    Vec L, U, alpha, beta, p0, q0, *pij, *qij;

    // Local: subproblem constant terms, dual gradient, dual hessian
    PetscScalar *b, *grad, *Hess;

    // Global: Old design variables
    Vec xo1, xo2;

    // Math helpers
    PetscErrorCode Factorize(PetscScalar* K, PetscInt nn);
    PetscErrorCode Solve(PetscScalar* K, PetscScalar* x, PetscInt nn);
    PetscScalar    Min(PetscScalar d1, PetscScalar d2);
    PetscScalar    Max(PetscScalar d1, PetscScalar d2);
    PetscInt       Min(PetscInt d1, PetscInt d2);
    PetscInt       Max(PetscInt d1, PetscInt d2);
    PetscScalar    Abs(PetscScalar d1);
};

#endif
