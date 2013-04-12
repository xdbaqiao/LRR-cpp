#include <stdio.h>
#include <string.h>
#include <iostream> 
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

static const double Tol = 1e-4;
static const int MaxIter = 1e6;
static const double Rho = 1.1;
static const long int MaxMu = 1e10;
static double Mu = 1e-6;


class LowRankRepresentation
{
    public:
        LowRankRepresentation(MatrixXd&, MatrixXd&, double){};
        vector<MatrixXd> result(MatrixXd&, MatrixXd&, double);
        MatrixXd solve_l1l2(MatrixXd&, double);
    private:
        VectorXd solve_l2(VectorXd&, double);
};
