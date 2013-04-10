#include "LowRankRepresentation.h"

int main(int argc, char*argv[])
{
    MatrixXd a(3, 4) ,b(3,4);
    vector<MatrixXd> ZE;
    a << 1,2,3,4,5,6,7,8,9,10,11,12;
    b << 13,14,15,16,17,18,19,20,21,22,23,24;

    double lambda=0.02;
    LowRankRepresentation lrr(a, b, lambda);
    ZE = lrr.result(a, b, lambda);
    cout << "Z = " << endl << ZE[0] <<endl;
    cout << "E = " << endl << ZE[1] <<endl;
    return 0;
}

vector<MatrixXd> LowRankRepresentation::result(MatrixXd& X, MatrixXd& A, double lambda)
{
    vector<MatrixXd> ZE;
    MatrixXd atx, J, Z, Y1, Y2, inv_a, tmp, U, V, xmaz, leq1, leq2;
    VectorXd sigma, sigma_tmp;
    int d, n, m, svp;
    double stopC, stopC_tmp1, stopC_tmp2;

    d = X.rows();
    n = X.cols();
    m = A.cols();

    atx = A.transpose() * X;
    inv_a = (A.transpose() *A + MatrixXd::Identity(m, m)).inverse();
    J = MatrixXd::Zero(m, n);
    Z = MatrixXd::Zero(m, n);
    Y1 = MatrixXd::Zero(d, n);
    Y2 = MatrixXd::Zero(m, n);
    MatrixXd E;
    E = MatrixXd::Zero(d, n);
//   SparseMatrix<double> E(d, n);
    int iter = 0;

    FullPivLU<MatrixXd> lu_decomp(Z);
    cout << "initial, rank=" << lu_decomp.rank() << endl;

    while(iter < MaxIter)
    {
        iter += 1;
        tmp = Z + Y2/Mu;
        JacobiSVD<MatrixXd> svd(tmp, ComputeThinU | ComputeThinV);
        sigma = svd.singularValues();
        //sigma = MatrixXd(sigma.asDiagonal());
        U = svd.matrixU();
        V = svd.matrixV();
        svp = (sigma.array() > 1/Mu).count();
        if(svp >= 1)
        {
            sigma_tmp = sigma.segment(0,svp)-(1/Mu) * VectorXd::Ones(svp) ;
            sigma = sigma_tmp;
        }
        else
        {
            svp = 1;
            sigma = VectorXd::Zero(1);
        }
        J = U.leftCols(svp) * MatrixXd(sigma.asDiagonal()) * V.leftCols(svp).transpose();
        Z = inv_a * (atx - A.transpose() * E + J + (A.transpose() * Y1 - Y2)/Mu);

        xmaz = X - A * Z;
        tmp = xmaz + Y1 / Mu;
        E = solve_l1l2(tmp, lambda/Mu);

        leq1 = xmaz - E;
        leq2 = Z - J;
        stopC_tmp1 = leq1.array().abs().maxCoeff();
        stopC_tmp2 = leq2.array().abs().maxCoeff();
        stopC = stopC_tmp1 < stopC_tmp2 ? stopC_tmp2:stopC_tmp1;

        if (iter==1 || (iter % 50)==0 || stopC<Tol)
        {
            FullPivLU<MatrixXd> lu_decomp(Z);
            //printf("iter %d, mu=%2.1e, rank=%d, stopALM=%2.3e\n", iter, Mu, lu_decomp.rank(), stopC);
            printf("iter %d, mu=%2.1e, stopALM=%2.3e\n", iter, Mu, stopC);
        }
        if (stopC<Tol)
        {
            cout << "LRR done." << endl;
            break;
        }
        else
        {
            Y1 = Y1 + Mu*leq1;
            Y2 = Y2 + Mu*leq2;
            Mu = MaxMu > Mu*Rho ? Mu*Rho:MaxMu;
        }
    }
    ZE.push_back(Z);
    ZE.push_back(E);
    return ZE;
}

MatrixXd LowRankRepresentation::solve_l1l2(MatrixXd& W, double lambda)
{
    MatrixXd E;
    VectorXd temp;
    int n, i;

    n = W.cols();
    E = W;

    for(i=0; i<n; i++)
    {
       temp = W.col(i);
       E.col(i) = solve_l2(temp, lambda);
    }
    return E;
}

VectorXd LowRankRepresentation::solve_l2(VectorXd& w, double lambda)
{
    double nw;
    VectorXd x;

    nw = w.norm();
    if(nw > lambda)
        x = (nw - lambda) * w / nw;
    else
        x = VectorXd::Zero(w.size(),1);
    return x;
}
