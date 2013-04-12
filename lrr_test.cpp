#include "mex.h"
#include "lrr_test.h"
#define Z_OUT plhs[0]
#define E_OUT plhs[1]
#define X_IN  prhs[0]
#define D_IN  prhs[1]
#define lambda_IN  prhs[2]

 /* Always include this */
void mexFunction(int nlhs, mxArray *plhs[],                       /* Output variables */
                          int nrhs, const mxArray *prhs[])             /* Input variables */
{

if(nrhs < 1 || nrhs > 3)                /* Check the number of arguments */
     mexErrMsgTxt("Wrong number of input arguments.");
else if(nlhs > 2)
     mexErrMsgTxt("Too many output arguments.");
/* if (nrhs == 2)                     
     lambda = 0.02;
else
     lambda = prhs[2];
*/

double *X, *D,*Z,*E,lambda;
int XM, XN,DM,DN, ZM,ZN,m, n;

lambda = mxGetScalar(lambda_IN);                         /* Get lambda */

X = mxGetPr(X_IN);                                                /* Get the pointer to  the data of X */
XM = mxGetM(X_IN);                                              /* Get the dimensions  of X */
XN = mxGetN(X_IN);

D = mxGetPr(D_IN);                                               
DM = mxGetM(D_IN);                                            
DN = mxGetN(D_IN);

ZM = DN;
ZN =  XN;
m=0,n=0;
//X[m + M*n]

MatrixXd x(XM, XN) ,d(DM,DN),z(ZM,ZN);

vector<MatrixXd> ze;

Z_OUT = mxCreateDoubleMatrix(DN, XN, mxREAL);
Z = mxGetPr(Z_OUT);                                             /* Get the pointer to  the data of Z */
E_OUT = mxCreateDoubleMatrix(DN, XN, mxREAL);
E = mxGetPr(E_OUT);

//cout << "ZM" << ZM << "ZN" << ZN << endl;;

for (n=0;n<XN;n++)
    for(m=0;m<XM;m++)
        x(m,n) = X[m+ZM*n];


for (n=0;n<DN;n++)
    for(m=0;m<DM;m++)
        d(m,n) = D[m+ZM*n];


LowRankRepresentation lrr(x, d, lambda);
ze = lrr.result(x, d, lambda);



for (n=0;n<ZN;n++)
    for(m=0;m<ZM;m++)
        Z[m+ZM*n] = ze[0](m,n);

for (n=0;n<ZN;n++)
    for(m=0;m<ZM;m++)
        E[m+ZM*n] = ze[1](m,n);

//Z = ZE[0]; 
//E = ZE[1];
        
        
//mexPrintf("Hello, world!\n");                                   /* Do something interesting */
return;
}

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
    MatrixXd atx, J, Z, Y1, Y2, inv_a, tmp, xmaz, leq1, leq2;
    VectorXd sigma_tmp;
    int d, n, m, svp, M, N, i, j;
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

//    FullPivLU<MatrixXd> lu_decomp(Z);
//    cout << "initial, rank=" << lu_decomp.rank() << endl;

    while(iter < MaxIter)
    {
        iter += 1;
        tmp = Z + Y2/Mu;
        M = tmp.rows();
        N = tmp.cols();
        double d_tmp[M*N];
        for(i=0; i<M; i++)
            for(j=0; j<N; j++)
                d_tmp[i + j*N] = tmp(i, j);

        mxArray *U_pro, *sigma_pro, *V_pro, *B;
        mxArray *ppLhs[3];
        B = mxCreateDoubleMatrix(M, N, mxREAL); 
        //ppRhs[0] = B;
        //ppRhs[1] = "econ";
        memcpy(mxGetPr(B), d_tmp, sizeof(double) * M * N);

        mexCallMATLAB(3, ppLhs, 1, &B, "svd");
        U_pro = ppLhs[0]; 
        sigma_pro = ppLhs[1]; 
        V_pro = ppLhs[2]; 

        M = mxGetM(U_pro);
        N = mxGetN(U_pro);
        MatrixXd U(M,N);
        for(i=0;i<M;i++)
            for(j=0; j<N; j++)
                U(i, j) = mxGetPr(U_pro)[i + j*N];

        VectorXd sigma(N);
        for(i=0; i<N; i++)
            sigma(i) = mxGetPr(sigma_pro)[i + i*N];

        N = mxGetN(V_pro);
        MatrixXd V(M,N);
        for(i=0; i<N; i++)
            for(j=0; j<M; j++)
                V(i, j) = mxGetPr(V_pro)[i + j*M];


        mxDestroyArray(U_pro);
        mxDestroyArray(sigma_pro); 
        mxDestroyArray(V_pro); 
        mxDestroyArray(B); 

//        JacobiSVD<MatrixXd> svd(tmp, ComputeThinU | ComputeThinV);
//        sigma = svd.singularValues();
//        //sigma = MatrixXd(sigma.asDiagonal());
//        U = svd.matrixU();
//        V = svd.matrixV();
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


        
//        if (iter==1 || (iter % 50)==0 || stopC<Tol)
//        {
//            FullPivLU<MatrixXd> lu_decomp(Z);
            //printf("iter %d, mu=%2.1e, rank=%d, stopALM=%2.3e\n", iter, Mu, lu_decomp.rank(), stopC);
 //           printf("iter %d, mu=%2.1e, stopALM=%2.3e\n", iter, Mu, stopC);
//           cout<< endl;
//        }
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
