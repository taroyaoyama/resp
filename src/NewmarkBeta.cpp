#include <Rcpp.h>
#include <RcppEigen.h>

//'
//'
//'@export

//[[Rcpp::export]]
Rcpp::List NewmarkBeta(
    Eigen::Map<Eigen::MatrixXd> dat, Eigen::Map<Eigen::MatrixXd> M,
    Eigen::Map<Eigen::MatrixXd> C, Eigen::Map<Eigen::MatrixXd> K,
    double dt, double beta = 1/4, double gamma = 1/2) {

    int len = dat.rows();
    int node = M.rows();

    Eigen::MatrixXd   Time = Eigen::MatrixXd::Zero(len,    1);
    Eigen::MatrixXd  inAcc = Eigen::MatrixXd::Zero(len,    1);
    Eigen::MatrixXd absAcc = Eigen::MatrixXd::Zero(len, node);
    Eigen::MatrixXd relAcc = Eigen::MatrixXd::Zero(len, node);
    Eigen::MatrixXd relVel = Eigen::MatrixXd::Zero(len, node);
    Eigen::MatrixXd relDis = Eigen::MatrixXd::Zero(len, node);

    inAcc(0, 0) = dat(0, 0);

    Eigen::VectorXd acc = Eigen::VectorXd::Constant(node, -dat(0, 0));
    Eigen::VectorXd vel = Eigen::VectorXd::Zero(node);
    Eigen::VectorXd dis = Eigen::VectorXd::Zero(node);
    Eigen::VectorXd tac = Eigen::VectorXd::Zero(node);

    double a1, a2, a3, a4, a5, a6, a7;
    Eigen::MatrixXd A8;

    a1 = 1 / (beta * dt);
    a2 = 1 / (beta * 2.0);
    a3 = gamma / beta;
    a4 = (gamma / (2 * beta) - 1) * dt;
    a5 = 1 / (beta * pow(dt, 2.0));
    a6 = gamma / beta / dt;
    a7 = (1 - gamma / 2 / beta) * dt;
    A8 = a5 * M.array() + a6 * C.array() + K.array();

    double t = 0.0;
    Eigen::VectorXd eye = Eigen::VectorXd::Ones(node);

    for (int i = 1; i < len; i++) {

        double dy = dat(i, 0) - dat(i-1, 0);
        Eigen::MatrixXd ddis = A8.inverse() * (
            M * (a1 * vel + a2 * acc - dy * eye) +
            C * (a3 * vel + a4 * acc) );
        Eigen::MatrixXd dvel = a6 * ddis - a3 * vel + a7 * acc;
        Eigen::MatrixXd dacc = a5 * ddis - a1 * vel - a2 * acc;

        t   += dt;
        dis += ddis;
        vel += dvel;
        acc += dacc;
        tac  = acc + Eigen::VectorXd::Constant(node, dat(i, 0));

          Time.row(i) = Eigen::VectorXd::Constant(1, t);
         inAcc.row(i) = Eigen::VectorXd::Constant(1, dat(i, 0));
        absAcc.row(i) = tac.transpose();
        relAcc.row(i) = acc.transpose();
        relVel.row(i) = vel.transpose();
        relDis.row(i) = dis.transpose();

    }

    Rcpp::List L = Rcpp::List::create(

        Rcpp::Named("Time")   = Time,
        Rcpp::Named("inAcc")  = inAcc,
        Rcpp::Named("absAcc") = absAcc,
        Rcpp::Named("relAcc") = relAcc,
        Rcpp::Named("relVel") = relVel,
        Rcpp::Named("relDis") = relDis

    );

    return L;
}
