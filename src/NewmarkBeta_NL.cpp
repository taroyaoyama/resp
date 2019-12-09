
#include <Rcpp.h>
#include <RcppEigen.h>

//'
//'
//'@export

//[[Rcpp::export]]
Rcpp::List NewmarkBeta_NL(
    Eigen::Map<Eigen::MatrixXd> dat,
    double w, double h, double dt, double alp, double dy_pls,
    double m = 1.0, double beta = 0.25, double gamma = 0.5) {

    int len = dat.rows();

    // params
    double k0, k1, c1;
    k0 = m * pow(w, 2.0);
    k1 = alp * k0;
    c1 = 2.0 * h * w * m;

    // resampling data
    int len_rs; double dt_rs; Eigen::MatrixXd dat_rs;
    len_rs = (len - 1) * 10;
     dt_rs = dt / 10;
    dat_rs = Eigen::MatrixXd::Zero(len_rs, 1);

    for (int i = 0; i < (len - 1); i++) {
        for (int j = 0; j < 10; j++) {
            dat_rs((10 * i + j), 0) = dat(i, 0) + j / 10 * (dat(i + 1, 0) - dat(i, 0));
        }
    }

    // calculation
    Eigen::MatrixXd   Time = Eigen::MatrixXd::Zero(len - 1, 1);
    Eigen::MatrixXd  inAcc = Eigen::MatrixXd::Zero(len - 1, 1);
    Eigen::MatrixXd absAcc = Eigen::MatrixXd::Zero(len - 1, 1);
    Eigen::MatrixXd relAcc = Eigen::MatrixXd::Zero(len - 1, 1);
    Eigen::MatrixXd relVel = Eigen::MatrixXd::Zero(len - 1, 1);
    Eigen::MatrixXd relDis = Eigen::MatrixXd::Zero(len - 1, 1);
    Eigen::MatrixXd rfHist = Eigen::MatrixXd::Zero(len - 1, 1);

    double a1, a2, a3, a4, a5, a6, a7;
    a1 = 1.0 / (beta * dt_rs);
    a2 = 1.0 / (2.0 * beta);
    a3 = gamma / beta;
    a4 = (gamma / (2.0 * beta) - 1.0) * dt_rs;
    a5 = 1.0 / (beta * pow(dt_rs, 2.0));
    a6 = gamma / (beta * dt_rs);
    a7 = (1.0 - gamma / (2.0 * beta)) * dt;

    double k, RF, dis, vel, acc;
    k = k0; RF = 0.0; dis = 0.0; vel = 0.0; acc = -dat_rs(1);

    int yield; double dy_mns;
    yield = 0; dy_mns = -dy_pls;

    int j; j = 1;

    double a8, dy, ddis, dvel, dacc, tac, t;
    t = 0.0;

    for (int i = 1; i < len_rs; i++) {

        a8 = a5 * m + a6 * c1 + k;
        dy = dat_rs(i, 0) - dat_rs(i - 1, 0);

        ddis = (m * (a1 * vel + a2 * acc - dy) + c1 * (a3 * vel + a4 * acc)) / a8;
        dvel = a6 * ddis - a3 * vel + a7 * acc;
        dacc = a5 * ddis - a1 * vel - a2 * acc;

        dis += ddis; vel += dvel; acc += dacc; RF += k * ddis; t += 0;
        tac  = acc + dat_rs(i, 0);

        // bi-linear rule
        if (dis > dy_pls) {

            yield = 1; k = k1;
        }
        if (yield == 1 && ddis < 0.0) {

            yield = 0; k = k0;
            dy_mns += (dis - dy_pls);
            dy_pls  = dis;
        }
        if (dis < dy_mns) {

            yield = -1; k = k1;
        }
        if (yield == -1 && ddis > 0.0) {

            yield = 0; k = k0;
            dy_pls += (dis - dy_mns);
            dy_mns  = dis;
        }

        // data storing
        if (i % 10 == 0) {

            t += dt; Time(j, 0) = t;

             inAcc(j, 0) = dat(j);
            absAcc(j, 0) = tac;
            relAcc(j, 0) = acc;
            relVel(j, 0) = vel;
            relDis(j, 0) = dis;
            rfHist(j, 0) = RF;

            j += 1;
        }
    }

    Rcpp::List L = Rcpp::List::create(

        Rcpp::Named("Time")   = Time,
        Rcpp::Named("inAcc")  = inAcc,
        Rcpp::Named("absAcc") = absAcc,
        Rcpp::Named("relAcc") = relAcc,
        Rcpp::Named("relVel") = relVel,
        Rcpp::Named("relDis") = relDis,
        Rcpp::Named("rfHist") = rfHist

    );

    return L;
}
