#pragma once
#include <RcppEigen.h>

Eigen::MatrixXd chol_Eigen(Eigen::Map<Eigen::MatrixXd> M);
Eigen::MatrixXd inv_Eigen(Eigen::Map<Eigen::MatrixXd> M);
