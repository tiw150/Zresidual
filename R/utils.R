if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Zresidual.coxph.frailty",
    "Zresidual.coxph",
    "Zresidual.survreg",

    "fit_survreg",
    "Ccoxcount1",
    "Ccoxcount2",
    "j",
    "is.outlier"
  ))
}
