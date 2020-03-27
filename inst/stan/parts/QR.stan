  matrix[n, dx] Q_ast;
  matrix[dx, dx] R_ast;
  matrix[dx, dx] R_inverse;
  if (dx) {
  Q_ast = qr_Q(x)[, 1:dx] * sqrt(n - 1);
  R_ast = qr_R(x)[1:dx, ] / sqrt(n - 1);
  R_inverse = inverse(R_ast);
  }

