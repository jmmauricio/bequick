void f_ini_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[0] = (p[0]*x[1]*y[1] - p[1]*x[0] + y[2])/p[0];
out[1] = (-p[0]*x[0]*y[1] - p[3]*y[1] - p[1]*x[1] + y[3])/p[0];
out[2] = -x[0] + u[1];
out[3] = -x[1] + y[4];
out[4] = (-p[10]*p[7]*p[8]*y[5] - 0.5*p[13]*p[11]*p[12]*pow(x[4], 2)*y[5] - p[7]*p[8]*sin(u[2]) + y[0]/(p[9]*p[6]))/p[8];
out[5] = x[4] - 0.0001*x[5];

}
void g_ini_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[0] = 1.5*p[2]*p[3]*x[1] - y[0];
out[1] = p[2]*y[6] - y[1];
out[2] = -p[5]*x[2] - p[4]*(-x[0] + u[1]) + p[0]*x[1]*y[1] + y[2];
out[3] = -p[5]*x[3] - p[4]*(-x[1] + y[4]) - p[0]*x[0]*y[1] - p[3]*y[1] + y[3];
out[4] = 1.5*p[2]*p[3]*y[4] - u[0];
out[5] = -y[5] - 1 + 2/(1 + exp(-p[14]*x[4]));
out[6] = -y[6] + x[4]/(p[9]*p[6]);

}
void f_run_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[0] = (p[0]*x[1]*y[1] - p[1]*x[0] + y[2])/p[0];
out[1] = (-p[0]*x[0]*y[1] - p[3]*y[1] - p[1]*x[1] + y[3])/p[0];
out[2] = -x[0] + u[1];
out[3] = -x[1] + y[4];
out[4] = (-p[10]*p[7]*p[8]*y[5] - 0.5*p[13]*p[11]*p[12]*pow(x[4], 2)*y[5] - p[7]*p[8]*sin(u[2]) + y[0]/(p[9]*p[6]))/p[8];
out[5] = x[4] - 0.0001*x[5];

}
void g_run_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[0] = 1.5*p[2]*p[3]*x[1] - y[0];
out[1] = p[2]*y[6] - y[1];
out[2] = -p[5]*x[2] - p[4]*(-x[0] + u[1]) + p[0]*x[1]*y[1] + y[2];
out[3] = -p[5]*x[3] - p[4]*(-x[1] + y[4]) - p[0]*x[0]*y[1] - p[3]*y[1] + y[3];
out[4] = 1.5*p[2]*p[3]*y[4] - u[0];
out[5] = -y[5] - 1 + 2/(1 + exp(-p[14]*x[4]));
out[6] = -y[6] + x[4]/(p[9]*p[6]);

}
void h_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[0] = p[5]*x[2] + p[4]*(-x[0] + u[1]);
out[1] = p[5]*x[3] + p[4]*(-x[1] + y[4]);
out[2] = u[0];
out[3] = sqrt(pow(y[2], 2) + pow(y[3], 2));
out[4] = sqrt(pow(x[0], 2) + pow(x[1], 2));
out[5] = y[6];
out[6] = y[0]/(p[9]*p[6]);
out[7] = 0.5*p[13]*p[11]*p[12]*pow(x[4], 2)*y[5];
out[8] = p[10]*p[7]*p[8]*y[5];
out[9] = y[0];
out[10] = y[6]*y[0];
out[11] = 3.6000000000000001*x[4];

}
void de_jac_ini_xy_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[1] = y[1];
out[7] = x[1];
out[13] = -y[1];
out[20] = (-p[0]*x[0] - p[3])/p[0];
out[56] = -1.0*p[13]*p[11]*p[12]*x[4]*y[5]/p[8];
out[63] = (-p[10]*p[7]*p[8] - 0.5*p[13]*p[11]*p[12]*pow(x[4], 2))/p[8];
out[105] = p[0]*y[1];
out[111] = p[0]*x[1];
out[117] = -p[0]*y[1];
out[124] = -p[0]*x[0] - p[3];
out[147] = 2*p[14]*exp(-p[14]*x[4])/pow(1 + exp(-p[14]*x[4]), 2);

}

void de_jac_ini_up_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[0] = -p[1]/p[0];
out[8] = 1.0/p[0];
out[14] = -p[1]/p[0];
out[22] = 1.0/p[0];
out[58] = 1/(p[9]*p[8]*p[6]);
out[79] = 1.5*p[2]*p[3];
out[103] = p[2];
out[104] = p[4];
out[106] = -p[5];
out[118] = p[4];
out[120] = -p[5];
out[127] = -p[4];
out[140] = 1.5*p[2]*p[3];
out[160] = 1/(p[9]*p[6]);

}

void de_jac_ini_num_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[26] = -1;
out[40] = -1;
out[49] = 1;
out[69] = 1;
out[70] = -0.0001;
out[84] = -1;
out[98] = -1;
out[112] = 1;
out[126] = 1;
out[154] = -1;
out[168] = -1;

}

void sp_jac_ini_xy_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[1] = y[1];
out[2] = x[1];
out[4] = -y[1];
out[6] = (-p[0]*x[0] - p[3])/p[0];
out[11] = -1.0*p[13]*p[11]*p[12]*x[4]*y[5]/p[8];
out[13] = (-p[10]*p[7]*p[8] - 0.5*p[13]*p[11]*p[12]*pow(x[4], 2))/p[8];
out[21] = p[0]*y[1];
out[23] = p[0]*x[1];
out[25] = -p[0]*y[1];
out[28] = -p[0]*x[0] - p[3];
out[32] = 2*p[14]*exp(-p[14]*x[4])/pow(1 + exp(-p[14]*x[4]), 2);

}

void sp_jac_ini_up_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[0] = -p[1]/p[0];
out[3] = 1.0/p[0];
out[5] = -p[1]/p[0];
out[7] = 1.0/p[0];
out[12] = 1/(p[9]*p[8]*p[6]);
out[16] = 1.5*p[2]*p[3];
out[19] = p[2];
out[20] = p[4];
out[22] = -p[5];
out[26] = p[4];
out[27] = -p[5];
out[30] = -p[4];
out[31] = 1.5*p[2]*p[3];
out[34] = 1/(p[9]*p[6]);

}

void sp_jac_ini_num_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[8] = -1;
out[9] = -1;
out[10] = 1;
out[14] = 1;
out[15] = -0.0001;
out[17] = -1;
out[18] = -1;
out[24] = 1;
out[29] = 1;
out[33] = -1;
out[35] = -1;

}

void de_jac_run_xy_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[1] = y[1];
out[7] = x[1];
out[13] = -y[1];
out[20] = (-p[0]*x[0] - p[3])/p[0];
out[56] = -1.0*p[13]*p[11]*p[12]*x[4]*y[5]/p[8];
out[63] = (-p[10]*p[7]*p[8] - 0.5*p[13]*p[11]*p[12]*pow(x[4], 2))/p[8];
out[105] = p[0]*y[1];
out[111] = p[0]*x[1];
out[117] = -p[0]*y[1];
out[124] = -p[0]*x[0] - p[3];
out[147] = 2*p[14]*exp(-p[14]*x[4])/pow(1 + exp(-p[14]*x[4]), 2);

}

void de_jac_run_up_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[0] = -p[1]/p[0];
out[8] = 1.0/p[0];
out[14] = -p[1]/p[0];
out[22] = 1.0/p[0];
out[58] = 1/(p[9]*p[8]*p[6]);
out[79] = 1.5*p[2]*p[3];
out[103] = p[2];
out[104] = p[4];
out[106] = -p[5];
out[118] = p[4];
out[120] = -p[5];
out[127] = -p[4];
out[140] = 1.5*p[2]*p[3];
out[160] = 1/(p[9]*p[6]);

}

void de_jac_run_num_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[26] = -1;
out[40] = -1;
out[49] = 1;
out[69] = 1;
out[70] = -0.0001;
out[84] = -1;
out[98] = -1;
out[112] = 1;
out[126] = 1;
out[154] = -1;
out[168] = -1;

}

void sp_jac_run_xy_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[1] = y[1];
out[2] = x[1];
out[4] = -y[1];
out[6] = (-p[0]*x[0] - p[3])/p[0];
out[11] = -1.0*p[13]*p[11]*p[12]*x[4]*y[5]/p[8];
out[13] = (-p[10]*p[7]*p[8] - 0.5*p[13]*p[11]*p[12]*pow(x[4], 2))/p[8];
out[21] = p[0]*y[1];
out[23] = p[0]*x[1];
out[25] = -p[0]*y[1];
out[28] = -p[0]*x[0] - p[3];
out[32] = 2*p[14]*exp(-p[14]*x[4])/pow(1 + exp(-p[14]*x[4]), 2);

}

void sp_jac_run_up_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[0] = -p[1]/p[0];
out[3] = 1.0/p[0];
out[5] = -p[1]/p[0];
out[7] = 1.0/p[0];
out[12] = 1/(p[9]*p[8]*p[6]);
out[16] = 1.5*p[2]*p[3];
out[19] = p[2];
out[20] = p[4];
out[22] = -p[5];
out[26] = p[4];
out[27] = -p[5];
out[30] = -p[4];
out[31] = 1.5*p[2]*p[3];
out[34] = 1/(p[9]*p[6]);

}

void sp_jac_run_num_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[8] = -1;
out[9] = -1;
out[10] = 1;
out[14] = 1;
out[15] = -0.0001;
out[17] = -1;
out[18] = -1;
out[24] = 1;
out[29] = 1;
out[33] = -1;
out[35] = -1;

}

void de_jac_trap_xy_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[1] = -0.5*Dt*y[1];
out[7] = -0.5*Dt*x[1];
out[13] = 0.5*Dt*y[1];
out[20] = -0.5*Dt*(-p[0]*x[0] - p[3])/p[0];
out[56] = 0.5*p[13]*Dt*p[11]*p[12]*x[4]*y[5]/p[8] + 1;
out[63] = -0.5*Dt*(-p[10]*p[7]*p[8] - 0.5*p[13]*p[11]*p[12]*pow(x[4], 2))/p[8];
out[105] = p[0]*y[1];
out[111] = p[0]*x[1];
out[117] = -p[0]*y[1];
out[124] = -p[0]*x[0] - p[3];
out[147] = 2*p[14]*exp(-p[14]*x[4])/pow(1 + exp(-p[14]*x[4]), 2);

}

void de_jac_trap_up_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[0] = 0.5*Dt*p[1]/p[0] + 1;
out[8] = -0.5*Dt/p[0];
out[14] = 0.5*Dt*p[1]/p[0] + 1;
out[22] = -0.5*Dt/p[0];
out[26] = 0.5*Dt;
out[40] = 0.5*Dt;
out[49] = -0.5*Dt;
out[58] = -0.5*Dt/(p[9]*p[8]*p[6]);
out[69] = -0.5*Dt;
out[70] = 5.0000000000000002e-5*Dt + 1;
out[79] = 1.5*p[2]*p[3];
out[103] = p[2];
out[104] = p[4];
out[106] = -p[5];
out[118] = p[4];
out[120] = -p[5];
out[127] = -p[4];
out[140] = 1.5*p[2]*p[3];
out[160] = 1/(p[9]*p[6]);

}

void de_jac_trap_num_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[28] = 1;
out[42] = 1;
out[84] = -1;
out[98] = -1;
out[112] = 1;
out[126] = 1;
out[154] = -1;
out[168] = -1;

}

void sp_jac_trap_xy_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[1] = -0.5*Dt*y[1];
out[2] = -0.5*Dt*x[1];
out[4] = 0.5*Dt*y[1];
out[6] = -0.5*Dt*(-p[0]*x[0] - p[3])/p[0];
out[13] = 0.5*p[13]*Dt*p[11]*p[12]*x[4]*y[5]/p[8] + 1;
out[15] = -0.5*Dt*(-p[10]*p[7]*p[8] - 0.5*p[13]*p[11]*p[12]*pow(x[4], 2))/p[8];
out[23] = p[0]*y[1];
out[25] = p[0]*x[1];
out[27] = -p[0]*y[1];
out[30] = -p[0]*x[0] - p[3];
out[34] = 2*p[14]*exp(-p[14]*x[4])/pow(1 + exp(-p[14]*x[4]), 2);

}

void sp_jac_trap_up_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[0] = 0.5*Dt*p[1]/p[0] + 1;
out[3] = -0.5*Dt/p[0];
out[5] = 0.5*Dt*p[1]/p[0] + 1;
out[7] = -0.5*Dt/p[0];
out[8] = 0.5*Dt;
out[10] = 0.5*Dt;
out[12] = -0.5*Dt;
out[14] = -0.5*Dt/(p[9]*p[8]*p[6]);
out[16] = -0.5*Dt;
out[17] = 5.0000000000000002e-5*Dt + 1;
out[18] = 1.5*p[2]*p[3];
out[21] = p[2];
out[22] = p[4];
out[24] = -p[5];
out[28] = p[4];
out[29] = -p[5];
out[32] = -p[4];
out[33] = 1.5*p[2]*p[3];
out[36] = 1/(p[9]*p[6]);

}

void sp_jac_trap_num_eval(double *out,double *x,double *y,double *u,double *p,double Dt){

out[9] = 1;
out[11] = 1;
out[19] = -1;
out[20] = -1;
out[26] = 1;
out[31] = 1;
out[35] = -1;
out[37] = -1;

}

