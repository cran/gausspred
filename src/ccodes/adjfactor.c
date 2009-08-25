double sum_double(int len,double A[len])
{  int i; double out=0;
   for(i = 0; i < len; i++) out += A[i];
   return(out);  
} 

void comp_adjfactor
       (int NF[1], int NL[1], double Q_F[NF[0]], double lambda[NL[0]], 
        double adjfactor[1]) 
{  
     double w_poi[NF[0]][NL[0]], adjfactor_F[NF[0]] ;
     int k,l, max_K;

     for(l = 0; l < NL[0]; l++) {
         w_poi[0][l] = exp( - lambda[l] );
         for(k = 1; k < NF[0]; k++) {
             w_poi[k][l] = w_poi[k-1][l] * lambda[l] / k; 
         }
     }
     //Rprintf("lim K = %d\n", lim_K);
     
     for(k = 0; k < NF[0]; k++) 
        adjfactor_F[k] = Q_F[k] * sum_double( NL[0], &w_poi[k][0] ) / NL[0];
     //lsw = log_sum_exp(K,ladjfactor);
     
     adjfactor[0] = sum_double(NF[0], &adjfactor_F[0]);
}

