#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct Tuple{
    double pe, f;
} Tuple;

Tuple compute_row(double r, double rc, double alpha, double m, double n){
    //calculate potential energy and force as a fn of r due to the short range ptl proposed by Wang
    //et al. (2019). sigma=1.
    static double eps = 1.0;

    double ri = 1./r;
    double ri2 = ri*ri;
    double rc2 = rc*rc;
    double rcri2 = rc2*ri2;
    double A = pow(r, -2* 1); 
    double B = pow(rc, 2*m) * pow(r, -2*m) - 1.;

    Tuple row;

    row.pe = alpha * eps * (pow(ri2,2 * m) - 1.) * pow(pow(rcri2,2 * m) - 1.,2*n) + eps; // Potential energy
    row.f = 2.*alpha * eps *m * pow(r,-2*m -1) * pow(B, 2*n -1) * (B + 2*n * pow(rc,2*m)*A); // force

    return row;
}

int main(int argc, char* argv[]){

    if(argc!=7){
        printf("Require 6 args: L (box length), resolution, rc (cutoff),m ,n,  filename. Exiting.\n");
        return 1;
    }

    double L = atof(argv[1]); // box length
    double res = atof(argv[2]); // resolution of potential table
    double rc = atof(argv[3]);
    double m = atof(argv[4]);
    double n = atof(argv[5]);

    int N = L/(2*res) + 1;

    double rmin = rc * pow(3./(1.+2.*rc*rc), 0.5);

    // get filename from argv[4] directly
    FILE* fp = fopen(argv[6], "w");
    fprintf(fp, "#UNITS: real\n\n");
    fprintf(fp, "colloid\n");
    fprintf(fp, "N %i\n\n", N);
    /* alpha is the prefactor that ensures the potential is continuous at rc.
     It is derived by setting the potential energy at r=rc to 0 and solving for alpha.*/
    double alpha = (pow(1 + 2* n, 2* n + 1)*pow(rc, 2*m)) / ((pow(2*n,2*n) * pow(pow(rc,2*m)-1, 2*n+1))) ;

    double r = 0.;
    int i=0;

    while(r<L/2.){
        r += res;
        i++;
        if (r>rmin){
            //write 0 for pe and force for all r beyond rc
            fprintf(fp, "%i %lf 0 0\n", i, r);
        }

        else{
            // calculate potential
            Tuple row = compute_row(r, rc, alpha,m,n);
            fprintf(fp, "%i %lf %lf %lf\n", i, r, row.pe, row.f);
        }
    }

    fclose(fp);
    return 0;
}