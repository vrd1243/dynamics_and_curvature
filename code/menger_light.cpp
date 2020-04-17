#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

int main(int argc, char **argv) {
    
    int i = 0;
     
    int tau_avg = 1;
    int n_neighbors = 0; 
    int length;// = atoi(argv[3]);
    int embed_dim;// = atoi(argv[4]);
    int tau_max;// = atoi(argv[5]);
    int tau_start = 1;
    ifstream f_in(argv[1]);
    ofstream f_out(argv[2]);
    ifstream f_params(argv[3]);

    f_params >> length;
    f_params >> embed_dim;
    f_params >> tau_max;
    f_params >> tau_avg;    
    f_params >> n_neighbors; 
    
    cout << length << endl;
    cout << embed_dim << endl;
    cout << tau_max << endl;
    cout << tau_avg << endl;
    cout << n_neighbors << endl;

    double series[length];
    double menger_sum[tau_max];
    double menger_variance[tau_max];
    
    while(1) {
	f_in >> series[i++];
	if (f_in.eof()) break;
    }
    double m = 0;
    for (int tau = tau_start; tau < tau_max; tau++) {
      int menger_length = length - 2*tau_avg - 2*n_neighbors - (embed_dim - 1)*tau;
      
      menger_sum[tau] = 0;
      double menger_square_sum = 0;

      for (int t = tau_avg + n_neighbors; t < length - (embed_dim - 1)*tau - tau_avg - n_neighbors; t++) {

          m = 0; 

	  for (int k = -n_neighbors; k < n_neighbors + 1; k++) {
              
	      double a=0, b=0, c=0; 
              double embed_point_A, embed_point_B, embed_point_C;
              
	      for (int d = 0; d < embed_dim; d++) {
                  embed_point_A = series[t - tau_avg + k + tau*d];	
                  embed_point_B = series[t + k + tau*d];	
                  embed_point_C = series[t + tau_avg + k + tau*d];	

		  a += (embed_point_B  - embed_point_C)*(embed_point_B  - embed_point_C);
		  b += (embed_point_A  - embed_point_C)*(embed_point_A  - embed_point_C);
		  c += (embed_point_B  - embed_point_A)*(embed_point_B  - embed_point_A);
	      }
	      a = sqrt(a);
	      b = sqrt(b);
	      c = sqrt(c);

	      /*double x_first = series[t-tau_avg+k], y_first = series[t+tau-tau_avg+k], z_first = series[t + 2*tau - tau_avg + k];
	      double x_second = series[t+k], y_second = series[t+tau+k], z_second = series[t + 2*tau + k];
	      double x_third = series[t+tau_avg+k], y_third = series[t+tau+tau_avg+k], z_third = series[t + 2*tau + tau_avg + k];
	      
	      double a = sqrt((x_first - x_second)*(x_first - x_second) + 
			     (y_first - y_second)*(y_first - y_second) +
			     (z_first - z_second)*(z_first - z_second));

	      double b = sqrt((x_third - x_second)*(x_third - x_second) + 
			     (y_third - y_second)*(y_third - y_second) +
			     (z_third - z_second)*(z_third - z_second));

	      double c = sqrt((x_first - x_third)*(x_first - x_third) + 
			     (y_first - y_third)*(y_first - y_third) +
			     (z_first - z_third)*(z_first - z_third));
              */
	      double s = (a + b + c)/2;
	      //if (isnan(sqrt(s*(s-a)*(s-b)*(s-c))))
	      //	cout << a << " " << b << " " << c << " " << s << " " << s*(s-a)*(s-b)*(s-c) << endl;
	      
	      if (((a*b*c) != 0) && (s*(s-a)*(s-b)*(s-c) >= 0))
		  m += 4*sqrt(s*(s-a)*(s-b)*(s-c)) / (a*b*c);
	     
	  }
	  menger_sum[tau] += m / (2*n_neighbors + 1);
	  menger_square_sum += (m / (2*n_neighbors + 1))*(m / (2*n_neighbors + 1));
      }
      
      double menger_mean = menger_sum[tau] / (menger_length);
      menger_sum[tau] = menger_mean;

      menger_variance[tau] = menger_square_sum / menger_length - menger_mean * menger_mean;	

    }

    for(int tau = tau_start; tau < tau_max; tau++) {
        f_out << menger_sum[tau] << " " << menger_variance[tau] << endl;  
    }
       
}
