#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"bandpass.h"
#include<assert.h>

#define pi 3.14159265358979

double* generatefunction(double minx, double maxx, int numpoints, double* values, int numvalues){
//Minf and maxf denote min and max freqs. Minx and maxx denote min and max limits of X.
//Recall: KERNEL(t) = lowf*sinc(lowf*t) - highf*sinc(highf*t);
	double middle = (maxx-minx)/2;
	double delta = (maxx - minx)/numpoints;	


	double* popcorn = malloc(numpoints * sizeof(double));
	assert(popcorn);
	int i = 0;
	double curr_x = 0-middle;
	for (popcorn; i < numpoints;)
	{
		int j;
		*(popcorn + i) = 0;
		for (j = 0; j < numvalues; j ++){
			*(popcorn + i) += sin(values[j] *pi* curr_x);
		}
		
		curr_x += delta;
		i ++;

	}

	return popcorn;
}


double* dirac(double minx, double maxx, int numpoints, double* values, int numvalues){
//Minf and maxf denote min and max freqs. Minx and maxx denote min and max limits of X.
//Recall: KERNEL(t) = lowf*sinc(lowf*t) - highf*sinc(highf*t);
	double middle = (maxx-minx)/2;
	double delta = (maxx - minx)/numpoints;	


	double* popcorn = malloc(numpoints * sizeof(double));
	assert(popcorn);
	int i = 0;
	double curr_x = 0-middle;
	for (popcorn; i < numpoints;)
	{
		int j;
		*(popcorn + i) = 0;
		for (j = 0; j < numvalues; j ++){

			if (i == (numpoints / 2)){
			*(popcorn + i) = 1;
			//printf("%f\n", curr_x);
			}
			else{
			*(popcorn + i) = 0;}
		}
		
		curr_x += delta;
		i ++;

	}

	return popcorn;
}

double sinc (double input)
{
	double pie = pi;
	//double pie = 1;
	double a = sin((pie)*input)/(input*(pie));
	if (input == 0){ return 1;}

	return a;
	
	
}

double* generatekernel(double minf, double maxf, double minx, double maxx, int numpoints){
//Minf and maxf denote min and max freqs. Minx and maxx denote min and max limits of X.
//Recall: KERNEL(t) = lowf*sinc(lowf*t) - highf*sinc(highf*t);
	double middle = (maxx-minx)/2;
	double delta = (maxx - minx)/numpoints;	


	double* popcorn = malloc(numpoints * sizeof(double));
	assert(popcorn);
	int i = 0;
	double curr_x = 0-middle;
	for (popcorn; i < numpoints;)
	{
		
		*(popcorn + i) = (maxf * sinc(maxf * curr_x)) - (minf * sinc(minf * curr_x));
		curr_x += delta;
		i ++;

	}

	return popcorn;
}

void convolve(double* input, double* kernel, double* output, int inputsize, int kernelsize, int outputsize)
//Convolve two arrays.
//Recall: (f*g)(n) = sum[f(m)g(n-m)] for all m.
/*
Note that size of convo product = INPUTSIZE + OUTPUTSIZE - 1.
*/
{
	int i = 0;
	for (i; i < outputsize; i++){
		double convvalue = 0;
		int inputstart = i - (kernelsize/2);
		int j = 0;
		for(j; j < kernelsize; j++)
			{if (inputstart + j >= 0 && inputstart + j < inputsize){
				convvalue += input[inputstart + j]*kernel[j];
			}
		}
		
		output[i] = convvalue;
	}

}

void convolve_complete(double* input, double* kernel, double* output, int inputsize, int kernelsize, int outputsize)
//Convolve two arrays.
//Recall: (f*g)(n) = sum[f(m)g(n-m)] for all m.
/*
Note that size of convo product = INPUTSIZE + OUTPUTSIZE - 1.
*/
{
	int t = 0; //t is the index of the product
	int tau; //obtain each t by integrating f(tau)g(t-tau). F is the kernel, G is input.
	double value; //value is the value of prod at t
	//The output array's size shall be INPUTSIZE + OUTPUTSIZE - 1.
	for (t; t < outputsize; t ++){
		tau = 0; //obtain each t by integrating f(tau)g(t-tau). F is the kernel, G is input.
		//Find the value of PROD at each t.
		value = 0;
		for (tau; tau < inputsize; tau++){
			/*In an ideal condition, we would integrate
			from plus infinity to minus infinity.
			In our case, though, we can only integrate
			when t - tau is valid. That is, t - tau
			lies between 0 and inputsize.*/
			if (t - tau >= 0 && t - tau < inputsize){
				value += (kernel[tau]*input[(t - tau)]);
			}
		}
		output[t] = value;
	}
	
}

int main(int argc, char** argv){
	char* warning = "\nUsage: bandpass <numpoints> <minfreq> <maxfreq> <filename> <input_filename> <kernel_filename>.\n";

	if (argc != 7){
		printf("%s", warning);
		return 0;
	}
	else{
		//A sample function. Sample rate is 1000hz for this function.
		int functnum = 10000;
		double stuff[12] = {
					0.1,
					1.0, 
					10.0, 
					5.0, 
					100.0
					};
		double* function = generatefunction(-10.00, 10.00, functnum, stuff, 5);		
		//double* function = dirac(-10, 10, functnum, stuff, 12);		
		
		
		
		//Make convo kernel. (WARNING: User should keep both min and max frequencies below the Nyquist frequency, or weird things happen.)
		//Rule of thumb: xmin and xmax are equal to that of function, and sample count is equal to sample.
		int count = atoi(argv[1]);
		double minfreq = atof(argv[2]);
		double maxfreq = atof(argv[3]);

		double* kernel = generatekernel(minfreq, maxfreq, -10.00, 10.00, count);
		
		

		//Convolve. Recall convo result size is INPUTSIZE + OUTPUTSIZE - 1.
		double* convo_result = malloc((count + functnum - 1)*sizeof(double));
		convolve(function, kernel, convo_result, functnum, count, (count + functnum - 1));
		
		int outputsize = (functnum + count) - 1;
		//double* convo_result = malloc(outputsize * sizeof(double));
		//convolve(function, kernel, convo_result, functnum, count, outputsize);



		//Export the results to a CSV.
		char* filename = argv[4];
		char* filename_kernel = argv[5];
		char* filename_input = argv[6];
		int i = 0;
		
/*
		FILE *f = fopen(filename, "w");
		for (
			i; 
			//i < outputsize; 
			i < count;
			i ++
		){
		fprintf(f, "%f\n", convo_result[i]);}
		fclose(f);
*/
		int startpt = ((count + functnum - 1)/2) - (functnum/2);
		FILE *f = fopen(filename, "w");
		//for (i; i < count; i ++){
		for (i; i < count; i ++){
		fprintf(f, "%f\n", convo_result[i]);}
		fclose(f);
		
		i = 0;
		FILE *g = fopen(filename_kernel, "w");
		for (i; i < count; i ++){
		fprintf(g, "%f\n", kernel[i]);}
		fclose(g);

		i = 0;
		FILE *h = fopen(filename_input, "w");
		for (i; i < functnum; i ++){
		fprintf(h, "%f\n", function[i]);}
		fclose(h);

	return 0;		
	}
}
