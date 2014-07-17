


struct kernel{
	//Convolution kernel. Contains the impulse response.
	int size;
	struct sample *next;
	
};

struct signal{
	int size;
	struct sample *next;
};

struct sample{
	//Singly linked list element containing a single sample.
	int nummy;
	struct sample *next;
	
};

double sinc(double input);//Just a basic sinc function.
double* generatekernel(double minf, double maxf, double minx, double maxx, int numpoints);//Make the sinc filter kernel.
double getsample (int n); //Access the n'th variable in the sample.
double convolution_singlevalue (struct kernel* mykernel, struct signal* mysignal); //Convolution for one value of t. A complete function will call this for each sample in signal.

/*
We need the following functions:
-Convolution function to convolve the signal with kernel
-Generator function to generate the convo kernel
-Blackman window function (?)
-Sinc function

*/
