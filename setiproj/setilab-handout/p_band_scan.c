#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"

#include <sched.h>    // for processor affinity
#include <unistd.h>   // unix standard apis
#include <pthread.h>  // pthread api



double avg_power(double *data, int num)
{
    int i;
    double ss;
    
    ss=0;
    for (i=0;i<num;i++) { 
	ss += data[i]*data[i];
    }
    
    return ss/num;
}

double max_of(double *data, int num)
{
    double m=data[0];
    int i;
    
    for (i=1;i<num;i++) { 
	if (data[i]>m) { m=data[i]; } 
    }
    return m;
}

double avg_of(double *data, int num)
{
    double s=0;
    int i;
    
    for (i=0;i<num;i++) { 
	s+=data[i];
    }
    return s/num;
}

void remove_dc(double *data, int num)
{
  int i;
  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (i=0;i<num;i++) {
    data[i] -= dc;
  }
}

signal *sig;
double Fs;
char sig_type;
char *sig_file;

double num_bands;
//double start, end;



int        num_threads;     // number of threads we will use
long        num_procs;       // number of processors we will use 
pthread_t *tid;             // array of thread ids
double    *band_power;     // array with the power of each band

struct worker_args{
	int new_num_bands;
	double  bandwidth;
	int filter_order;
	
}arguments;

void *worker_band_scan(void *arg)
{
			
	 
	double filter_coeffs[arguments.filter_order+1];
	long myid = (long)arg;
	int mystart, myend;
	int  band;

	
	cpu_set_t set;
  	CPU_ZERO(&set);
  	CPU_SET(myid%num_procs,&set);
	
	if (sched_setaffinity(0,sizeof(set),&set)<0) { // do it
    		perror("Can't setaffinity");  // hopefully doesn't fail
   		 exit(-1);
  	}

	//remove_dc(sig->data,sig->num_samples);

	mystart=myid*arguments.new_num_bands; // #bands/#threads     //which band/bands this thread has to filter
	myend=(myid+1)*arguments.new_num_bands;
	
	if(num_bands<=myend){
		myend=num_bands;
	}
	//printf("mystart%ld:%d \n",myid,mystart);
	//printf("myend%ld:%d \n",myid,myend);
	
	for (band=mystart;band<myend;band++) { 
		
       		 generate_band_pass(sig->Fs, band*arguments.bandwidth+0.0001, (band+1)*arguments.bandwidth-0.0001,arguments.filter_order, filter_coeffs);

		hamming_window(arguments.filter_order,filter_coeffs);

	// Convolve
		convolve_and_compute_power(sig->num_samples,sig->data,arguments.filter_order,filter_coeffs,&(band_power[band]));  //band_power[band] stores the power of the band

		
	}
	
	
	
	pthread_exit(NULL);
}

int main(int argc, char *argv[]){
	
	if (argc != 8){
		printf("call p_band_scan with 8 \n");
	}
	sig_type = toupper(argv[1][0]);
	sig_file = argv[2];
	Fs = atof(argv[3]);

	arguments.filter_order = atoi(argv[4]);

	num_bands = atoi(argv[5]);
	num_threads=atoi(argv[6]);      // number of threads
  	num_procs=atoi(argv[7]);        // numer of processors to use
	printf("Argc: %d\n", argc);

	assert(Fs>0.0);
	assert(arguments.filter_order>0 && !(arguments.filter_order & 0x1));
	assert(num_bands>0);
	unsigned long long tstart, tend;

	tstart = get_cycle_count();

	switch (sig_type) {
	case 'T':
	    sig = load_text_format_signal(sig_file);
	    break;

	case 'B':
	    sig = load_binary_format_signal(sig_file);
	    break;

	case 'M':
	    sig = map_binary_format_signal(sig_file);
	    break;
	    
	default:
	    printf("Unknown signal type\n");
	    return -1;
    	}
    
    if (!sig) { 
	printf("Unable to load or map file\n");
	return -1;
    }

    sig->Fs=Fs;

	double Fc=(sig->Fs)/2;	
		
   	arguments.bandwidth = Fc / num_bands;	//bandwidth
	arguments.new_num_bands=(num_bands+num_threads-1)/num_threads;
	printf("%d \n",arguments.new_num_bands);

	tid = (pthread_t *) malloc(sizeof(pthread_t)*num_threads);

	band_power = (double *) malloc(sizeof(double)*num_bands);
	
	remove_dc(sig->data,sig->num_samples);
	double signal_power;
    	signal_power = avg_power(sig->data,sig->num_samples);

    	printf("signal average power:     %lf\n", signal_power);
	long i;
	long rc;
	

//-----------------------------------------------------------

  for (i=0;i<num_threads;i++) {
    rc=pthread_create( &(tid[i]), // thread id gets put here
		       NULL,      // use default attributes
		       worker_band_scan,    // thread will begin in this function
		       (void*) i  // we'll give it i as the argument
	             );
    if (rc!=0) { 
      perror("Failed to start thread");
      exit(-1);
    }
	//printf("create tid[%ld]\n",i);
  }

  // now we will join all the threads
  for (i=0;i<num_threads;i++) {
	
    rc=pthread_join(tid[i],NULL);   // 
	//printf("join tid[%ld]\n",i);
    if (rc!=0) { 
	perror("join failed");
	exit(-1);
    }
	
  }

	tend=get_cycle_count();

#define MAXWIDTH 40

#define THRESHOLD 2.0

#define ALIENS_LOW   50000.0
#define ALIENS_HIGH  150000.0

//print power bands 
	int x;
	double start=-1;
	double end=-1;
	int wow=0;
	//double max_band_power = max_of(band_power,num_bands);
	double avg_band_power = avg_of(band_power,num_bands);
	for(x=0;x<num_bands;x++){
		double band_low = x*arguments.bandwidth+0.0001;
     		double band_high = (x+1)*arguments.bandwidth-0.0001;
		printf("%f   ",band_power[x]);

		if ( (band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
	 	  (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) { 

	// band of interest

			if (band_power[x] > THRESHOLD * avg_band_power) { 
	 		 	printf("(WOW)");
				wow=1;
	 		 	if (start<0) { start=x*arguments.bandwidth+0.0001; }
	  			end = (x+1)*arguments.bandwidth-0.0001;
			}
		}
		printf(" \n");
	}
	printf("seconds %lf \n",cycles_to_seconds(tend-tstart));
	if (wow==1) { 
		printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n",start,end,(end+start)/2.0);
    	} else {
		printf("no aliens\n");
    	}
	free_signal(sig);

	return 0;
}
