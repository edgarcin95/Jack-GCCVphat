//make: gcc -o GCCVphat GCCVphatV4.c -ljack -lfftw3 -lm
//Edgar Rigoberto Santos martín
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <jack/jack.h>
#include <complex.h>
#include <fftw3.h>

double complex *i_fft, *i_time, *o_fft, *o_time, *i2_fft, *i2_time, *o2_fft, *o2_time;
fftw_plan i_forward, o_inverse, i2_forward, o2_inverse;

jack_port_t *input_port;
jack_port_t *input_port2;
jack_port_t *output_port;
jack_port_t *output_port2;
jack_client_t *client;

double sample_rate, *hann;
float dist;

int jack_callback (jack_nframes_t nframes, void *arg){
	jack_default_audio_sample_t *in, *out, *in2, *out2, xMag=0, yMag=0;
	int i;
	
	in = (jack_default_audio_sample_t *)jack_port_get_buffer (input_port2, nframes);
	in2 = (jack_default_audio_sample_t *)jack_port_get_buffer (input_port, nframes);
	out = (jack_default_audio_sample_t *)jack_port_get_buffer (output_port, nframes);
	out2 = (jack_default_audio_sample_t *)jack_port_get_buffer (output_port2, nframes);
	
	float xMean = 0, yMean = 0;
	
	for(i=0; i<nframes; ++i){
		xMean += in[i];
		yMean += in2[i];
	}
	
	if(abs(xMean)>=1 || abs(yMean)>=1){
		//Obteniendo el promedio
		xMean /= nframes;
		yMean /= nframes;
		
		//A cada señal original se le resta el promedio y se aplica la ventana de Hann antes de FFT
		for(i=0; i<nframes; ++i){
			i_time[i] = (in[i] - xMean)*hann[i];
			i2_time[i]= (in2[i]- yMean)*hann[i];
		}
		
		fftw_execute(i_forward);
		fftw_execute(i2_forward);

		for(i=0; i<nframes; ++i){ //aplicamos la transformada de fase
			o_fft[i] = i_fft[i] * conj(i2_fft[i]) / (cabs(i_fft[i]*conj(i2_fft[i])));
		}
		
		//Regresando al dominio del tiempo
		fftw_execute(o_inverse);
		double ccvp[nframes];
		for(i=0; i<nframes; i++){//fftw3 requiere normalizar su salida real de esta manera
			ccvp[i] = (creal(o_time[i]))/nframes;
		}
		
		//Encontrando el desfase y el ángulo
		int desfase=0, idx=0;
		float anglerad, max=0;
		for(i=0; i<nframes; ++i){
			if(ccvp[i]>max) {
				max = ccvp[i];
				idx = i;
			}
		}
		if(idx < nframes/2){
			desfase = idx;
		}
		else{
			desfase = idx-nframes;
		}
		anglerad = asin(((desfase/sample_rate) * (344)) / dist);
		//printf ("desfase = %d\n", desfase);
		printf ("Ángulo: %f°\n", anglerad*180/M_PI);
	}
	
	for(i = 0; i < nframes; ++i){
		out[i] = in[i];
		out2[i] = in2[i];
	}
	
	return 0;
}

void jack_shutdown (void *arg){
	exit (1);
}


int main (int argc, char *argv[]) {
	if(argc < 2){
		printf("Debes darme la distancia entre micrófonos.\n");
		exit(1);
	}
	dist = atof(argv[1]);
	
	const char *client_name = "GCCVphat";
	jack_options_t options = JackNoStartServer;
	jack_status_t status;
	int i;
	
	/* open a client connection to the JACK server */
	client = jack_client_open (client_name, options, &status);
	if (client == NULL){
		/* if connection failed, say why */
		printf ("jack_client_open() failed, status = 0x%2.0x\n", status);
		if (status & JackServerFailed) {
			printf ("Unable to connect to JACK server.\n");
		}
		exit (1);
	}
	
	/* if connection was successful, check if the name we proposed is not in use */
	if (status & JackNameNotUnique){
		client_name = jack_get_client_name(client);
		printf("Warning: other agent with our name is running, `%s' has been assigned to us.\n",client_name);
	}
	
	/* tell the JACK server to call 'jack_callback()' whenever there is work to be done. */
	jack_set_process_callback (client, jack_callback, 0);
	
	
	/* tell the JACK server to call 'jack_shutdown()' if it ever shuts down,
	   either entirely, or if it just decides to stop calling us. */
	jack_on_shutdown (client, jack_shutdown, 0);
	
	
	/* display the current sample rate. */
	printf ("Sample rate: %d\n", jack_get_sample_rate (client));
	printf ("Window size: %d\n", jack_get_buffer_size (client));
	sample_rate = (double)jack_get_sample_rate(client);
	int nframes = jack_get_buffer_size (client);
	
	//preparing FFTW3 buffers
	i_fft = (double complex *) fftw_malloc(sizeof(double complex) * nframes);
	i2_fft = (double complex *) fftw_malloc(sizeof(double complex) * nframes);
	i_time = (double complex *) fftw_malloc(sizeof(double complex) * nframes);
	i2_time = (double complex *) fftw_malloc(sizeof(double complex) * nframes);
	o_fft = (double complex *) fftw_malloc(sizeof(double complex) * nframes);
	o2_fft = (double complex *) fftw_malloc(sizeof(double complex) * nframes);
	o_time = (double complex *) fftw_malloc(sizeof(double complex) * nframes);
	o2_time = (double complex *) fftw_malloc(sizeof(double complex) * nframes);
	
	i_forward = fftw_plan_dft_1d(nframes, i_time, i_fft , FFTW_FORWARD, FFTW_MEASURE);
	i2_forward = fftw_plan_dft_1d(nframes, i2_time, i2_fft , FFTW_FORWARD, FFTW_MEASURE);
	o_inverse = fftw_plan_dft_1d(nframes, o_fft , o_time, FFTW_BACKWARD, FFTW_MEASURE);
	o2_inverse = fftw_plan_dft_1d(nframes, o2_fft , o2_time, FFTW_BACKWARD, FFTW_MEASURE);
	
	hann = (double *) malloc(sizeof(double) * nframes);
	for(i=0; i<nframes; ++i){
		hann[i] = 0.5 - 0.5*cos((2*M_PI*i)/(nframes-1));
		//printf("han %f\t", hann[i]);
	}
	
	/* create the agent input port */
	input_port = jack_port_register (client, "input1", JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);
	input_port2 = jack_port_register (client, "input2", JACK_DEFAULT_AUDIO_TYPE,JackPortIsInput, 0);
	
	/* create the agent output port */
	output_port = jack_port_register (client, "output1",JACK_DEFAULT_AUDIO_TYPE,JackPortIsOutput, 0);
	output_port2 = jack_port_register (client, "output2",JACK_DEFAULT_AUDIO_TYPE,JackPortIsOutput, 0);
	
	/* check that both ports were created succesfully */
	if ((input_port == NULL) || (output_port == NULL)) {
		printf("Could not create agent ports. Have we reached the maximum amount of JACK agent ports?\n");
		exit (1);
	}
	
	if ((input_port2 == NULL) || (output_port2 == NULL)) {
		printf("Could not create agent ports. Have we reached the maximum amount of JACK agent ports?\n");
		exit (1);
	}
	/* Tell the JACK server that we are ready to roll.
	   Our jack_callback() callback will start running now. */
	if (jack_activate (client)) {
		printf ("Cannot activate client.");
		exit (1);
	}
	
	printf ("Agent activated.\n");
	
	printf ("Connecting ports... ");
	 
	/* Assign our input port to a server output port*/
	// Find possible output server port names
	const char **serverports_names;
	/*serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsOutput);
	if (serverports_names == NULL) {
		printf("No available physical capture (server output) ports.\n");
		exit (1);
	}
	// Connect the first available to our input port
	if (jack_connect (client, serverports_names[0], jack_port_name (input_port))) {
		printf("Cannot connect input port 1.\n");
		exit (1);
	}
	if (jack_connect (client, serverports_names[1], jack_port_name (input_port2))) {
		printf("Cannot connect input port 2.\n");
		exit (1);
	}
	// free serverports_names variable for reuse in next part of the code
	free (serverports_names);*/
	
	/* Assign our output port to a server input port*/
	// Find possible input server port names
	serverports_names = jack_get_ports (client, NULL, NULL, JackPortIsPhysical|JackPortIsInput);
	if (serverports_names == NULL) {
		printf("No available physical playback (server input) ports.\n");
		exit (1);
	}
	// Connect the first available to our output port
	if (jack_connect (client, jack_port_name (output_port), serverports_names[0])) {
		printf ("Cannot connect output port 1.\n");
		exit (1);
	}
	if (jack_connect (client, jack_port_name (output_port2), serverports_names[1])) {
		printf ("Cannot connect output port 2.\n");
		exit (1);
	}
	// free serverports_names variable, we're not going to use it again
	free (serverports_names);
	
	
	printf ("done.\n");
	/* keep running until stopped by the user */
	sleep (-1);
	
	
	/* this is never reached but if the program
	   had some other way to exit besides being killed,
	   they would be important to call.
	*/
	jack_client_close (client);
	exit (0);
}
