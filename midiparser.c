//
//  midiparser.c
//  Adapted from helloring.c by Steve Philbert
//
//  Created by Gavin Baker on 1/27/22.
//

//to compile: gcc -o midi midiparser.c -lportaudio -lportmidi -lfftw3f -lm

#include <stdio.h>
#include <string.h>
#include <portaudio.h>
#include <portmidi.h>
#include <porttime.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

PaStream *audioStream;

/*
 for the harmonizer function
 */

//autocorrelation ffts and stuff
fftwf_plan autocorr_fft;
fftwf_plan autocorr_ifft;

fftwf_complex autocorr_spectrum[(2048/2) + 1];
fftwf_complex autocorr_out_complex[(2048/2) + 1];
float autocorr_frame_padded[2048];
float autocorr_spec_mags[1025];
float autocorr_spec_mags_padded[2048];
float autocorr_out[1024];

//define globals: fft plans, arrays

fftwf_plan fft; //spectrum
fftwf_plan ifft; //audio out

//these need to be malloc()ed *Done!
fftwf_complex* spectrum;
fftwf_complex* output_complex;

float spectrum_mags[(44100/2)+1];

//set binsize=100
float local_maxes[441];
float local_maxes_interp[44100];

float zeros[44100];

float output_mags[(44100/2)+1];
float padded_out[44100];
float padded_in[44100];

/*
 end globals for harmonizer function
 */

/*
 COLA arrays
 */
float input_concat[1024];
float output_concat[1024];
//float earlyhalf[512];
float latehalf[512];
float upslope_window[512];
float downslope_window[512];

//assume filled with -1 for empty cell
float midinotes[8];

float mtof(int midi){
    float step1 = midi - 69;
    float step2 = step1/12;
    return pow(2,step2) * 440;
    //return pow(2,((midi - 69)/12)) * 440;
}

//insert
void insertnote(int note){
    //printf("%i\n",note);
    for(int i=0;i<8;i++){
        if(midinotes[i] == -1){
            midinotes[i] = mtof(note);
            break;
        }
    }
}

//delete
void removenote(int note){
    float freq = mtof(note);
    for(int i=0;i<8;i++){
        if(midinotes[i] == freq){
            midinotes[i] = -1;
            //break;
        }
    }
}

//align left
void alignleft(){
    for(int i=0;i<8;i++){
        if(midinotes[i] == -1){
            for(int j=i;j<8;j++){
                if(midinotes[j] != -1){
                    midinotes[i] = midinotes[j];
                    //midipositions[i] = midipositions[j];
                    midinotes[j] = -1;
                    //midipositions[j] = 0;
                    break;
                }
            }
        }
    }
}

float complex_mag(fftwf_complex input){
    float real = creal(input);
    float imag = cimag(input);
    float mag = sqrt(pow(real,2) + pow(imag,2));
    return mag;
}

float complex_phase(fftwf_complex input){
    float real = creal(input);
    float imag = cimag(input);
    float phase = atan2(imag,real);
    return phase;
}

fftwf_complex polar_to_cart(float mag, float phase){
    float real = mag * cos(phase);
    float imag = mag * sin(phase);
    return (fftwf_complex) {real,imag};
}

void zero_pad(float* input_ary, float* zeros_ary, float* output_ary, int input_size, int tot_size){
    memcpy(output_ary, input_ary, input_size * sizeof(float));
    memcpy(output_ary+input_size, zeros_ary, (tot_size-input_size)*sizeof(float));
}

//input is n_bins long
void interp_local_maxes(float* input, float* output, int binsize, int n_bins, int input_size){
    float prev_peak = 0;
    float next_peak = 0;
    float binsize_float = (float) binsize;
    //int counter = 1;
    for(int i=0;i<n_bins-1;i++){
        prev_peak = input[i];
        next_peak = input[i+1];
        output[i*binsize] = input[i];
        for(int j=1;j<=binsize;j++){
            output[(i*binsize)+j] = prev_peak + ((next_peak - prev_peak)/binsize_float)*((float) j);
        }
    }
    //do the last one, bring it down to 0
    prev_peak = next_peak;
    next_peak = 0;
    for(int j=1;j<=binsize;j++){
        //prev_peak + ((next_peak - prev_peak)/peak_dist)*counter;
        
        //just for safety...
        if(((n_bins-1)*binsize)+j>input_size){
            break;
        }
        
        output[((n_bins-1)*binsize)+j] = prev_peak + ((next_peak - prev_peak)/binsize_float)*((float) j);
    }
}

//find the max out of binsize inputs
float* get_local_maxes(float* input, float* output, int binsize, int inputsize){
    int n_bins = inputsize/binsize;
    for(int i=0;i<n_bins;i++){
        float local_max = 0;
        for(int j=0;j<binsize;j++){
            if(input[(i*binsize)+j] > local_max){
                local_max = input[(i*binsize)+j];
            }
        }
        output[i] = local_max;
    }
    return output;
}

//repitches voice to harmonics of midi-inputted notes
void harmonizer(float input[1024], float output[1024]){
    zero_pad(input,zeros,padded_in,1024,44100);
    
    fftwf_execute(fft);
    
    //now we have spectrum
    
    /*
     This is the new autocorrelation method, see
     stackoverflow.com/questions/4583950/cepstral-analysis-for-pitch-detection
     
     fftwf_complex* autocorr_spectrum;
     fftwf_complex* autocorr_out_complex;
     float autocorr_frame_padded[2048];
     float autocorr_spec_mags[1024];
     float autocorr_spec_mags_padded[2048];
     float autocorr_out[1025]
     */
    
    zero_pad(input,zeros,autocorr_frame_padded,1024,2048);
    fftwf_execute(autocorr_fft);
    for(int i=0;i<(2048/2)+1;i++){
        autocorr_spec_mags[i] = pow((complex_mag(autocorr_spectrum[i])),2);
    }
    autocorr_spec_mags[0] = 0;
    zero_pad(autocorr_spec_mags,zeros,autocorr_spec_mags_padded,(2048/2)+1,2048);
    fftwf_execute(autocorr_ifft);
    float div = complex_mag(autocorr_out_complex[0]);
    autocorr_out[0] = 1;
    for(int i=1;i<(2048/2)+1;i++){
        autocorr_out[i] = (complex_mag(autocorr_out_complex[i]))/div;
    }
    
    //now we want to find the fs/50 to fs/400 peak in autocorr_out
    //here's where we would change the frequency range to tailor to a different voice input
    int max_samps = 44100/50; //(bigger n, lower freq)
    int min_samps = 44100/400; //(smaller n, higher freq)
    
    float max_auto_val = -1;
    int max_auto_index = 0;
    for(int i=min_samps;i<max_samps;i++){
        if(autocorr_out[i]>max_auto_val){
            max_auto_val = autocorr_out[i];
            max_auto_index = i;
        }
    }
    
    float input_freq = 44100/((float)max_auto_index);
    
    /*
     now we have the pitch of the input signal-
     want to find the local maxima array of the spectrum
     */
    
    //find the magnitudes of the spectrum values
    for(int i=0;i<(44100/2)+1;i++){
        spectrum_mags[i] = complex_mag(spectrum[i]);
    }
    
    get_local_maxes(spectrum_mags,local_maxes,100,(44100/2)+1);
    interp_local_maxes(local_maxes,local_maxes_interp,100,220,(44100/2)+1);
    //printf("%f\n",local_maxes_interp[300]);
    
    /*
     we now have a local max array that we can pull from to get output mags
     */
    
    memset(output_mags,0,(44100/2)+1 *sizeof(float));
    
    //float midi_freq = 200;
    for(int j=0;j<8;j++){
        if(midinotes[j] == -1){
            break;
        }
        float freq_diff = input_freq - midinotes[j];
        
        //max 200 harmonics
        for(int i=1;i<200;i++){
            //if we're dealing with infinitys
            if(max_auto_index == 0){
                freq_diff = 0;
            }
            //find the corresponding bin
            int spec_bin = i*(midinotes[j] + freq_diff);
            
            //if we will look for the magnitude of an invalid bin, break
            //spec_bin+1 is so we can place in the bin above and below as well.
            if(spec_bin+1>(44100/2)+1){
                break;
            }
            //if we will place that magnitude in an invalid bin, break
            //+1 is so we can place above and below
            if((i*midinotes[j])+1>(44100/2)+1){
                break;
            }
            //else put it in
            output_mags[i*((int) midinotes[j])] = local_maxes_interp[spec_bin];
            
            //also insert into one below and above
            output_mags[(i*((int) midinotes[j])) -1] = (1/2) * local_maxes_interp[spec_bin];
            output_mags[(i*((int) midinotes[j])) +1] = (1/2) * local_maxes_interp[spec_bin];
        }
    }
    
    for(int i=0;i<(44100/2)+1;i++){
        output_complex[i] = polar_to_cart(output_mags[i],complex_phase(spectrum[i]));
    }
    
    //printf("%f: %f\n",freq, max_auto_val);
    
    fftwf_execute(ifft);
    
    memcpy(output,padded_out,1024*sizeof(float));
    
    for(int i=0;i<1024;i++){
        output[i] = output[i]/44100;
    }
}



//a bit of reverb?

static int midiparserCallback( const void *inputBuffer, void *outputBuffer,
                        unsigned long framesPerBuffer,
                        const PaStreamCallbackTimeInfo* timeInfo,
                        PaStreamCallbackFlags statusFlags,
                        void *userData )
{
    //init stuff
    float *in = (float*)inputBuffer, *out = (float*)outputBuffer;
    
    float input_max = 0;
    for(int i=0;i<512;i++){
        if(in[i]>input_max){
            input_max = in[i];
        }
    }
    
    //copy the back half to the front half
    memcpy(input_concat,input_concat+512,512*sizeof(float));
    //copy the input into the back half
    memcpy(input_concat+512,in,512*sizeof(float));
    
    //takes in the input data and writes (for now) to the output buffer. can change what it writes to.
    //requires 1024 samples in and out!
    harmonizer(input_concat,output_concat);
    
    //COLA latehalf and first half of output_concat
    float output_max = 0;
    for(int i=0;i<512;i++){
        out[i] = (latehalf[i] * downslope_window[i]) + (output_concat[i] * upslope_window[i]);
        if(out[i]>output_max){
            output_max = out[i];
        }
    }
    float scale_factor = input_max/output_max;
    for(int i=0;i<512;i++){
        out[i] = out[i] * scale_factor;
    }
    //save second half of output_concat as latehalf
    memcpy(latehalf,output_concat+512,512*sizeof(float));
    
    return paContinue;
}

//portaudio fxn
void pa_init(){
    int i,id;
    const PaDeviceInfo *info;
    const PaHostApiInfo *hostapi;
    
    PaStreamParameters outputParameters, inputParameters;
    Pa_Initialize();
    
    //input parameters
    for (i=0;i < Pa_GetDeviceCount(); i++) {
        info = Pa_GetDeviceInfo(i);         /* get information from current device */
        hostapi = Pa_GetHostApiInfo(info->hostApi); /*get info from curr. host API */
        
        if (info->maxInputChannels > 0)         /* if curr device supports input */
            printf("%d: [%s] %s (input)\n",i, hostapi->name, info->name );
        
    }
    
    printf("\nType AUDIO input device number: ");
    scanf("%d", &id);                   /* get the output device id from the user */
    info = Pa_GetDeviceInfo(id);       /* get chosen device information structure */
    hostapi = Pa_GetHostApiInfo(info->hostApi);         /* get host API structure */
    printf("Opening AUDIO input device [%s] %s\n", hostapi->name, info->name);
    
    inputParameters.device = id;                             /* chosen device id */
    inputParameters.channelCount = 1;                              /* mono input */
    inputParameters.sampleFormat = paFloat32;    /* 32 bit floating point output */
    inputParameters.suggestedLatency = info->defaultLowOutputLatency;/* set default */
    inputParameters.hostApiSpecificStreamInfo = NULL;        /* no specific info */
    
    //output parameters
    for (i=0;i < Pa_GetDeviceCount(); i++) {
        info = Pa_GetDeviceInfo(i);         /* get information from current device */
        hostapi = Pa_GetHostApiInfo(info->hostApi); /*get info from curr. host API */
        
        if (info->maxOutputChannels > 0)         /* if curr device supports output */
        printf("%d: [%s] %s (output)\n",i, hostapi->name, info->name );
        
    }
    
    printf("\nType AUDIO output device number: ");
    scanf("%d", &id);                   /* get the output device id from the user */
    info = Pa_GetDeviceInfo(id);       /* get chosen device information structure */
    hostapi = Pa_GetHostApiInfo(info->hostApi);         /* get host API structure */
    printf("Opening AUDIO output device [%s] %s\n", hostapi->name, info->name);
    
    outputParameters.device = id;                             /* chosen device id */
    outputParameters.channelCount = 1;                             /* mono output */
    outputParameters.sampleFormat = paFloat32;    /* 32 bit floating point output */
    outputParameters.suggestedLatency = info->defaultLowOutputLatency;/* set default */
    outputParameters.hostApiSpecificStreamInfo = NULL;        /* no specific info */
    
    Pa_OpenStream(               /* open the PaStream object and get its address */
                  &audioStream,      /* get the address of the portaudio stream object */
                  &inputParameters,              /* provide input parameters */
                  &outputParameters, /* provide output parameters */
                  44100,     /* set sampling rate */
                  512,   /* set frames per buffer ********** NOTE- Was 1024 ***********/
                  paClipOff,         /* set no clip */
                  midiparserCallback,    /* provide the callback function address */
                  NULL );            /* provide no data for the callback */
    
    //printf("stream opened, but not started\n");
    Pa_StartStream(audioStream); /* start the callback mechanism */
    printf("Stream started");
}

//portaudio fxn
void terminate_stuff()
{
    Pa_StopStream( audioStream );    /* stop the callback mechanism */
    Pa_CloseStream( audioStream );   /* destroy the audio stream object */
    Pa_Terminate();                  /* terminate portaudio */
}

int main(){
    
    //init a portaudio stream
    //populate global arrays
    memset(zeros,0,44100*sizeof(float));
    
    //init the COLA arrays
    //memset(earlyhalf,0,512*sizeof(float));
    memset(latehalf,0,512*sizeof(float));
    memset(input_concat,0,1024*sizeof(float));
    memset(output_concat,0,1024*sizeof(float));
    float step = 1/512;
    for(int i=0;i<512;i++){
        upslope_window[i] = step*i;
        downslope_window[i] = 1 - (step*i);
    }
    
    //setup the harmonizer spectrum arrays
    spectrum = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * (44100/2)+1);
    output_complex = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * (44100/2)+1);
    
    //autocorrelation stuff
    autocorr_fft = fftwf_plan_dft_r2c_1d(2048,autocorr_frame_padded,autocorr_spectrum,FFTW_MEASURE);
    autocorr_ifft = fftwf_plan_dft_r2c_1d(2048,autocorr_spec_mags_padded,autocorr_out_complex,FFTW_MEASURE);
    
    //init the fft plans
    fft = fftwf_plan_dft_r2c_1d(44100,padded_in,spectrum,FFTW_MEASURE);
    ifft = fftwf_plan_dft_c2r_1d(44100,output_complex,padded_out,FFTW_MEASURE);
    printf("FFT's Measured\n");
    
    for(int a=0;a<8;a++){
        midinotes[a] = -1;
    }
    
    int cnt, i, dev;
    PmError retval;
    const PmDeviceInfo *info;
    PmEvent msg[32];
    PortMidiStream *mstream;
    Pm_Initialize();
    cnt = Pm_CountDevices();
    if(cnt){
        for(i=0;i<cnt;i++){
            info = Pm_GetDeviceInfo(i);
            if(info->input){
                printf("%d %s \n", i, info->name);
            }
        }
        
        printf("choose device: ");
        scanf("%d", &dev);
        Pt_Start(1,NULL,NULL);
        retval = Pm_OpenInput(&mstream, dev, NULL, 512L, NULL, NULL);
        
        if(retval != pmNoError){
            printf("error: %s \n", Pm_GetErrorText(retval));
        }
        else{
            
            //portaudio stuff
            pa_init();
            //printf("hello!");
            //maybe move this stuff into the PA callback? if it works here it's fine, we might have issues with input coming too fast
            //if we're not careful if it gets moved though... If this works, keep it.
            
            //def make a different interrupt
            while(Pt_Time(NULL)<40000){
                if(Pm_Poll(mstream)){
                    cnt = Pm_Read(mstream, msg, 32);
                    for(i=0;i<cnt;i++){
                        
                        if(Pm_MessageStatus(msg[i].message) == 144){
                            insertnote(Pm_MessageData1(msg[i].message));
                            //alignleft();
                        }
                        else if(Pm_MessageStatus(msg[i].message) == 128){
                            removenote(Pm_MessageData1(msg[i].message));
                            alignleft();
                        }
                        
                        //print out the array if we want to do that, otherwise the array gets maintained backstage i guess
                        //for(int k=0;k<8;k++){
                        //    printf("%d ",midinotes[k]);
                        //}
                        //printf("\n");
                        
                    }
                }
            }
        }
        Pm_Close(mstream);
    }
    else{
        printf("Nothing to see here!");
    }
    // portaudio stuff
    terminate_stuff();
    
    Pm_Terminate();
    return 0;
}
