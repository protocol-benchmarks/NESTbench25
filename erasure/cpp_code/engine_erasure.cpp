#include "engine_erasure.h"

#include <random>
static std::mt19937_64 rng_engine; // default-constructed engine
static std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
static std::normal_distribution<double> normal_dist;  // global object

//information erasure using an underdamped oscillator
//Compiles as C++ for I/O and stdlib compatibility, but uses mostly C-style logic

//Build instructions (via Makefile):
//  make standalone    # compiles as stand-alone executable 'sim' (requires uncommenting main())
//  make library       # compiles as static library 'libengine_erasure.a'

//Note:
//  - Only functions declared in engine_erasure.h are exposed from the library
//  - External code should include the header: #include "engine_erasure.h"
//  - To link against the library in your program:
//      g++ your_program.c -L. -lengine_erasure -o your_program

// NOTE ON INDENTATION:
// This file uses minimal or no indentation as a matter of personal preference
// If you'd prefer conventional formatting, most editors can reindent automatically:
//   - clang-format:   clang-format -i engine_erasure.c
//   - astyle:         astyle --style=kr --indent=spaces=4 engine_erasure.c
//   - VS Code:        Right-click → Format Document
//   - Emacs:          M-x indent-region
//   - Vim:            gg=G
//
// Feel free to reindent locally as needed.

#include <cmath>
#include <ctime>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unistd.h>

//minimal std imports
using std::ios;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;

static char st[256];

//mmm, pi
static const double pi=4.0*atan(1.0);

//model parameters
static double eff_0=1090; //in Hz
static double quality=7.0;
static double omega_0=(2.0*pi*eff_0)/(1e6); //in reciprocal microseconds
static double cycle_time=1e6/eff_0; //1/f_0, in microseconds

//time parameters (microseconds)
static double tau;
static double timestep=0.1;

static const int number_of_protocol_timepoints=1000;

static double trajectory_time=0.5*cycle_time;
static double trajectory_time_quiescent=1.0*cycle_time; //wait for one period at end to check stability
static double protocol_update_time=trajectory_time/(1.0*number_of_protocol_timepoints);
static int tau_last_protocol_update;

//control parameters
static const int number_of_control_parameters=2;
static double cee_instantaneous[number_of_control_parameters];
static double cee[number_of_protocol_timepoints][number_of_control_parameters];
static double cee_boundary[number_of_control_parameters][2] = {{0.0, 0.0},{5.0, 5.0}};

//microscopic variables
static double position;
static double velocity;

//thermodynamic parameters
static double work;
static double heat;
static double mean_work;
static double mean_heat;
static double mean_prob_success;

//for visualization
static const int number_of_pictures=150;
static const int n_traj_pictures=100000;

static int traj_number=0;
static int q_record_trajectory=0;

static int tau_last_picture_update;
static double position_time[n_traj_pictures][number_of_pictures];
static double picture_update_time=(trajectory_time+trajectory_time_quiescent)/(1.0*number_of_pictures);

//functions void
static void initialize(void);
void equilibrate(void);
void final_answer(void); //100 instances of the order parameter calculated using 10^4 trajectories
void langevin_step(void);
void load_protocol(void);//reads input_control_parameters.dat
void run_trajectory(void);//runs one trajectory
void output_protocol(void);
void simulation_step(void);
void update_potential(void);
void visualize_protocol(void); //visualize the effect of the current protocol
void load_default_protocol(void); //loads default protocol
void output_trajectory_data(void);
void output_potential(int time_slice);
void output_histogram_position(int time_slice);
void make_work_histogram(double *work_values,int n_traj);

//functions int
static int sign(double x);

//functions double
static double potential(void);
static double heaviside(double x);
static double gauss_rv(double sigma);
double calculate_order_parameter(int n_traj); //computes order parameter over n_traj trajectories
double calculate_custom_order_parameter(int n_traj); //for worked example -- does "curriculum learning", i.e. produces an intermediate order parameter until "real" order parameter is large enough


void initialize(void){

//clean up
snprintf(st,sizeof(st),"rm report_*.dat");
cout << st << endl;
system(st);

//seed RN generator
rng_engine.seed(static_cast<uint64_t>(time(NULL)));

}

//
// MAIN_BLOCK_START
// To compile as a library, comment out everything between here, and DOWN THERE vvv
//

int main(void){

initialize();
load_default_protocol();
visualize_protocol();

return 0;

}

//
// To compile as a library, comment out everything betweehn here, and UP THERE ^^^
// MAIN_BLOCK_END
//


double gauss_rv(double sigma){

    normal_dist = std::normal_distribution<double>(0.0, sigma);
    return normal_dist(rng_engine);
    
}


void equilibrate(void){

//PDF(x) Gaussian with variance kT/k = nm
//PDF(v) Gaussian with variance kT/m = omega_0^2 (nm)^2

//initialize position
double initial_state=1.0;
if(drand48()<0.5){initial_state=0.0;}
position=(2.0*initial_state-1.0)*cee_boundary[1][0]+gauss_rv(1.0);

//initialize velocity
velocity=gauss_rv(omega_0);

//reset counters
tau=0.0;
work=0.0;
heat=0.0;
tau_last_protocol_update=-1;
tau_last_picture_update=-1;

//set potential
for(int i=0;i<number_of_control_parameters;i++){cee_instantaneous[i]=cee_boundary[i][0];}

}


void langevin_step(void){

double e1,e2;
double c0=cee_instantaneous[0];
double c1=cee_instantaneous[1];

//initial energy
e1=potential();

//parameters
double lambda=exp(-timestep*omega_0/quality);

//gradient term
double grad=(position-sign(position-c0)*c1);
grad=grad*quality*omega_0; //Q omega_0= k/gamma

//position update
position=position+velocity*timestep;

//velocity update (nm/ms)
velocity=lambda*velocity-(1.0-lambda)*grad+sqrt(1.0-lambda*lambda)*omega_0*gauss_rv(1.0);

//final energy
e2=potential();

//heat increment
heat+=e2-e1;

//increment time
tau+=timestep;

}


double potential(void){

double c0=cee_instantaneous[0];
double c1=cee_instantaneous[1];
double p1=position-sign(position-c0)*c1;
double q1=0.5*p1*p1;

//origin and continuity terms
q1+=2.0*c0*c1*heaviside(position-c0);
q1-=2.0*c0*c1*(1.0-heaviside(c0));

return (q1);

}

int sign(double x){
    return (x > 0.0) - (x < 0.0);
}


double heaviside(double x){

if(x<0.0){return (0.0);}
else{
if(x>0.0){return (1.0);}
else{return (0.5);}
}

}

void update_potential(void){

int protocol_time_index=(int) (tau/protocol_update_time);
if((protocol_time_index > tau_last_protocol_update) && (tau <= trajectory_time)){

//reset counter
tau_last_protocol_update=protocol_time_index;

double e1,e2;

//initial energy
e1=potential();

//update protocol
if(protocol_time_index>=number_of_protocol_timepoints){protocol_time_index=number_of_protocol_timepoints-1;}

for(int i=0;i<number_of_control_parameters;i++){cee_instantaneous[i]=cee[protocol_time_index][i];}

//final energy
e2=potential();

//work
work+=e2-e1;


}}

void simulation_step(void){

update_potential();
langevin_step();
output_trajectory_data();

}

void run_trajectory(void){

double e1,e2;

equilibrate();

while(tau<trajectory_time){simulation_step();}

//final-time potential
e1=potential();
for(int i=0;i<number_of_control_parameters;i++){cee_instantaneous[i]=cee_boundary[i][1];}
e2=potential();
work+=e2-e1;

//quiescent period to check stability of memory
while(tau<trajectory_time+trajectory_time_quiescent){simulation_step();}

//cout << work << " " << heat << " " << position << " " << velocity/omega_0 << endl;

}

void load_default_protocol(void){

int i;
double t1;

for(i=0;i<number_of_protocol_timepoints;i++){

t1=i/(number_of_protocol_timepoints-1.0);

if(t1<0.5){

cee[i][0]=cee_boundary[0][0];
cee[i][1]=cee_boundary[1][0]*(1.0-2.0*t1);

}
else{

cee[i][0]=10.0;
cee[i][1]=2.0*cee_boundary[1][0]*(t1-0.5);

}

}}


void load_protocol(void){

//protocol for 0^+ < t < tf^- (boundary conditions for t=0 and t=tf enforced elsewhere)

int i,j;
snprintf(st, sizeof(st), "input_control_parameters.dat");
ifstream in;

in.open((st));
if(!in){
cerr << "Error: could not open " << (st) << endl;
exit(1);
}

for(i=0;i<number_of_protocol_timepoints;i++){
for(j=0;j<number_of_control_parameters;j++){
in >> cee[i][j];

if(in.fail()){
cerr << "Error reading data from " << st << " at line " << i << ", column " << j << endl;
exit(1);

}}}

in.close();

}


double calculate_order_parameter(int n_traj){

//seed RN generator
rng_engine.seed(static_cast<uint64_t>(time(NULL)));

int i;

double op;
mean_prob_success=0.0;

for(i=0;i<n_traj;i++){

run_trajectory();
if(position<0){mean_prob_success+=1.0/(1.0*n_traj);}

}

op=1.0-mean_prob_success;

return (op);

}


double calculate_custom_order_parameter(int n_traj){

//seed RN generator
rng_engine.seed(static_cast<uint64_t>(time(NULL)));

int i;

double op;
double rms_position_offset=0.0;
mean_prob_success=0.0;

for(i=0;i<n_traj;i++){

run_trajectory();
rms_position_offset+=(position+cee_boundary[1][1])*(position+cee_boundary[1][1])/(1.0*n_traj);
if(position<0){mean_prob_success+=1.0/(1.0*n_traj);}

}

rms_position_offset=sqrt(rms_position_offset);

if(mean_prob_success<0.9){op=rms_position_offset+1.0;}
else{op=1.0-mean_prob_success;}

return (op);

}

void final_answer(void){

snprintf(st,sizeof(st),"report_answer.dat");
ofstream out1(st,ios::out);


int i;
int n_traj=10000;
int n_samples=100;

double op;
double op_local=0.0;
double op2_local=0.0;

double av,var;

for(i=0;i<n_samples;i++){

op=calculate_order_parameter(n_traj);

op_local+=op;
op2_local+=op*op;

av=op_local/(i+1.0);
var=op2_local/(i+1.0)-av*av;
if(i==0){var=0.0;}

out1 << "order parameter = " << av << " ± " << sqrt(var/(i+1.0)) << ", n_samples = " << i+1 << endl;

}
}

void visualize_protocol(void){

int i;

//for work histogram
double *work_values = (double *)malloc(n_traj_pictures * sizeof(double));
if(work_values == NULL){fprintf(stderr, "Error: could not allocate memory for work_values\n");exit(1);}

mean_work=0.0;
mean_heat=0.0;
mean_prob_success=0.0;

q_record_trajectory=1;
traj_number=0;
for(i=0;i<n_traj_pictures;i++){

run_trajectory();

//for histogram
work_values[i]=work;

//for averages
mean_work+=work/(1.0*n_traj_pictures);
mean_heat+=heat/(1.0*n_traj_pictures);
if(position<0){mean_prob_success+=1.0/(1.0*n_traj_pictures);}

traj_number++;

}
q_record_trajectory=0;

//position histograms
for(i=0;i<number_of_pictures;i++){output_histogram_position(i);}

//work histogram
make_work_histogram(work_values,n_traj_pictures);
free(work_values);

//ppotential and Boltzmann weight
for(i=0;i<number_of_pictures;i++){output_potential(i);}

//mean values
cout << "<W>=" << mean_work << ", <Q> = " << mean_heat << ", P= " << mean_prob_success <<endl;

//run python movie and picture maker
snprintf(st,sizeof(st),"python movie.py");
cout << st << endl;
system(st);


}

void output_trajectory_data(void){
int picture_index=(int) (tau/picture_update_time);
if((q_record_trajectory==1) && (picture_index>tau_last_picture_update)){

//reset counter
tau_last_picture_update=picture_index;
if(picture_index>=number_of_pictures){picture_index=number_of_pictures-1;}

position_time[traj_number][picture_index]=position;

}}


void output_histogram_position(int time_slice){

snprintf(st,sizeof(st),"report_position_time_%d.dat",time_slice);
ofstream out1(st,ios::out);

int i;
int nbin=0;
const int bins=50;

double histo[bins];
double histo_width;
double maxpos=0.0;
double minpos=0.0;

//position recorded elswehere
for(i=0;i<bins;i++){histo[i]=0.0;}

for(i=0;i<n_traj_pictures;i++){

if(i==0){maxpos=position_time[i][time_slice];minpos=position_time[i][time_slice];}
else{

if(position_time[i][time_slice]>maxpos){maxpos=position_time[i][time_slice];}
if(position_time[i][time_slice]<minpos){minpos=position_time[i][time_slice];}
}}

//record width
histo_width=(maxpos-minpos)/(1.0*bins);

//safeguard
if(histo_width<1e-6){histo_width=1e-6;}

for(i=0;i<n_traj_pictures;i++){

nbin=(int) (1.0*bins*(position_time[i][time_slice]-minpos)/(maxpos-minpos));
if(nbin==bins){nbin--;}

histo[nbin]+=1.0/(1.0*n_traj_pictures);

}

double hmax=0.0;
for(i=0;i<bins;i++){
if(histo[i]>hmax){hmax=histo[i];}
}

//output
double x1;
for(i=0;i<bins;i++){

if(histo[i]>0.5/(1.0*n_traj_pictures)){
x1=maxpos*i/(1.0*bins)+minpos*(1.0-i/(1.0*bins))+0.5*(maxpos-minpos)/(1.0*bins);

out1 << x1 << " " << histo[i]/histo_width << endl;

}
}

}



void output_potential(int time_slice){

snprintf(st,sizeof(st),"report_boltz_time_%d.dat",time_slice);
ofstream out1(st,ios::out);

snprintf(st,sizeof(st),"report_potential_time_%d.dat",time_slice);
ofstream out2(st,ios::out);

int i;
const int n_points=1000;

double e1;
double e_min=0;
double e_values[n_points];

double x1=-2.0*cee_boundary[1][0];
double x2=2.0*cee_boundary[1][0];
double zed=0.0;
double delta_x;

//store potential parameters
double cee_instantaneous_holder[number_of_control_parameters];
for(i=0;i<number_of_control_parameters;i++){cee_instantaneous_holder[i]=cee_instantaneous[i];}


//time slice
double t1=picture_update_time*time_slice;


if(time_slice==0){
for(i=0;i<number_of_control_parameters;i++){cee_instantaneous[i]=cee_boundary[i][0];}
}
else{
if(t1>=trajectory_time){
for(i=0;i<number_of_control_parameters;i++){cee_instantaneous[i]=cee_boundary[i][1];}
}
else{

int protocol_time_slice=(int) (t1/(1.0*protocol_update_time));
for(i=0;i<number_of_control_parameters;i++){cee_instantaneous[i]=cee[protocol_time_slice][i];

}}}

//enforce final-time potential
if(time_slice==number_of_pictures-1){
for(i=0;i<number_of_control_parameters;i++){cee_instantaneous[i]=cee_boundary[i][1];}
}


//record position
double position_holder=position;

for(i=0;i<n_points;i++){

position=x1+(x2-x1)*i/(n_points-1.0);
e1=potential();

out2 << position << " " << e1 << endl;

}

//boltzmann weight
delta_x=(x2-x1)/(1.0*n_points);
for(i=0;i<n_points;i++){

position=x1+(x2-x1)*i/(1.0*n_points-1.0);
e1=potential();

if(i==0){e_min=e1;}
else{if(e1<e_min){e_min=e1;}}

e_values[i]=e1;

}

//calculate Z
for(i=0;i<n_points;i++){

e_values[i]-=e_min;
zed+=delta_x*exp(-e_values[i]);

}

//plot point
for(i=0;i<n_points;i++){

position=x1+(x2-x1)*i/(1.0*n_points-1.0);
out1 << position << " " << exp(-e_values[i])/zed << endl;

}

//close files
out1.close();
out2.close();

//reset registers
position=position_holder;
for(i=0;i<number_of_control_parameters;i++){cee_instantaneous[i]=cee_instantaneous_holder[i];}


}

void output_protocol(void){

int i,j;
double t1;
ofstream out[number_of_control_parameters];

for(j=0;j<number_of_control_parameters;j++){
snprintf(st,sizeof(st),"report_control_parameters_%d.dat",j);
out[j].open(st,ios::out);
}

//initial values
for(j=0;j<number_of_control_parameters;j++){out[j] << 0 << " " << cee_boundary[j][0] << endl;}

for(i=0;i<number_of_protocol_timepoints;i++){

t1=i/(number_of_protocol_timepoints-1.0);
for(j=0;j<number_of_control_parameters;j++){out[j] << t1 << " " << cee[i][j] << endl;}

}

//final values
for(j=0;j<number_of_control_parameters;j++){out[j] << 1 << " " << cee_boundary[j][1] << endl;}

//close output files
for(j=0;j<number_of_control_parameters;j++){out[j].close();}

}

void make_work_histogram(double *work_values,int n_traj) {

    const int bins = 100;
    double histo[bins] = {0.0};
    double work_min = work_values[0];
    double work_max = work_values[0];

    // find min and max
    for (int i = 1; i < n_traj; i++) {
        if (work_values[i] < work_min) work_min = work_values[i];
        if (work_values[i] > work_max) work_max = work_values[i];
    }

    double bin_width = (work_max - work_min) / bins;
    if (bin_width < 1e-6) bin_width = 1e-6; // safeguard

    // fill histogram
    for (int i = 0; i < n_traj; i++) {
        int bin = (int)((work_values[i] - work_min) / (work_max - work_min) * bins);
        if (bin == bins) bin--; // include upper edge
        histo[bin] += 1.0 / (n_traj * bin_width);
    }

    // output
    snprintf(st, sizeof(st), "report_work_histogram.dat");
    ofstream out(st, ios::out);
    for (int i = 0; i < bins; i++) {
        double center = work_max * i / (double)bins + work_min * (1.0 - i / (double)bins) + 0.5 * (work_max - work_min) / bins;
        if (histo[i] > 0) {
            out << center << " " << histo[i] << endl;
        }
    }
    out.close();
}
