#include "engine_abp.h"

#include <random>
static std::mt19937_64 rng_engine; // default-constructed engine
static std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
static std::normal_distribution<double> normal_dist;  // global object

//abp in a trap
//Compiles as C++ for I/O and stdlib compatibility, but uses mostly C-style logic

//to compile as static library: g++ -c engine_abp.c -o engine_abp.o; ar rcs libengine_abp.a engine_abp.o
//note that only the functions in engine_abp.h are "exposed"

//external code must use the header #include "engine_abp.h"
//compile your program with g++ your_program.c -L. -lengine_abp -o your_program

//to compile as stand-alone code, restore main function (see below), and compile as e.g.
//g++ -Wall -o sim engine_abp.c -lm -O

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

//time parameters (microseconds)
static double tau;
static double timestep=0.001;

static const int number_of_protocol_timepoints=1000;

static double trajectory_time=0.44;
static double protocol_update_time=trajectory_time/(1.0*number_of_protocol_timepoints);
static int tau_last_protocol_update;

//control parameters (kappa, lambda)
static const int number_of_control_parameters=2;
static double cee_instantaneous[number_of_control_parameters];
static double cee[number_of_protocol_timepoints][number_of_control_parameters];
static double cee_boundary[number_of_control_parameters][2] = {{4.0,4.0},{0.0, 5.0}};

//microscopic variables
static double angle;
static double position[2];

//thermodynamic parameters
static double work;
static double mean_work;
static double delta_global;

//for visualization
static const int number_of_pictures=150;
static const int n_traj_pictures=100000;

static int traj_number=0;
static int q_record_trajectory=0;

static int tau_last_picture_update;
static double** position_time = nullptr;
static double picture_update_time=trajectory_time/(1.0*number_of_pictures);

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


//functions double
static double potential(void);
static double gauss_rv(double sigma);
double calculate_order_parameter(int n_traj); //computes order parameter over n_traj trajectories
double histogram_l2_error(double *r_values, int n_traj);


void initialize(void){
//clean up
snprintf(st,sizeof(st),"rm report_*.dat");
cout << st << endl;
system(st);

//seed RN generator
rng_engine.seed(static_cast<uint64_t>(time(NULL)));

}

//restore main to compile as stand-alone code


/*
int main(void){

initialize();
load_default_protocol();
output_protocol();
//visualize_protocol();
final_answer();

return 0;

}
*/

double gauss_rv(double sigma){

    normal_dist = std::normal_distribution<double>(0.0, sigma);
    return normal_dist(rng_engine);
    
}


void equilibrate(void){

//initialize position
for(int i=0;i<2;i++){position[i]=gauss_rv(1.0/sqrt(cee_boundary[0][0]));}

//initialize angle
angle=2.0*pi*uniform_dist(rng_engine);

//reset counters
tau=0.0;
work=0.0;
tau_last_protocol_update=-1;
tau_last_picture_update=-1;

//set potential
for(int i=0;i<number_of_control_parameters;i++){cee_instantaneous[i]=cee_boundary[i][0];}

}


void langevin_step(void){

//double e1,e2;
double kappa=cee_instantaneous[0];
double lambda=cee_instantaneous[1];

//initial energy
//e1=potential();

//swim terms
double swim_x=lambda*cos(angle)*timestep;
double swim_y=lambda*sin(angle)*timestep;

//gradient terms
double grad_x=-kappa*position[0]*timestep;
double grad_y=-kappa*position[1]*timestep;

//noise terms
double noise_x=sqrt(2.0*timestep)*gauss_rv(1.0);
double noise_y=sqrt(2.0*timestep)*gauss_rv(1.0);

//position updates
position[0]+=swim_x+grad_x+noise_x;
position[1]+=swim_y+grad_y+noise_y;

//angle update (wrap to [0,2pi)
angle+=sqrt(2.0*timestep)*gauss_rv(1.0);
angle=fmod(angle,2.0*pi);
if(angle<0){angle+=2.0*pi;}

//final energy
//e2=potential();

//increment time
tau+=timestep;

}


double potential(void){

double kappa=cee_instantaneous[0];
double q1=0.5*kappa*(position[0]*position[0]+position[1]*position[1]);

return (q1);

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

}

void load_default_protocol(void){

for(int i=0;i<number_of_protocol_timepoints;i++){

cee[i][0]=cee_boundary[0][0];
cee[i][1]=cee_boundary[1][1];

}

}


void load_protocol(void){

//protocol for 0^+ < t < tf^- (boundary conditions for t=0 and t=tf enforced elsewhere)

int i,j;
double dummy; //dummy variable

snprintf(st, sizeof(st), "input_control_parameters.dat");
ifstream in;

in.open((st));
if(!in){
cerr << "Error: could not open " << (st) << endl;
exit(1);
}

for(i=0;i<number_of_protocol_timepoints;i++){
for(j=0;j<number_of_control_parameters;j++){

in >> dummy;

if(j==0){
//kappa limits
if(dummy<1.0){dummy=1.0;}
if(dummy>7.0){dummy=7.0;}
}
else{
//lambda limits
if(dummy<0.0){dummy=0.0;}
if(dummy>11.0){dummy=11.0;}
}

cee[i][j]=dummy;

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
std::vector<double> r_values(n_traj);

double op;
double delta;
mean_work=0.0;

//record values of r(tf)
for(i=0;i<n_traj;i++){

run_trajectory();
r_values[i]=sqrt(position[0]*position[0]+position[1]*position[1]);
mean_work+=work/(1.0*n_traj);

}

//calculate ||P(r(tf)-P_target(r)||
delta=histogram_l2_error(r_values.data(),n_traj);
delta_global=delta;

op=delta+0.001*mean_work;


return (op);

}


void final_answer(void){

snprintf(st,sizeof(st),"report_answer.dat");
ofstream out1(st,ios::out);

int n_traj=1000000;

double op=calculate_order_parameter(n_traj);
out1 << " order parameter = " << op << ", Delta = " << delta_global << ", mean_work = " << mean_work << endl;

}


void visualize_protocol(void){
int i;

// Allocate position_time dynamically
position_time = new double*[n_traj_pictures];
for (int i = 0; i < n_traj_pictures; ++i) {
    position_time[i] = new double[number_of_pictures]();
}

//for work histogram
double *work_values = (double *)malloc(n_traj_pictures * sizeof(double));
if(work_values == NULL){fprintf(stderr, "Error: could not allocate memory for work_values\n");exit(1);}

mean_work=0.0;

q_record_trajectory=1;
traj_number=0;
for(i=0;i<n_traj_pictures;i++){

run_trajectory();

//for histogram
work_values[i]=work;

//for averages
mean_work+=work/(1.0*n_traj_pictures);

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
cout << "<W>=" << mean_work << endl;

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

position_time[traj_number][picture_index]=sqrt(position[0]*position[0]+position[1]*position[1]);

}}


void output_histogram_position(int time_slice){

//Plots G(r) = r P(r), so that \int dr G(r) = 1

snprintf(st,sizeof(st),"report_position_time_%d.dat",time_slice);
ofstream out1(st,ios::out);

int i;
int nbin=0;
const int bins=100;

double histo[bins];
double histo_width;
double maxpos=5;
double minpos=0.001;

//position recorded elswehere
for(i=0;i<bins;i++){histo[i]=0.0;}

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
x1 = minpos + (i + 0.5) * (maxpos - minpos) / bins;

out1 << x1 << " " << histo[i]/histo_width << endl;

}
}

}



void output_potential(int time_slice){

//Boltzmann distribution is G(R) = r P(r), so that \int dr G(r) = 1

snprintf(st,sizeof(st),"report_boltz_time_%d.dat",time_slice);
ofstream out1(st,ios::out);

snprintf(st,sizeof(st),"report_potential_time_%d.dat",time_slice);
ofstream out2(st,ios::out);

int i;
const int n_points=1000;

double e1;
double e_min=0;
double e_values[n_points];

double r1;
double x1=0.001;
double x2=5.0/sqrt(2.0);
double zed=0.0;
double delta_r;

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
double position_holder[2];
position_holder[0]=position[0];
position_holder[1]=position[1];

for(i=0;i<n_points;i++){

r1=x1+(x2-x1)*i/(n_points-1.0);
position[0]=r1/sqrt(2.0);
position[1]=r1/sqrt(2.0);
e1=potential();

out2 << sqrt(position[0]*position[0]+position[1]*position[1]) << " " << e1 << endl;

}

//boltzmann weight
delta_r=(x2-x1)/(1.0*n_points);
for(i=0;i<n_points;i++){

r1=x1+(x2-x1)*i/(n_points-1.0);
position[0]=r1/sqrt(2.0);
position[1]=r1/sqrt(2.0);

e1=potential();

if(i==0){e_min=e1;}
else{if(e1<e_min){e_min=e1;}}

e_values[i]=e1;

}

//calculate Z
for(i=0;i<n_points;i++){

r1=x1+(x2-x1)*i/(n_points-1.0);
position[0]=r1/sqrt(2.0);
position[1]=r1/sqrt(2.0);

e_values[i]-=e_min;
zed+=delta_r*r1*exp(-e_values[i]);

}

//plot point
for(i=0;i<n_points;i++){

r1=x1+(x2-x1)*i/(n_points-1.0);
position[0]=r1/sqrt(2.0);
position[1]=r1/sqrt(2.0);

out1 << r1 << " " << r1*exp(-e_values[i])/zed << endl;

}

//close files
out1.close();
out2.close();

//reset registers
position[0]=position_holder[0];
position[1]=position_holder[1];
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



double histogram_l2_error(double *r_values, int n_traj){

    const int bins = 100;
    const double minpos = 0.001;
    const double maxpos = 5.0;
    const double histo_width = (maxpos - minpos) / bins;

    double histo[bins] = {0.0};

    // Fill histogram
    for (int i = 0; i < n_traj; i++) {
        int bin = (int)((r_values[i] - minpos) / (maxpos - minpos) * bins);
        if (bin < 0 || bin >= bins) continue;
        histo[bin] += 1.0 / (n_traj * histo_width);
    }

    // Read target histogram
    snprintf(st, sizeof(st), "target_p_r.dat");
    ifstream in_target(st);
    if (!in_target) {
        cerr << "Error: could not open " << st << endl;
        exit(1);
    }

    std::vector<std::pair<double, double> > target;
    double x, y;
    while (in_target >> x >> y) {
        target.emplace_back(x, y);
    }
    in_target.close();

    // Compare only entries with matching positions
    double l2_norm = 0.0;
    int matched = 0;
    for (int i = 0; i < bins; i++) {
        double x1 = minpos + (i + 0.5) * histo_width;
        for (size_t j = 0; j < target.size(); ++j) {
            double xt = target[j].first;
            double yt = target[j].second;
            if (fabs(xt - x1) < 1e-8) {
                double diff = histo[i] - yt;
                l2_norm += diff * diff;
                matched++;
                break;
            }
        }
    }

    if(matched==0){return 100.0;}
    return (l2_norm/matched);
}
