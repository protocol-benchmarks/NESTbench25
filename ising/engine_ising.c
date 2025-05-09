#include "engine_ising.h"

#include <random>
static std::mt19937_64 rng_engine;
static std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
static std::normal_distribution<double> normal_dist;

//ising model undergoing state change
//Compiles as C++ for I/O and stdlib compatibility, but uses mostly C-style logic

//Build instructions (via Makefile):
//  make standalone    # compiles as stand-alone executable 'sim' (requires uncommenting main())
//  make library       # compiles as static library 'libengine_ising.a'

//Note:
//  - Only functions declared in engine_ising.h are exposed from the library
//  - External code should include the header: #include "engine_ising.h"
//  - To link against the library in your program:
//      g++ your_program.c -L. -lengine_ising -o your_program

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

//model parameters
static const int lx=32;
static const int number_of_sites=lx*lx;

//time parameters (microseconds)
static int tau;
static const int number_of_protocol_timepoints=1000;

static int trajectory_time=100*number_of_sites; //100 sweeps
static double protocol_update_time=(1.0*trajectory_time)/(1.0*number_of_protocol_timepoints);
static int tau_last_protocol_update;

//control parameters (temp, mag. field)
static const int number_of_control_parameters=2;
static double cee_instantaneous[number_of_control_parameters];
static double cee[number_of_protocol_timepoints][number_of_control_parameters];
static double cee_boundary[number_of_control_parameters][2]={{0.70,0.70},{-1.0,1.0}};

//microscopic variables
static double magnetization;
static double internal_nrg;
static int spin[number_of_sites];
static int neighbors[number_of_sites][4];
static int spin_sequence[number_of_sites];

//thermodynamic parameters
static double entprod;
static double mean_entprod;
static double mean_magnetization;


//for visualization
static const int number_of_pictures=150;
static const int n_traj_pictures=10000;

static int traj_number=0;
static int q_record_trajectory=0;

static int tau_last_picture_update;
static double magnetization_time[n_traj_pictures][number_of_pictures];
static double picture_update_time=(1.0*trajectory_time)/(1.0*number_of_pictures);

//functions void
static void initialize(void);
void final_answer(void); //100 instances of the order parameter calculated using 10^4 trajectories
void glauber_step(void);
void load_protocol(void);//reads input_control_parameters.dat
void run_trajectory(void);//runs one trajectory
void output_protocol(void);
void simulation_step(void);
void update_potential(void);
void initialize_lattice(void);
void visualize_protocol(void); //visualize the effect of the current protocol
void load_default_protocol(void); //loads default protocol
void set_initial_conditions(void);
void output_trajectory_data(void);
void output_lattice(int picture_index);
void output_histogram_magnetization(int time_slice);
void make_entprod_histogram(double *entprod_values,int n_traj);

//functions int
int choose_spin(int t1);
int node_x_coord(int i);
int node_y_coord(int i);

//functions double
double delta_nrg(int chosen_site);
double calculate_order_parameter(int n_traj); //computes order parameter over n_traj trajectories

void initialize(void){

//clean up
snprintf(st,sizeof(st),"rm report_*.dat");
cout << st << endl;
system(st);

//seed RN generator
rng_engine.seed(static_cast<uint64_t>(time(NULL)));

//initialize lattice
initialize_lattice();

}

//restore main to compile as stand-alone code

/*
int main(void){

//seed the RN generator
initialize();

//load default protocol
load_default_protocol();

//write current protocol to file
output_protocol();

cout << calculate_order_parameter(10000) << endl;

//calculate averages and distributions
visualize_protocol();

return 0;

}
*/

void set_initial_conditions(void){

int i;

//set potential
for(i=0;i<number_of_control_parameters;i++){cee_instantaneous[i]=cee_boundary[i][0];}

double magnetic_field=cee_instantaneous[1];

//lattice
magnetization=-1.0;
internal_nrg=-2.0*number_of_sites-magnetic_field*number_of_sites;
for(i=0;i<number_of_sites;i++){spin[i]=-1;}

//reset counters
tau=0;
entprod=0.0;
tau_last_protocol_update=-1;
tau_last_picture_update=-1;


}


void glauber_step(void){

//set temperature
double beta=1.0/cee_instantaneous[0];

//choose spin (sequential)
int chosen_site = choose_spin(tau);

//glauber flip
double de=delta_nrg(chosen_site);
double glauber=exp(-beta*de);
glauber=glauber/(glauber+1.0);

double prob1=uniform_dist(rng_engine);


if(prob1<glauber){

entprod-=beta*de;
internal_nrg+=de;
magnetization-=2.0*spin[chosen_site]/(1.0*number_of_sites);
spin[chosen_site]=-spin[chosen_site];

}

//increment time
tau++;

}


void update_potential(void){

int protocol_time_index=(int) (1.0*tau/(1.0*protocol_update_time));
if((protocol_time_index > tau_last_protocol_update) && (tau <= trajectory_time)){

//reset counter
tau_last_protocol_update=protocol_time_index;

//update protocol
if(protocol_time_index>=number_of_protocol_timepoints){protocol_time_index=number_of_protocol_timepoints-1;}

internal_nrg+=cee_instantaneous[1]*magnetization*number_of_sites;
for(int i=0;i<number_of_control_parameters;i++){cee_instantaneous[i]=cee[protocol_time_index][i];}
internal_nrg-=cee_instantaneous[1]*magnetization*number_of_sites;

}}

void simulation_step(void){

update_potential();
glauber_step();
output_trajectory_data();

}

void run_trajectory(void){

double e1,e2;

set_initial_conditions();
e1=internal_nrg; //set above

while(tau<trajectory_time){simulation_step();}

//final-time potential
internal_nrg+=cee_instantaneous[1]*magnetization*number_of_sites;
for(int i=0;i<number_of_control_parameters;i++){cee_instantaneous[i]=cee_boundary[i][1];}
internal_nrg-=cee_instantaneous[1]*magnetization*number_of_sites;

//endpoint
e2=internal_nrg;
entprod+=(e2-e1)/(cee_instantaneous[0]);

}

void load_default_protocol(void){

//protocol for 0^+ < t < tf^- (boundary conditions for t=0 and t=tf enforced elsewhere)

int i,j;
double dummy1;

snprintf(st, sizeof(st), "default_protocol.dat");
ifstream in;

in.open((st));
if(!in){
cerr << "Error: could not open " << (st) << endl;
exit(1);
}

for(i=0;i<number_of_protocol_timepoints;i++){
for(j=0;j<number_of_control_parameters;j++){

in >> dummy1;

//temperature safeguard
if((j==0) && (dummy1<1e-3)){dummy1=1e-3;}

cee[i][j]=dummy1;

if(in.fail()){
cerr << "Error reading data from " << st << " at line " << i << ", column " << j << endl;
exit(1);

}}}

in.close();

}


void load_protocol(void){

//protocol for 0^+ < t < tf^- (boundary conditions for t=0 and t=tf enforced elsewhere)

int i,j;
double dummy1;

snprintf(st, sizeof(st), "input_control_parameters.dat");
ifstream in;

in.open((st));
if(!in){
cerr << "Error: could not open " << (st) << endl;
exit(1);
}

for(i=0;i<number_of_protocol_timepoints;i++){
for(j=0;j<number_of_control_parameters;j++){

in >> dummy1;

//temperature safeguard
if((j==0) && (dummy1<1e-3)){dummy1=1e-3;}

cee[i][j]=dummy1;

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
double op=0.0;

//initialize lattice
initialize_lattice();

for(i=0;i<n_traj;i++){

run_trajectory();
op+=entprod/(1.0*n_traj);

}

return (op);

}


void final_answer(void){

snprintf(st,sizeof(st),"report_answer.dat");
ofstream out1(st,ios::out);

int i;
int n_traj=1000;
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

//for ep histogram
double *entprod_values = (double *)malloc(n_traj_pictures * sizeof(double));
if(entprod_values == NULL){fprintf(stderr, "Error: could not allocate memory for entprod_values\n");exit(1);}

mean_entprod=0.0;
mean_magnetization=0.0;

q_record_trajectory=1;
traj_number=0;
for(i=0;i<n_traj_pictures;i++){

run_trajectory();

//for histogram
entprod_values[i]=entprod;

//for averages
mean_entprod+=entprod/(1.0*n_traj_pictures);
mean_magnetization+=magnetization/(1.0*n_traj_pictures);

traj_number++;

}
q_record_trajectory=0;

//magnetization histograms
//for(i=0;i<number_of_pictures;i++){output_histogram_magnetization(i);}

//entprod histogram
make_entprod_histogram(entprod_values,n_traj_pictures);
free(entprod_values);

//mean values
cout << " <sigma>= " << mean_entprod << ", <m> = " << mean_magnetization << endl;

//run python movie and picture maker
snprintf(st,sizeof(st),"python movie.py");
cout << st << endl;
system(st);


}

void output_trajectory_data(void){
int picture_index=(int) (1.0*tau/(1.0*picture_update_time));
if((q_record_trajectory==1) && (picture_index>tau_last_picture_update)){

//reset counter
tau_last_picture_update=picture_index;
if(picture_index>=number_of_pictures){picture_index=number_of_pictures-1;}

magnetization_time[traj_number][picture_index]=magnetization;

//lattice
output_lattice(picture_index);

}}


void output_histogram_magnetization(int time_slice){

snprintf(st,sizeof(st),"report_magnetization_time_%d.dat",time_slice);
ofstream out1(st,ios::out);

int i;
int nbin=0;
const int bins=50;

double histo[bins];
double histo_width;
double maxpos=0.0;
double minpos=0.0;

//magnetization recorded elswehere
for(i=0;i<bins;i++){histo[i]=0.0;}

for(i=0;i<n_traj_pictures;i++){

if(i==0){maxpos=magnetization_time[i][time_slice];minpos=magnetization_time[i][time_slice];}
else{

if(magnetization_time[i][time_slice]>maxpos){maxpos=magnetization_time[i][time_slice];}
if(magnetization_time[i][time_slice]<minpos){minpos=magnetization_time[i][time_slice];}
}}

//record width
histo_width=(maxpos-minpos)/(1.0*bins);

//safeguard
if(histo_width<1e-6){histo_width=1e-6;}

for(i=0;i<n_traj_pictures;i++){

nbin=(int) (1.0*bins*(magnetization_time[i][time_slice]-minpos)/(maxpos-minpos));
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



void output_protocol(void){

int i,j;
double t1;
ofstream out[number_of_control_parameters],out_para;

for(j=0;j<number_of_control_parameters;j++){
snprintf(st,sizeof(st),"report_control_parameters_%d.dat",j);
out[j].open(st,ios::out);
}

snprintf(st,sizeof(st),"report_control_parameters_parameteric.dat");
out_para.open(st,ios::out);

//initial values
for(j=0;j<number_of_control_parameters;j++){out[j] << 0 << " " << cee_boundary[j][0] << endl;}
out_para << cee_boundary[0][0] << " " << cee_boundary[1][0] << endl;

for(i=0;i<number_of_protocol_timepoints;i++){

t1=i/(number_of_protocol_timepoints-1.0);
for(j=0;j<number_of_control_parameters;j++){out[j] << t1 << " " << cee[i][j] << endl;}
out_para << cee[i][0] << " " << cee[i][1] << endl;


}

//final values
for(j=0;j<number_of_control_parameters;j++){out[j] << 1 << " " << cee_boundary[j][1] << endl;}
out_para << cee_boundary[0][1] << " " << cee_boundary[1][1] << endl;


//close output files
for(j=0;j<number_of_control_parameters;j++){out[j].close();}
out_para.close();

}

void make_entprod_histogram(double *entprod_values,int n_traj) {

    const int bins = 100;
    double histo[bins] = {0.0};
    double entprod_min = entprod_values[0];
    double entprod_max = entprod_values[0];

    // find min and max
    for (int i = 1; i < n_traj; i++) {
        if (entprod_values[i] < entprod_min) entprod_min = entprod_values[i];
        if (entprod_values[i] > entprod_max) entprod_max = entprod_values[i];
    }

    double bin_width = (entprod_max - entprod_min) / bins;
    if (bin_width < 1e-6) bin_width = 1e-6; // safeguard

    // fill histogram
    for (int i = 0; i < n_traj; i++) {
        int bin = (int)((entprod_values[i] - entprod_min) / (entprod_max - entprod_min) * bins);
        if (bin == bins) bin--; // include upper edge
        histo[bin] += 1.0 / (n_traj * bin_width);
    }

    // output
    snprintf(st, sizeof(st), "report_entprod_histogram.dat");
    ofstream out(st, ios::out);
    for (int i = 0; i < bins; i++) {
        double center = entprod_max * i / (double)bins + entprod_min * (1.0 - i / (double)bins) + 0.5 * (entprod_max - entprod_min) / bins;
        if (histo[i] > 0) {
            out << center << " " << histo[i] << endl;
        }
    }
    out.close();
}



int node_x_coord(int i){

return (i % lx);

}

int node_y_coord(int i){

return (i / lx);

}

int choose_spin(int t1){

return spin_sequence[t1 % number_of_sites];
// return spin_sequence[(int)(uniform_dist(rng_engine)*number_of_sites)];

}


void initialize_lattice(void){

int i;
int x1,y1;

for(i=0;i<number_of_sites;i++){

x1=node_x_coord(i);
y1=node_y_coord(i);

//left
if(x1==0){neighbors[i][0]=i+(lx-1);}
else{neighbors[i][0]=i-1;}

//right
if(x1==lx-1){neighbors[i][1]=i-(lx-1);}
else{neighbors[i][1]=i+1;}

//down
if(y1==0){neighbors[i][2]=i+lx*(lx-1);}
else{neighbors[i][2]=i-lx;}

//up
if(y1==lx-1){neighbors[i][3]=i-lx*(lx-1);}
else{neighbors[i][3]=i+lx;}

}
  
  
  /*
for(i=0;i<number_of_sites;i++){
cout << "node " << i << " (x,y) " << node_x_coord(i) << " " << node_y_coord(i) << endl;
for(int j=0;j<4;j++){
cout << " neighbor " << neighbors[i][j] << " (x,y) " << node_x_coord(neighbors[i][j]) << " " << node_y_coord(neighbors[i][j]) << endl;
}}
exit(2);
*/
  
// Populate spin_sequence with alternating row-parity pattern
int index = 0;
for (int phase = 0; phase < 2; ++phase) {
for (int row = 0; row < lx; ++row) {
int parity = (row % 2 == 0) ? 0 : 1;
if (phase % 2 == 1) parity = 1 - parity;
for (int col = 0; col < lx / 2; ++col) {
int x = (parity == 0) ? 2 * col : 2 * col + 1;
int y = row;
spin_sequence[index++] = y * lx + x;

}}}

}


double delta_nrg(int chosen_site){

int i;

double s2;
double de=0.0;
double s1=1.0*spin[chosen_site];
double magnetic_field = cee_instantaneous[1];

//coupling
for(i=0;i<4;i++){

s2=1.0*spin[neighbors[chosen_site][i]];
de+=2.0*s1*s2;

}

//field
de+=2.0*magnetic_field*s1;

return (de);

}

void output_lattice(int picture_index){
if(traj_number==0){

snprintf(st,sizeof(st),"report_lattice_time_%d.dat",picture_index);
ofstream out1(st,ios::out);

for(int y=lx-1;y>=0;y--){
for(int x=0;x<lx;x++){
int index=y*lx+x;
out1 << spin[index] << " ";
}
out1 << endl;
}

}}
