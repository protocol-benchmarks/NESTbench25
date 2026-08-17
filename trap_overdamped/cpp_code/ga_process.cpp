
//import simulation engine
#include "engine_trap_overdamped.h"

#include <random>
static std::mt19937_64 rng_engine; // default-constructed engine
static std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
static std::normal_distribution<double> normal_dist;  // global object

//compile with
//g++ ga_process.c -L. -lengine_trap_overdamped -o sim

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

char st[256];

//mmm, pi
//const double pi=4.0*atan(1.0);

//production
const int q_production_run=1; //flag

//pre-trained?
const int q_read_start_net=0; //flag

//GA variables
int iteration;
int generation;
const int number_of_trajectories=10000;
double phi[number_of_trajectories];
double sigma_mutate=0.02;
double mean_phi;

//lookup table
const int number_of_input_parameters=1;
const int number_of_control_parameters=1;
const int number_of_protocol_timepoints=1000;

//neural-net parameters
const int depth=4; //stem of mushroom net, >=1
const int width=4;
const int width_final=10;

double protocol_offsets[number_of_control_parameters]={0};

//specific net parameters
const int number_of_nets=1;
const int number_inputs[number_of_nets]={number_of_input_parameters};
const int number_outputs[number_of_nets]={number_of_control_parameters};

//net parameter count
const int net_params0=number_inputs[0]*width+width*width*(depth-1)+width*width_final+depth*width+width_final+(width_final+1)*number_outputs[0];
const int number_of_net_params[number_of_nets]={net_params0};
const int number_of_net_parameters=number_of_net_params[0];

double inputs[number_of_input_parameters];
double outputs[number_of_control_parameters];
double* net_parameters = (double*) malloc(number_of_net_parameters * sizeof(double));


//functions void
void ga(void);
void read_net(void);
void initialize(void);
void output_net(void);
void mutate_net(void);
void jobcomplete(int i);
void initialize_net(void);
void run_net(int net_id);
void make_lookup_table(void);

//functions double
double gauss_rv(double sigma);

int main(void){

//simulation loop
int q_loop=1;

while(q_loop==1){

initialize();
ga();

if(q_production_run==0){q_loop=0;}

}

return 0;

}


void initialize(void){

//clean up
snprintf(st,sizeof(st),"rm report_*");
cout << st << endl;
cout << endl;
system(st);

ifstream infile0("input_parameters.dat", ios::in);
while (!infile0.eof ()){infile0 >> iteration >> generation;}

if(generation>1){
snprintf(st,sizeof(st),"rm net_*_gen_%d.dat",generation-2);
cout << st << endl;
cout << endl;
system(st);
}

//seed RN generator
rng_engine.seed(static_cast<uint64_t>(time(NULL)));

//initialize net
if(generation==0){

if(q_read_start_net==0){initialize_net();}
else{

snprintf(st,sizeof(st),"cp start_net.dat net_in_gen_0.dat");
system(st);

read_net();

}}
else{read_net();}

}

void read_net(void){

int i;

snprintf(st,sizeof(st),"net_in_gen_%d.dat",generation);
ifstream infile(st, ios::in);

for(i=0;i<number_of_net_parameters;i++){infile >> net_parameters[i];}

}

void output_net(void){

int i;

//parameter file
snprintf(st,sizeof(st),"net_out_gen_%d.dat",generation);
ofstream out_net(st,ios::out);

for(i=0;i<number_of_net_parameters;i++){out_net << net_parameters[i] << " ";}

}

void mutate_net(void){

int i;
for(i=0;i<number_of_net_parameters;i++){net_parameters[i]+=gauss_rv(sigma_mutate);}

}


void initialize_net(void){

int i;
for(i=0;i<number_of_net_parameters;i++){net_parameters[i]=gauss_rv(0.0);}

}


double gauss_rv(double sigma){

    normal_dist = std::normal_distribution<double>(0.0, sigma);
    return normal_dist(rng_engine);
    
}


void ga(void){

double op;

//mutate net
if((q_production_run==1) && (iteration!=0)){mutate_net();}

//make protocol lookup table
make_lookup_table();

//call engine
load_protocol();

//visual check
output_protocol();

//calculate order parameter
op=calculate_order_parameter(10000);

//visualize protocol
//visualize_protocol();

//answer
//final_answer();

//output evolutionary order parameter
snprintf(st,sizeof(st),"report_phi.dat");
ofstream output_phi(st,ios::out);
output_phi << op << endl;

//output net
output_net();

//jobcomplete
jobcomplete(1);

if (q_production_run == 1) {
    // Monitor jobcomplete.dat every 2 seconds.
    int flag = 1;  // Initialize flag to a non-zero value to enter the loop.
    while (flag != 0) {
        sleep(2);  // Pause for 2 seconds.
        std::ifstream jobFile("jobcomplete.dat");
        if (jobFile && (jobFile >> flag)) {
            // Successfully read flag value from the file.
    }}}


}

void jobcomplete(int i){

 //snprintf(st,sizeof(st),"rm jobcomplete.dat");
 //system(st);
 
 snprintf(st,sizeof(st),"jobcomplete.dat");
 ofstream output_job(st,ios::out);
 output_job << i << endl;
 output_job.close();

}

void make_lookup_table(void){

int i,j;

double t1;


snprintf(st,sizeof(st),"input_control_parameters.dat");
ofstream out1(st,ios::out);

for(i=0;i<number_of_protocol_timepoints;i++){

t1=i/(number_of_protocol_timepoints-1.0);

inputs[0]=t1;
run_net(0);

//out1 << t1 << " ";
for(j=0;j<number_of_control_parameters;j++){out1 << outputs[j]+protocol_offsets[j] << " ";}
out1 << endl;

}}

void run_net(int net_id){

//inputs fixed externally and outputs used externally

int i,j,k;

//param offset
int pid=0;

double mu=0.0;
double sigma=0.0;
double delta=1e-4;

double hidden_node[width][depth];
double hidden_node_final[width_final];
const int number_of_inputs=number_inputs[net_id];
const int number_of_outputs=number_outputs[net_id];

//surface layer
for(i=0;i<width;i++){
hidden_node[i][0]=net_parameters[pid];pid++;
for(j=0;j<number_of_inputs;j++){hidden_node[i][0]+=inputs[j]*net_parameters[pid];pid++;}
}

//layer norm (pre-activation)
mu=0.0;sigma=0.0;
for(i=0;i<width;i++){mu+=hidden_node[i][0];sigma+=hidden_node[i][0]*hidden_node[i][0];}
mu=mu/width;sigma=sigma/width;
sigma=sqrt(sigma-mu*mu)+delta;
for(i=0;i<width;i++){hidden_node[i][0]=(hidden_node[i][0]-mu)/sigma;}


//activation
for(i=0;i<width;i++){hidden_node[i][0]=tanh(hidden_node[i][0]);}

//stem layers
for(j=1;j<depth;j++){
for(i=0;i<width;i++){
hidden_node[i][j]=net_parameters[pid];pid++;
for(k=0;k<width;k++){hidden_node[i][j]+=hidden_node[k][j-1]*net_parameters[pid];pid++;}
}

//layer norm (pre-activation)
mu=0.0;sigma=0.0;
for(i=0;i<width;i++){mu+=hidden_node[i][j];sigma+=hidden_node[i][j]*hidden_node[i][j];}
mu=mu/width;sigma=sigma/width;
sigma=sqrt(sigma-mu*mu)+delta;
for(i=0;i<width;i++){hidden_node[i][j]=(hidden_node[i][j]-mu)/sigma;}


//activation
for(i=0;i<width;i++){hidden_node[i][j]=tanh(hidden_node[i][j]);}
}

//final layer
for(i=0;i<width_final;i++){
hidden_node_final[i]=net_parameters[pid];pid++;
for(j=0;j<width;j++){hidden_node_final[i]+=hidden_node[j][depth-1]*net_parameters[pid];pid++;}
}

//layer norm (pre-activation)
mu=0.0;sigma=0.0;
for(i=0;i<width_final;i++){mu+=hidden_node_final[i];sigma+=hidden_node_final[i]*hidden_node_final[i];}
mu=mu/width_final;sigma=sigma/width_final;
sigma=sqrt(sigma-mu*mu)+delta;
for(i=0;i<width_final;i++){hidden_node_final[i]=(hidden_node_final[i]-mu)/sigma;}


//activation
for(i=0;i<width_final;i++){hidden_node_final[i]=tanh(hidden_node_final[i]);}

//outputs
for(i=0;i<number_of_outputs;i++){
outputs[i]=net_parameters[pid];pid++;
for(j=0;j<width_final;j++){outputs[i]+=hidden_node_final[j]*net_parameters[pid];pid++;}
}

}
