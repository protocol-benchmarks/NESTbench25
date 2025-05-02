
#include <random>
static std::mt19937_64 rng_engine; // default-constructed engine
static std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

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

using namespace std;
char st[256];

long int time_last_launch;
long int time_limit = 4*60*60*24;

const int n_traj = 50;
const int parents = 10;

int parent_list[parents];
int iteration = 0;
int generation = 0;

double phi_value[n_traj];

//other measured parameters
const int number_of_measured_parameters=1;
double measured_parameters[number_of_measured_parameters];

//functions void
void setup(void);
void relaunch(void);
void jobcomplete(int i);
void determine_parents(void);
void ssystem(const std::string& st);

//functions int
int q_jobs_done(void);

//functions bool
bool q_file_exist(const char *fileName);

int main(void){

//seed RN generator
rng_engine.seed(static_cast<uint64_t>(time(NULL)));

//setup
setup();

//GA
while(2>1){

if(q_jobs_done()==1){
determine_parents();
relaunch();
}

sleep(5);

}}


void setup(void){

int n;
 
//store data
snprintf(st,sizeof(st),"mkdir genome_and_phenotype");
ssystem(st);
//loop through parameter cases
for(n=0;n<n_traj;n++){

//directory
snprintf(st,sizeof(st),"mkdir iteration_%d",n);
ssystem(st);
//parameters for etna
snprintf(st,sizeof(st),"swarm_%d.sh",n);
ofstream soutput(st,ios::out);

soutput << "#!/usr/bin/env bash" << endl;
snprintf(st,sizeof(st),"#SBATCH --job-name=swarm_%d",n);
soutput << st << endl;
soutput << "#SBATCH --partition=lr7" << endl;
soutput << "#SBATCH --account=nano" << endl;
soutput << "#SBATCH --qos=lr_normal" << endl;
snprintf(st,sizeof(st),"#SBATCH --nodes=%d",1);
soutput << st << endl;
soutput << "#SBATCH --time=36:00:00" << endl;
//soutput << "#SBATCH --mem=1000mb" << endl;
soutput << "#SBATCH --ntasks=1" << endl;
snprintf(st,sizeof(st),"./swarm_%d",n);
soutput << st << endl;
soutput.close();
 
//copy executable
snprintf(st,sizeof(st),"mv swarm_%d.sh iteration_%d",n,n);
ssystem(st);

ofstream output("input_parameters.dat",ios::out);
output << n << " " << generation << endl;

snprintf(st,sizeof(st),"cp input_parameters.dat iteration_%d",n);
ssystem(st);

//copy code
snprintf(st,sizeof(st),"cp ga_process.c iteration_%d",n);
ssystem(st);

snprintf(st,sizeof(st),"cp *engine* iteration_%d",n);
ssystem(st);

//copy start net
snprintf(st,sizeof(st),"cp start_net.dat iteration_%d",n);
ssystem(st);

//copy file
jobcomplete(0);
snprintf(st,sizeof(st),"mv jobcomplete.dat iteration_%d",n);
ssystem(st);

//compile code
snprintf(st,sizeof(st),
    "cd iteration_%d; "
    "g++ -c engine_trap_overdamped.c -o engine_trap_overdamped.o; "
    "ar rcs libengine_trap_overdamped.a engine_trap_overdamped.o; "
    "g++ -Wall -mcmodel=medium -o swarm_%d ga_process.c -L. -lengine_trap_overdamped -lm -O",
    n, n);
ssystem(st);
cout << st << endl;
 
snprintf(st,sizeof(st),"cd iteration_%d; sbatch swarm_%d.sh",n,n); ssystem(st);
//pause
sleep(1);

}

time_last_launch=time(NULL);

}



void jobcomplete(int i){

 //snprintf(st,sizeof(st),"rm jobcomplete.dat");
 //ssystem(st);
 snprintf(st,sizeof(st),"jobcomplete.dat");
 ofstream output_job(st,ios::out);
 output_job << i << endl;
 output_job.close();

}


int q_jobs_done(void){

int i;
int q1=0;
int flag=0;
int counter=0;

for(i=0;i<n_traj;i++){

snprintf(st,sizeof(st),"iteration_%d/jobcomplete.dat",i);
ifstream infile(st, ios::in);
infile >> flag;

if(flag==1){counter++;}
flag=0;

}

cout << " gen " << generation << " counter " << counter << " time elapsed " << time(NULL)-time_last_launch << endl;

if(counter==n_traj){q1=1;}
else{q1=0;}

//safeguard
if((time(NULL)-time_last_launch>time_limit) && (counter>n_traj/2)){cout << " time limit exceeded " << endl;q1=1;}

return (q1);

}


bool q_file_exist(const char *fileName){
  
  std::ifstream infile(fileName);
  return infile.good();
 
 }

void determine_parents(void){

//int k;
int i,j;
int p1=0;
int flag;
int top_parent=-1;

double phi_min=1e5;
double phi_dummy=1e5;

//output files
snprintf(st,sizeof(st),"output_phi_min.dat");
ofstream output_phi_min(st,ios::app);


//read phi values
for(i=0;i<n_traj;i++){

snprintf(st,sizeof(st),"iteration_%d/jobcomplete.dat",i);
ifstream infile0(st, ios::in);
infile0 >> flag;

if(flag==1){

//input file 1
snprintf(st,sizeof(st),"iteration_%d/report_phi.dat",i);
if(q_file_exist(st)==1){
ifstream infile(st, ios::in);
while (!infile.eof()){infile >> phi_value[i];}}
else{phi_value[i]=phi_dummy; cout << " job " << i << " did not run " << endl;} //node failed

//input file 2
//snprintf(st,sizeof(st),"iteration_%d/report_order_parameters_gen_%d.dat",i,generation);
//if(q_file_exist(st)==1){
//ifstream infile2(st, ios::in);
//for(j=0;j<number_of_measured_parameters;j++){infile2 >> measured_parameters[j];}
//}

}}

//compute ordered list of parents
for(j=0;j<parents;j++){

for(i=0;i<n_traj;i++){

if(i==0){phi_min=phi_value[i];p1=i;}
else{
if(phi_value[i]<=phi_min){phi_min=phi_value[i];p1=i;}}}

//log and output
parent_list[j]=p1;
if(j==0){

top_parent=p1;
output_phi_min << generation << " " << phi_min << endl;

//for(k=0;k<number_of_measured_parameters;k++){

//snprintf(st,sizeof(st),"output_parameter_%d.dat",k);
//ofstream out1(st,ios::app);
//out1 << generation << " " << measured_parameters[k] << endl;
//out1.close();

//}

}

//flag
//output_phi_list << p1 << " " << phi_value[p1] << endl;
phi_value[p1]=phi_dummy;

}

//record data from top parent
snprintf(st,sizeof(st),"cp iteration_%d/report*.dat genome_and_phenotype/.",top_parent);
ssystem(st);
//network used to produce config
snprintf(st,sizeof(st),"cp iteration_%d/net_out_gen_%d.dat genome_and_phenotype/.",top_parent,generation);
ssystem(st);
//clean up
if(generation>1){
snprintf(st,sizeof(st),"rm genome_and_phenotype/net_out_gen_%d.dat",generation-2);
ssystem(st);
}

//increment generation counter
generation++;

}


void relaunch(void){

int i;
int parent_choice;

//loop through parameter cases
for(i=0;i<n_traj;i++){
  
if(i==0){parent_choice=parent_list[0];}
else{
parent_choice= (int) (drand48()*parents);
parent_choice=parent_list[parent_choice];
}
 
//input parameters file
ofstream output("input_parameters.dat",ios::out);
output << i << " " << generation << endl;

//copy file
snprintf(st,sizeof(st),"cp input_parameters.dat iteration_%d",i);
ssystem(st);

//copy net
snprintf(st,sizeof(st),"cp iteration_%d/net_out_gen_%d.dat iteration_%d/net_in_gen_%d.dat",parent_choice,generation-1,i,generation);
ssystem(st);
//cout << st << endl;

//clean up
if(generation>2){snprintf(st,sizeof(st),"cd iteration_%d; rm net_*_gen_%d.dat",i,generation-2);ssystem(st);}

//copy file
jobcomplete(0);
snprintf(st,sizeof(st),"mv jobcomplete.dat iteration_%d",i);
ssystem(st);
//clean up
//snprintf(st,sizeof(st),"cd iteration_%d; rm slurm*out",i); ssystem(st);
//launch
//snprintf(st,sizeof(st),"cd iteration_%d; sbatch swarm_%d.sh",i,i); ssystem(st);
//pause
//sleep(2);
 
}

time_last_launch=time(NULL);

}

void ssystem(const std::string& st) {

   int status=system(st.c_str());
   
   if(status!=0){cout << " error: " << st << " status " << status << endl;}

}
