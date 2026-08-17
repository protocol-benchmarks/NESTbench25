#ifndef ENGINE_ISING_H
#define ENGINE_ISING_H

#ifdef __cplusplus
extern "C" {
#endif

//functions exposed from engine_ising.c.

void final_answer(void); //100 instances of the order parameter calculated using 10^4 trajectories
void load_protocol(void);//reads input_control_parameters.dat
void output_protocol(void);
void output_optimal_protocol(void);
void visualize_protocol(void); //visualize the effect of the current protocol
void load_default_protocol(void); //loads default protocol
double calculate_order_parameter(int n_traj); //computes order parameter over n_traj trajectories


//run-time configuration (added in response to referee report)
void load_protocol_from(const char *filename); //reads protocol from a user-specified file
void set_trajectory_time(double tf); //set trajectory time; derived intervals are recomputed
void set_boundary_conditions(const char *spec); //comma-separated list c0_i,c0_f,c1_i,c1_f,...
void final_answer_n(long n_traj,long n_samples); //final_answer with user-specified sampling (0 = engine default)

#ifdef __cplusplus
}
#endif

#endif // ENGINE_ISING_H
