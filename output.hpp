#ifndef DYNEARTHSOL3D_OUTPUT_HPP
#define DYNEARTHSOL3D_OUTPUT_HPP

#include "array2d.hpp"

class Output
{
private:
    const std::string &modelname;
    const double start_time;
    const int average_interval;
    const bool has_marker_output;
    int frame;

    // stuffs for averging fields
    double time0;
    array_t coord0;
    tensor_t strain0;
    tensor_t stress_avg;
    double_vec delta_plstrain_avg;

    void write_info(const Variables& var, double dt);

public:
    Output(const Param& param, double start_time, int start_frame);
    ~Output();
    void write(const Variables& var, bool is_averaged=true);
    void write_checkpoint(const Param& param, const Variables& var);
    void average_fields(Variables& var);

};


#endif
