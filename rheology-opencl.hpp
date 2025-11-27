#ifndef DYNEARTHSOL3D_RHEOLOGY_OPENCL_HPP
#define DYNEARTHSOL3D_RHEOLOGY_OPENCL_HPP

#include <string>
#include <vector>
#include "constants.hpp"
#include "parameters.hpp"

// Forward declaration of Variables struct
struct Variables;

// OpenCL context management for stress update
class StressUpdateGPU {
public:
    // Initialize OpenCL with the given device type (GPU, CPU, or ALL)
    StressUpdateGPU(const std::string& device_type = "ALL");
    
    // Cleanup OpenCL resources
    ~StressUpdateGPU();
    
    // Check if OpenCL initialization was successful
    bool isInitialized() const;
    
    // Returns true if running on GPU, false for CPU fallback
    bool isUsingGPU() const;
    
    // Get device name and info
    std::string getDeviceInfo() const;
    
    // Process stress update for all elements using OpenCL
    void processStressUpdate(
        const Variables& var,
        tensor_t& stress, 
        double_vec& stressyy, 
        double_vec& thermal_stress, 
        double_vec& dP,
        tensor_t& strain, 
        tensor_t& elastic_strain, 
        double_vec& plstrain, 
        double_vec& delta_plstrain, 
        double_vec& dtemp,
        tensor_t& strain_rate, 
        double_vec& power, 
        double_vec& tenergy, 
        double_vec& venergy, 
        double_vec& denergy
    );

private:
    // OpenCL handles - implemented in rheology-opencl.cpp
    struct OpenCLImpl;
    OpenCLImpl* impl;
    
    // Whether initialization was successful
    bool initialized;
    
    // Whether we're using GPU or CPU
    bool using_gpu;
    
    // Compile and initialize OpenCL kernel
    bool initializeKernel();
    
    // Helper to convert error codes to string
    std::string getOpenCLErrorString(int error) const;
};

// Global instance that can be used throughout the application
extern StressUpdateGPU* gpu_stress_update;

// Initialize GPU acceleration if supported
bool initializeGPUAcceleration(const std::string& device_type = "ALL");

// Cleanup GPU resources
void finalizeGPUAcceleration();

#endif // DYNEARTHSOL3D_RHEOLOGY_OPENCL_HPP 