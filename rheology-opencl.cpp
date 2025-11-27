#ifdef WITH_OPENCL
#include <CL/cl.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstring>

#include "rheology-opencl.hpp"
#include "parameters.hpp"

// Global instance that can be used throughout the application
StressUpdateGPU* gpu_stress_update = nullptr;

bool initializeGPUAcceleration(const std::string& device_type) {
    if (gpu_stress_update == nullptr) {
        try {
            gpu_stress_update = new StressUpdateGPU(device_type);
            if (!gpu_stress_update->isInitialized()) {
                delete gpu_stress_update;
                gpu_stress_update = nullptr;
                return false;
            }
            std::cout << "GPU acceleration initialized: " << gpu_stress_update->getDeviceInfo() << std::endl;
            return true;
        }
        catch (const std::exception& e) {
            std::cerr << "Failed to initialize GPU acceleration: " << e.what() << std::endl;
            delete gpu_stress_update;
            gpu_stress_update = nullptr;
            return false;
        }
    }
    return gpu_stress_update->isInitialized();
}

void finalizeGPUAcceleration() {
    delete gpu_stress_update;
    gpu_stress_update = nullptr;
}

#ifdef WITH_OPENCL
// Implementation details using OpenCL
struct StressUpdateGPU::OpenCLImpl {
    cl_platform_id platform;
    cl_device_id device;
    cl_context context;
    cl_command_queue queue;
    cl_program program;
    cl_kernel kernel;
    
    // Buffers for element data
    cl_mem input_buffer;
    cl_mem output_buffer;
    
    size_t global_work_size;
    size_t local_work_size;
    
    // Constructor initializes to null
    OpenCLImpl() : platform(nullptr), device(nullptr), context(nullptr), 
                  queue(nullptr), program(nullptr), kernel(nullptr),
                  input_buffer(nullptr), output_buffer(nullptr),
                  global_work_size(0), local_work_size(0) {}
    
    // Clean up resources
    ~OpenCLImpl() {
        if (kernel) clReleaseKernel(kernel);
        if (program) clReleaseProgram(program);
        if (input_buffer) clReleaseMemObject(input_buffer);
        if (output_buffer) clReleaseMemObject(output_buffer);
        if (queue) clReleaseCommandQueue(queue);
        if (context) clReleaseContext(context);
    }
};

// Get string representation of OpenCL error codes
std::string StressUpdateGPU::getOpenCLErrorString(int error) const {
    switch (error) {
        case CL_SUCCESS: return "CL_SUCCESS";
        case CL_DEVICE_NOT_FOUND: return "CL_DEVICE_NOT_FOUND";
        case CL_DEVICE_NOT_AVAILABLE: return "CL_DEVICE_NOT_AVAILABLE";
        case CL_COMPILER_NOT_AVAILABLE: return "CL_COMPILER_NOT_AVAILABLE";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case CL_OUT_OF_RESOURCES: return "CL_OUT_OF_RESOURCES";
        case CL_OUT_OF_HOST_MEMORY: return "CL_OUT_OF_HOST_MEMORY";
        case CL_PROFILING_INFO_NOT_AVAILABLE: return "CL_PROFILING_INFO_NOT_AVAILABLE";
        case CL_MEM_COPY_OVERLAP: return "CL_MEM_COPY_OVERLAP";
        case CL_IMAGE_FORMAT_MISMATCH: return "CL_IMAGE_FORMAT_MISMATCH";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case CL_BUILD_PROGRAM_FAILURE: return "CL_BUILD_PROGRAM_FAILURE";
        case CL_MAP_FAILURE: return "CL_MAP_FAILURE";
        case CL_INVALID_VALUE: return "CL_INVALID_VALUE";
        case CL_INVALID_DEVICE_TYPE: return "CL_INVALID_DEVICE_TYPE";
        case CL_INVALID_PLATFORM: return "CL_INVALID_PLATFORM";
        case CL_INVALID_DEVICE: return "CL_INVALID_DEVICE";
        case CL_INVALID_CONTEXT: return "CL_INVALID_CONTEXT";
        case CL_INVALID_QUEUE_PROPERTIES: return "CL_INVALID_QUEUE_PROPERTIES";
        case CL_INVALID_COMMAND_QUEUE: return "CL_INVALID_COMMAND_QUEUE";
        case CL_INVALID_HOST_PTR: return "CL_INVALID_HOST_PTR";
        case CL_INVALID_MEM_OBJECT: return "CL_INVALID_MEM_OBJECT";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        case CL_INVALID_IMAGE_SIZE: return "CL_INVALID_IMAGE_SIZE";
        case CL_INVALID_SAMPLER: return "CL_INVALID_SAMPLER";
        case CL_INVALID_BINARY: return "CL_INVALID_BINARY";
        case CL_INVALID_BUILD_OPTIONS: return "CL_INVALID_BUILD_OPTIONS";
        case CL_INVALID_PROGRAM: return "CL_INVALID_PROGRAM";
        case CL_INVALID_PROGRAM_EXECUTABLE: return "CL_INVALID_PROGRAM_EXECUTABLE";
        case CL_INVALID_KERNEL_NAME: return "CL_INVALID_KERNEL_NAME";
        case CL_INVALID_KERNEL_DEFINITION: return "CL_INVALID_KERNEL_DEFINITION";
        case CL_INVALID_KERNEL: return "CL_INVALID_KERNEL";
        case CL_INVALID_ARG_INDEX: return "CL_INVALID_ARG_INDEX";
        case CL_INVALID_ARG_VALUE: return "CL_INVALID_ARG_VALUE";
        case CL_INVALID_ARG_SIZE: return "CL_INVALID_ARG_SIZE";
        case CL_INVALID_KERNEL_ARGS: return "CL_INVALID_KERNEL_ARGS";
        case CL_INVALID_WORK_DIMENSION: return "CL_INVALID_WORK_DIMENSION";
        case CL_INVALID_WORK_GROUP_SIZE: return "CL_INVALID_WORK_GROUP_SIZE";
        case CL_INVALID_WORK_ITEM_SIZE: return "CL_INVALID_WORK_ITEM_SIZE";
        case CL_INVALID_GLOBAL_OFFSET: return "CL_INVALID_GLOBAL_OFFSET";
        case CL_INVALID_EVENT_WAIT_LIST: return "CL_INVALID_EVENT_WAIT_LIST";
        case CL_INVALID_EVENT: return "CL_INVALID_EVENT";
        case CL_INVALID_OPERATION: return "CL_INVALID_OPERATION";
        case CL_INVALID_GL_OBJECT: return "CL_INVALID_GL_OBJECT";
        case CL_INVALID_BUFFER_SIZE: return "CL_INVALID_BUFFER_SIZE";
        case CL_INVALID_MIP_LEVEL: return "CL_INVALID_MIP_LEVEL";
        case CL_INVALID_GLOBAL_WORK_SIZE: return "CL_INVALID_GLOBAL_WORK_SIZE";
        default: return "UNKNOWN_ERROR";
    }
}

StressUpdateGPU::StressUpdateGPU(const std::string& device_type) 
    : impl(new OpenCLImpl()), initialized(false), using_gpu(false)
{
    cl_int err;
    cl_uint num_platforms;
    
    // Get platform
    err = clGetPlatformIDs(0, nullptr, &num_platforms);
    if (err != CL_SUCCESS || num_platforms == 0) {
        std::cerr << "No OpenCL platforms found." << std::endl;
        return;
    }
    
    std::vector<cl_platform_id> platforms(num_platforms);
    err = clGetPlatformIDs(num_platforms, platforms.data(), nullptr);
    if (err != CL_SUCCESS) {
        std::cerr << "Failed to get platform IDs: " << getOpenCLErrorString(err) << std::endl;
        return;
    }
    
    // Find first platform with a suitable device
    cl_device_type device_type_flag = CL_DEVICE_TYPE_ALL;
    if (device_type == "GPU") {
        device_type_flag = CL_DEVICE_TYPE_GPU;
    } else if (device_type == "CPU") {
        device_type_flag = CL_DEVICE_TYPE_CPU;
    }
    
    bool found_device = false;
    for (cl_platform_id platform : platforms) {
        cl_uint num_devices;
        err = clGetDeviceIDs(platform, device_type_flag, 0, nullptr, &num_devices);
        if (err == CL_SUCCESS && num_devices > 0) {
            std::vector<cl_device_id> devices(num_devices);
            err = clGetDeviceIDs(platform, device_type_flag, num_devices, devices.data(), nullptr);
            if (err == CL_SUCCESS) {
                impl->platform = platform;
                impl->device = devices[0]; // Choose first device
                found_device = true;
                
                // Check if it's a GPU
                cl_device_type actual_type;
                clGetDeviceInfo(impl->device, CL_DEVICE_TYPE, sizeof(cl_device_type), &actual_type, nullptr);
                using_gpu = (actual_type & CL_DEVICE_TYPE_GPU) != 0;
                
                break;
            }
        }
    }
    
    if (!found_device) {
        std::cerr << "No suitable OpenCL device found for type: " << device_type << std::endl;
        return;
    }
    
    // Create context
    cl_context_properties context_props[] = {
        CL_CONTEXT_PLATFORM, reinterpret_cast<cl_context_properties>(impl->platform),
        0
    };
    
    impl->context = clCreateContext(context_props, 1, &impl->device, nullptr, nullptr, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Failed to create OpenCL context: " << getOpenCLErrorString(err) << std::endl;
        return;
    }
    
    // Create command queue
    impl->queue = clCreateCommandQueue(impl->context, impl->device, 0, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "Failed to create command queue: " << getOpenCLErrorString(err) << std::endl;
        return;
    }
    
    // Proceed with kernel initialization
    initialized = initializeKernel();
}

bool StressUpdateGPU::initializeKernel() {
    // Read OpenCL kernel source code here
    // We will implement this in a separate function later
    
    // For now, return false to indicate initialization failure
    // This will be replaced with actual kernel initialization
    return false;
}

StressUpdateGPU::~StressUpdateGPU() {
    delete impl;
}

bool StressUpdateGPU::isInitialized() const {
    return initialized;
}

bool StressUpdateGPU::isUsingGPU() const {
    return using_gpu;
}

std::string StressUpdateGPU::getDeviceInfo() const {
    if (!initialized) {
        return "Not initialized";
    }
    
    cl_char vendor[1024] = {0};
    cl_char name[1024] = {0};
    
    clGetDeviceInfo(impl->device, CL_DEVICE_VENDOR, sizeof(vendor), vendor, nullptr);
    clGetDeviceInfo(impl->device, CL_DEVICE_NAME, sizeof(name), name, nullptr);
    
    std::stringstream ss;
    ss << vendor << " " << name << " (" << (using_gpu ? "GPU" : "CPU") << ")";
    return ss.str();
}

void StressUpdateGPU::processStressUpdate(
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
    double_vec& denergy)
{
    // This will be implemented later with actual OpenCL kernel execution
    std::cerr << "GPU stress update not fully implemented yet." << std::endl;
}

#else // !WITH_OPENCL

// Fallback implementation when OpenCL is not available
StressUpdateGPU::StressUpdateGPU(const std::string& device_type) : impl(nullptr), initialized(false), using_gpu(false) {
    std::cerr << "OpenCL support not compiled in. Build with -DWITH_OPENCL=ON to enable GPU acceleration." << std::endl;
}

StressUpdateGPU::~StressUpdateGPU() {}

bool StressUpdateGPU::isInitialized() const { return false; }

bool StressUpdateGPU::isUsingGPU() const { return false; }

std::string StressUpdateGPU::getDeviceInfo() const { return "OpenCL not available"; }

bool StressUpdateGPU::initializeKernel() { return false; }

void StressUpdateGPU::processStressUpdate(
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
    double_vec& denergy)
{
    // Do nothing, OpenCL not available
}

std::string StressUpdateGPU::getOpenCLErrorString(int error) const { 
    return "OpenCL not available"; 
}

#endif // WITH_OPENCL 