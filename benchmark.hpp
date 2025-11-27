#ifndef DYNEARTHSOL3D_BENCHMARK_HPP
#define DYNEARTHSOL3D_BENCHMARK_HPP

#include <chrono>
#include <string>
#include <map>
#include <vector>
#include <stack>
#include <iostream>
#include <fstream>

// High-resolution timer for benchmarking
class Timer {
private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    
    TimePoint start_time;
    bool running;
    
public:
    Timer() : running(false) {}
    
    void start() {
        start_time = Clock::now();
        running = true;
    }
    
    double elapsed() const {
        if (!running) return 0.0;
        auto end_time = Clock::now();
        std::chrono::duration<double> diff = end_time - start_time;
        return diff.count();
    }
    
    void stop() {
        running = false;
    }
};

// Statistics for a timed section
struct TimingStats {
    std::string name;
    double total_time;
    double min_time;
    double max_time;
    int count;
    
    TimingStats() 
        : name(""), total_time(0.0), min_time(1e100), max_time(0.0), count(0) {}
    
    TimingStats(const std::string& n) 
        : name(n), total_time(0.0), min_time(1e100), max_time(0.0), count(0) {}
    
    void add(double time) {
        total_time += time;
        min_time = std::min(min_time, time);
        max_time = std::max(max_time, time);
        count++;
    }
    
    double average() const {
        return count > 0 ? total_time / count : 0.0;
    }
};

// Benchmark manager for collecting timing data
class BenchmarkManager {
private:
    std::map<std::string, TimingStats> stats;
    std::stack<std::pair<std::string, std::chrono::time_point<std::chrono::high_resolution_clock>>> active_sections;
    bool enabled;
    
public:
    BenchmarkManager() : enabled(true) {}
    
    void enable() { enabled = true; }
    void disable() { enabled = false; }
    bool is_enabled() const { return enabled; }
    
    void start_section(const std::string& name) {
        if (!enabled) return;
        active_sections.push({name, std::chrono::high_resolution_clock::now()});
    }
    
    void end_section() {
        if (!enabled || active_sections.empty()) return;
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto start_pair = active_sections.top();
        active_sections.pop();
        
        std::chrono::duration<double> diff = end_time - start_pair.second;
        double elapsed = diff.count();
        
        // Add to stats
        if (stats.find(start_pair.first) == stats.end()) {
            stats[start_pair.first] = TimingStats(start_pair.first);
        }
        stats[start_pair.first].add(elapsed);
    }
    
    void print_summary(std::ostream& os = std::cout) const {
        if (stats.empty()) {
            os << "No timing data collected.\n";
            return;
        }
        
        os << "\n=== Performance Summary ===\n";
        os << "Section                          Count      Total(s)    Avg(s)     Min(s)     Max(s)\n";
        os << "--------------------------------------------------------------------------------\n";
        
        for (const auto& pair : stats) {
            const TimingStats& s = pair.second;
            char buffer[256];
            snprintf(buffer, sizeof(buffer), 
                    "%-30s %6d %12.4f %10.4f %10.4f %10.4f\n",
                    s.name.c_str(), s.count, s.total_time, 
                    s.average(), s.min_time, s.max_time);
            os << buffer;
        }
        os << "================================================================================\n";
    }
    
    void write_csv(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open " << filename << " for writing\n";
            return;
        }
        
        file << "Section,Count,Total_s,Average_s,Min_s,Max_s\n";
        for (const auto& pair : stats) {
            const TimingStats& s = pair.second;
            file << s.name << "," << s.count << "," << s.total_time << ","
                 << s.average() << "," << s.min_time << "," << s.max_time << "\n";
        }
        
        file.close();
        std::cout << "Timing data written to " << filename << "\n";
    }
    
    const std::map<std::string, TimingStats>& get_stats() const {
        return stats;
    }
};

// RAII helper for automatic section timing
class ScopedTimer {
private:
    BenchmarkManager& manager;
    std::string section_name;
    
public:
    ScopedTimer(BenchmarkManager& mgr, const std::string& name) 
        : manager(mgr), section_name(name) {
        manager.start_section(section_name);
    }
    
    ~ScopedTimer() {
        manager.end_section();
    }
};

// Global benchmark manager instance
extern BenchmarkManager g_benchmark;

// Convenience macros
#define BENCHMARK_SECTION(name) ScopedTimer _scoped_timer_##__LINE__(g_benchmark, name)
#define BENCHMARK_START(name) g_benchmark.start_section(name)
#define BENCHMARK_END() g_benchmark.end_section()

#endif // DYNEARTHSOL3D_BENCHMARK_HPP
