import cProfile
import pstats
import tracemalloc
import sys
import time
import json
from datetime import datetime
from pathlib import Path
from memory_profiler import profile
from functools import wraps

def create_report_directory():
    # Create a timestamped directory for this profiling session
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_dir = Path(f'profiling_results_{timestamp}')
    report_dir.mkdir(exist_ok=True)
    return report_dir

def write_memory_trace(snapshot, report_dir):
    # Write detailed memory allocation information
    memory_file = report_dir / 'memory_trace.txt'
    with open(memory_file, 'w') as f:
        stats = snapshot.statistics('lineno')
        f.write("=== Memory Allocation Trace ===\n\n")
        for stat in stats[:50]:  # Top 50 memory allocations
            f.write(f"{stat}\n")
            
def write_profile_stats(stats, report_dir):
    # Write detailed profiling statistics
    stats_file = report_dir / 'profile_stats.prof'
    stats.dump_stats(str(stats_file))
    
    # Write human-readable profile stats
    readable_stats = report_dir / 'profile_stats.txt'
    with open(readable_stats, 'w') as f:
        stats.stream = f  # Redirect output to file
        f.write("=== Function Timing Profile ===\n\n")
        stats.print_stats()
        f.write("\n=== Callers Information ===\n\n")
        stats.print_callers()

def write_summary_json(summary_data, report_dir):
    # Write summary metrics in JSON format for easy parsing
    summary_file = report_dir / 'summary_metrics.json'
    with open(summary_file, 'w') as f:
        json.dump(summary_data, f, indent=4)

def comprehensive_profile(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        # Create directory for this profiling session
        report_dir = create_report_directory()
        
        # Start memory tracking
        tracemalloc.start()
        
        # Profile execution time
        profiler = cProfile.Profile()
        profiler.enable()
        
        # Track initial memory state
        initial_mem = tracemalloc.get_traced_memory()
        start_time = time.time()
        
        try:
            # Execute the function
            result = func(*args, **kwargs)
            
            # Collect execution metrics
            execution_time = time.time() - start_time
            current_mem, peak_mem = tracemalloc.get_traced_memory()
            
            # Take memory snapshot
            snapshot = tracemalloc.take_snapshot()
            tracemalloc.stop()
            
            # Disable profiler and collect stats
            profiler.disable()
            stats = pstats.Stats(profiler).sort_stats('cumulative')
            
            # Prepare summary data
            summary_data = {
                'function_name': func.__name__,
                'timestamp': datetime.now().isoformat(),
                'execution_time_seconds': execution_time,
                'memory_metrics': {
                    'initial_memory_mb': initial_mem[0] / 1024 / 1024,
                    'final_memory_mb': current_mem / 1024 / 1024,
                    'peak_memory_mb': peak_mem / 1024 / 1024
                },
                'status': 'success'
            }
            
            # Write all results to files
            write_memory_trace(snapshot, report_dir)
            write_profile_stats(stats, report_dir)
            write_summary_json(summary_data, report_dir)
            
            # Write a human-readable summary
            with open(report_dir / 'summary.txt', 'w') as f:
                f.write(f"Profiling Summary for {func.__name__}\n")
                f.write(f"Timestamp: {summary_data['timestamp']}\n")
                f.write(f"Total execution time: {execution_time:.2f} seconds\n")
                f.write(f"Peak memory usage: {peak_mem / 1024 / 1024:.1f} MB\n")
                f.write(f"\nDetailed results can be found in:\n")
                f.write(f"- Memory trace: memory_trace.txt\n")
                f.write(f"- Profile stats: profile_stats.txt\n")
                f.write(f"- Machine-readable metrics: summary_metrics.json\n")
            
            return result
            
        except Exception as e:
            # Log error information if something goes wrong
            error_data = {
                'function_name': func.__name__,
                'timestamp': datetime.now().isoformat(),
                'error_type': type(e).__name__,
                'error_message': str(e),
                'status': 'error'
            }
            
            with open(report_dir / 'error_log.txt', 'w') as f:
                f.write(f"Error during profiling: {str(e)}\n")
            write_summary_json(error_data, report_dir)
            raise
            
    return wrapper

