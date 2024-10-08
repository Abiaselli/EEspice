import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
import argparse
 
def run_executable(executable_path):
    subprocess.run(executable_path, shell=True)
 
def main():
    parser = argparse.ArgumentParser(description="Run an executable N times in parallel and measure total execution time.")
    parser.add_argument("executable", help="Path to the executable file")
    parser.add_argument("-n", type=int, default=1, help="Number of times to run the executable (default: 1)")
    args = parser.parse_args()
 
    executable_path = args.executable
    n = args.n
 
    start_time = time.time()
 
    with ThreadPoolExecutor(max_workers=n) as executor:
        futures = [executor.submit(run_executable, executable_path) for _ in range(n)]
        for future in futures:
            future.result()
 
    end_time = time.time()
    total_time = end_time - start_time
 
    print(f"Executed {executable_path} {n} times in parallel")
    print(f"Total execution time: {total_time:.2f} seconds")
 
if __name__ == "__main__":
    main()